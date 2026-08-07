#!/usr/bin/env python
"""limma-style moderated-t confirmation of a condition-linked canonical DECREASE
for one gene in one PRISM dataset, with a batch covariate.

Reuses the miner's sample<->condition<->batch mapping (build_sample_conditions) so
the grouping is identical to what `mine` reported. Empirical-Bayes variance
shrinkage (Smyth 2004) is fit across ALL peptides; the target gene's canonical
peptides are then read out. Prints the batch x condition balance table first, per
the project methodology rule (a disease contrast confounded with batch is not
interpretable even after correction).

Contrast: (condition == target_level)  vs  (all other levels), adjusting for the
covariate (batch/Plate) as fixed effects.

CLI:  python limma_confirm.py --dataset D --gene G --condition C --target L [--covariate Plate]
"""
import argparse, os, re, sys
from pathlib import Path
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from scipy import special, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cryptic_remine import (load_registry, build_sample_conditions,
                            canonical_peptides_for_gene, peptide_key)


# ---------- eBayes (shared with the iNeuron script) ----------
def trigamma_inverse(x):
    x = np.asarray(x, float); out = np.empty_like(x)
    for i, xi in np.ndenumerate(x):
        if not np.isfinite(xi): out[i] = np.nan; continue
        if xi > 1e7: out[i] = 1.0/np.sqrt(xi); continue
        if xi < 1e-6: out[i] = 1.0/xi; continue
        y = 0.5 + 1.0/xi
        for _ in range(50):
            tri = special.polygamma(1, y)
            d = tri*(1-tri/xi)/special.polygamma(2, y); y += d
            if -d/y < 1e-8: break
        out[i] = y
    return out


def fit_fdist(s2, df1):
    s2 = np.asarray(s2, float); ok = np.isfinite(s2) & (s2 > 0)
    z = np.log(s2[ok]); e = z - special.digamma(df1/2.0) + np.log(df1/2.0)
    emean = e.mean(); evar = e.var(ddof=1) * ok.sum()/(ok.sum()-1)
    adj = evar - special.polygamma(1, df1/2.0)
    if adj > 0:
        d0 = 2.0*float(trigamma_inverse(np.array([adj]))[0])
        s02 = np.exp(emean + special.digamma(d0/2.0) - np.log(d0/2.0))
    else:
        d0, s02 = np.inf, np.exp(emean)
    return d0, s02


def bh(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank; q[idx] = prev = min(prev, p[idx]*n/i)
    return q


def run_confirm(db_dir, dataset, gene, condition, target, covariate=None,
                matrix="peptides_log2_internal.parquet", outdir=None, verbose=True):
    reg = load_registry(db_dir)
    ds = next((d for d in reg["datasets"] if d["name"] == dataset), None)
    if ds is None:
        raise SystemExit(f"dataset {dataset} not in registry")
    prism = Path(ds["prism_dir"])
    frame = build_sample_conditions(prism, ds)          # indexed by matrix column
    if condition not in frame.columns:
        raise SystemExit(f"condition '{condition}' not found; have {list(frame.columns)}")

    sub = frame[frame[condition].notna()].copy()
    sub["is_target"] = (sub[condition].astype(str) == str(target)).astype(int)
    if sub["is_target"].sum() == 0:
        raise SystemExit(f"no samples with {condition}=={target}")

    cov = None
    if covariate and covariate in sub.columns and sub[covariate].nunique(dropna=True) > 1:
        cov = covariate
        bal = pd.crosstab(sub[condition].astype(str), sub[cov].astype(str))
        if verbose:
            print(f"\n[balance] {condition} x {cov} (watch for confounding):")
            print(bal.to_string())
    elif verbose:
        print(f"[note] no usable covariate ({covariate!r}) - unadjusted contrast")

    n_t = int(sub["is_target"].sum()); n_o = int((sub["is_target"] == 0).sum())
    if verbose:
        print(f"\n{gene} @ {dataset}: {condition}=={target}  (n_target={n_t} vs n_other={n_o})")

    # design matrix: intercept + covariate dummies + is_target (last column = effect)
    cols = [pd.Series(1.0, index=sub.index, name="intercept")]
    if cov:
        d = pd.get_dummies(sub[cov].astype(str), prefix=cov, drop_first=True).astype(float)
        cols += [d[c] for c in d.columns]
    cols.append(sub["is_target"].astype(float).rename("is_target"))
    X = pd.concat(cols, axis=1)
    Xv = X.values; p = Xv.shape[1]; j = p - 1        # index of is_target coef

    # load matrix (peptides x samples), keep our samples; fall back to the dataset's
    # declared matrix when the requested one is absent (e.g. Skyline-imported dirs)
    if not (prism / matrix).exists():
        matrix = ds.get("peptide_matrix", "corrected_peptides.parquet")
        if verbose:
            print(f"[note] requested matrix absent; using ds matrix '{matrix}'")
    pf = pq.ParquetFile(prism / matrix)
    key = peptide_key(pf.schema_arrow.names)
    scols = [c for c in sub.index if c in pf.schema_arrow.names]
    dfm = pf.read(columns=[key] + scols).to_pandas().set_index(key)
    Y = dfm[scols].apply(pd.to_numeric, errors="coerce")
    # log2 if the matrix looks linear (internal is already log2)
    if np.nanmedian(Y.values) > 1000:
        Y = np.log2(Y.where(Y > 0))
    Xv = X.loc[scols].values                          # align design to matrix col order
    good = np.isfinite(Y.values).all(1)
    Yg = Y.values[good]; peps = Y.index.values[good]
    n = Xv.shape[0]; dfres = n - p

    pinv = np.linalg.pinv(Xv)                          # p x n
    B = pinv @ Yg.T                                    # p x npep
    resid = Yg.T - Xv @ B
    rss = (resid**2).sum(0); s2 = rss/dfres
    XtX_inv = np.linalg.pinv(Xv.T @ Xv); c_jj = XtX_inv[j, j]
    logFC = B[j]                                       # target vs other (adjusted)
    d0, s02 = fit_fdist(s2, dfres)
    if np.isinf(d0):
        s2_post = np.full_like(s2, s02); df_tot = np.full_like(s2, 1e6)
    else:
        s2_post = (d0*s02 + dfres*s2)/(d0+dfres); df_tot = np.full_like(s2, dfres+d0)
    se = np.sqrt(s2_post * c_jj); t = logFC/se
    pval = 2.0*stats.t.sf(np.abs(t), df_tot); q = bh(pval)
    if verbose:
        modf = dfres + (0 if np.isinf(d0) else d0)
        print(f"peptides tested: {good.sum()}  residual df={dfres}  prior d0={d0:.2f}  moderated df={modf:.1f}")

    res = pd.DataFrame({"peptide": peps, "logFC_target_minus_other": logFC,
                        "mod_t": t, "mod_p": pval, "BH": q}).set_index("peptide")

    cps, uni, srctag = canonical_peptides_for_gene(prism, gene, ds)
    g = res[res.index.isin(cps)]
    if len(g) == 0:
        print(f"** {gene}: no canonical peptides in matrix (src={srctag})"); return None
    pv = np.sort(g["mod_p"].values); m = len(pv)
    simes = float(np.min(pv*m/np.arange(1, m+1)))
    ndown = int(((g["BH"] < 0.05) & (g["logFC_target_minus_other"] < 0)).sum())
    summary = dict(gene=gene, dataset=dataset, condition=condition, target=target,
                   covariate=cov, n_target=n_t, n_other=n_o, n_canon_pep=len(g),
                   n_sig_down_BH05=ndown, median_logFC=float(g["logFC_target_minus_other"].median()),
                   min_BH=float(g["BH"].min()), simes_gene_p=simes)
    if verbose:
        print(f"\n=== {gene}: median logFC={summary['median_logFC']:+.3f} "
              f"({'DOWN' if summary['median_logFC']<0 else 'UP'} in {target})  "
              f"Simes p={simes:.4g}  peptide-level sig-down(BH<.05)={ndown}/{len(g)}")
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        tag = f"{gene}__{dataset}__{condition}_{target}".replace(" ", "")
        g.reset_index().to_csv(os.path.join(outdir, tag + "_peptides.csv"), index=False)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db-dir", default=r"G:\Manuscripts\LatimerCryptic\cryptic_db")
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--covariate", default=None)
    ap.add_argument("--matrix", default="peptides_log2_internal.parquet")
    a = ap.parse_args()
    out = os.path.join(a.db_dir, "remine_runs", "limma")
    run_confirm(Path(a.db_dir), a.dataset, a.gene, a.condition, a.target,
                covariate=a.covariate, matrix=a.matrix, outdir=out)


if __name__ == "__main__":
    main()
