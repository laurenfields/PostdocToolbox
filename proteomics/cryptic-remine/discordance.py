#!/usr/bin/env python
"""Proteome-wide PeCorA-style peptide-discordance scan.

For every protein (>= MIN_PEP peptides), every peptide is tested for a
condition x peptide interaction (batch-adjusted): does this peptide's disease
response differ from its protein's other peptides? That within-protein
discordance is the stoichiometry footprint of a cryptic/proteoform event and is
detectable even when the cryptic peptide itself was never in the search space.

Dense (Skyline peak-boundary imputed) matrices make single-dataset hits
unreliable, so the payoff is the CROSS-DATASET reproducible table: the same
stripped peptide, same direction, discordant in >= 2 datasets.

CLI:
  # scan one dataset against a disease contrast
  python cryptic_remine.py discordance --dataset huad_smtg_cleandx \
      --condition "Cognitive Status" --target Dementia [--covariate batch]
  # after scanning several, build the reproducible cross-dataset table
  python cryptic_remine.py discordance --aggregate
"""
import os, re, glob
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import pyarrow.parquet as pq

from cryptic_remine import (load_registry, build_sample_conditions, peptide_key,
                            resolve_cols)

MIN_PEP, MIN_OBS, MIN_GRP, BH_SIG = 4, 8, 6, 0.10


def _strip_mod(p):
    return re.sub(r"[^A-Za-z]", "", re.sub(r"\[[^\]]*\]|\([^)]*\)", "", str(p))).upper()


def _gene_peptides(P, ds):
    """gene -> [peptides] via protein_groups.csv, else transition-level merged_data."""
    pg = P / ds.get("protein_groups", "protein_groups.csv")
    if pg.exists():
        d = pd.read_csv(pg); out = {}
        for _, r in d.iterrows():
            ap = str(r.get("AllPeptides", ""))
            if ap and ap.lower() != "nan":
                out[str(r["LeadingGeneName"])] = ap.split(";")
        return out
    import pyarrow.dataset as pads
    md = P / ds.get("merged_data", "merged_data.parquet")
    if not md.exists():
        return {}
    dss = pads.dataset(str(md), format="parquet")
    r = resolve_cols(dss.schema.names, ["gene", "peptide"])
    if not r["gene"] or not r["peptide"]:
        return {}
    t = dss.to_table(columns=[r["gene"], r["peptide"]]).to_pandas().dropna()
    t.columns = ["g", "p"]
    return {g: sorted(set(s["p"])) for g, s in t.groupby("g") if len(set(s["p"])) >= MIN_PEP}


def scan_dataset(db_dir, name, cond_col, target, covariate="batch"):
    """Run the discordance scan for one dataset; write and return the results table."""
    db_dir = Path(db_dir)
    reg = load_registry(db_dir)
    ds = next((d for d in reg["datasets"] if d["name"] == name), None)
    if ds is None:
        raise SystemExit(f"dataset {name} not in registry")
    P = Path(ds["prism_dir"])
    frame = build_sample_conditions(P, ds)
    if cond_col not in frame.columns:
        print(f"[{name}] no condition '{cond_col}'; have {list(frame.columns)}"); return None
    sub = frame[frame[cond_col].notna()].copy()
    cond = (sub[cond_col].astype(str) == str(target)).astype(float)
    if cond.sum() < MIN_GRP or (len(cond) - cond.sum()) < MIN_GRP:
        print(f"[{name}] target '{target}' group too small (n={int(cond.sum())})"); return None
    batch = sub[covariate].astype(str) if covariate in sub else pd.Series("b0", index=sub.index)

    matrix = ds.get("peptide_matrix", "corrected_peptides.parquet")
    pf = pq.ParquetFile(P / matrix); key = peptide_key(pf.schema_arrow.names)
    scols = [c for c in sub.index if c in pf.schema_arrow.names]
    cond = cond.loc[scols].values
    blv = pd.get_dummies(batch.loc[scols], drop_first=True).astype(float).values
    M = pq.ParquetFile(P / matrix).read(columns=[key] + scols).to_pandas().set_index(key)
    M = np.log2(M[scols].apply(pd.to_numeric, errors="coerce").where(lambda x: x > 0)).T

    rows = []
    for g, peps in _gene_peptides(P, ds).items():
        cols = [p for p in peps if p in M.columns]
        if len(cols) < MIN_PEP:
            continue
        Y = M[cols].values
        keep = np.isfinite(Y).sum(0) >= MIN_OBS
        cols = [c for c, k in zip(cols, keep) if k]; Y = Y[:, keep]
        if len(cols) < MIN_PEP:
            continue
        yl, cl, pep_idx, samp_idx = [], [], [], []
        for j in range(len(cols)):
            m = np.isfinite(Y[:, j])
            yl.append(Y[m, j]); cl.append(cond[m])
            pep_idx.append(np.full(m.sum(), j)); samp_idx.append(np.where(m)[0])
        yl = np.concatenate(yl); cl = np.concatenate(cl)
        pep_idx = np.concatenate(pep_idx); samp_idx = np.concatenate(samp_idx)
        Bl = blv[samp_idx] if blv.shape[1] else np.empty((len(yl), 0))
        base = np.column_stack([np.ones(len(yl)), Bl, cl])
        for j in range(len(cols)):
            is_t = (pep_idx == j).astype(float)
            if is_t.sum() < MIN_OBS or (len(is_t) - is_t.sum()) < MIN_OBS:
                continue
            X = np.column_stack([base, is_t, cl * is_t])
            try:
                inv = np.linalg.pinv(X.T @ X); beta = inv @ (X.T @ yl)
                resid = yl - X @ beta; dfree = len(yl) - X.shape[1]
                if dfree < 5:
                    continue
                s2 = (resid @ resid) / dfree; se = np.sqrt(s2 * inv[-1, -1])
                tstat = beta[-1] / se; p = 2 * stats.t.sf(abs(tstat), dfree)
                rows.append((g, cols[j], _strip_mod(cols[j]), beta[-1], p))
            except Exception:
                continue
    df = pd.DataFrame(rows, columns=["gene", "peptide", "stripped", "interaction_log2", "p"]).dropna()
    if df.empty:
        print(f"[{name}] no testable peptides"); return None
    df["BH"] = (df["p"] * len(df) / df["p"].rank()).clip(upper=1.0)
    out = db_dir / "discordance"; out.mkdir(exist_ok=True)
    df.to_csv(out / f"discordance_{name}.csv", index=False)
    print(f"[{name}] {cond_col}={target} | {len(df)} peptide tests, "
          f"{(df['BH'] < BH_SIG).sum()} discordant (BH<{BH_SIG}) -> {out.name}/discordance_{name}.csv")
    return df


def aggregate(db_dir, ad_regions=None):
    """Cross-dataset reproducible discordant peptides (same peptide+direction, >=2 datasets)."""
    db_dir = Path(db_dir); out = db_dir / "discordance"
    ad_regions = set(ad_regions or [])
    recs = []
    for f in glob.glob(str(out / "discordance_*.csv")):
        nm = os.path.basename(f)[len("discordance_"):-4]
        d = pd.read_csv(f)
        s = d[d["BH"] < BH_SIG][["gene", "stripped", "interaction_log2", "BH"]].copy()
        s["dataset"] = nm; s["sign"] = np.sign(s["interaction_log2"]); recs.append(s)
    if not recs:
        print("no discordant peptides in any dataset"); return None
    allsig = pd.concat(recs)
    agg = []
    for (pep, sign), grp in allsig.groupby(["stripped", "sign"]):
        dsets = sorted(grp["dataset"].unique())
        if len(dsets) >= 2:
            agg.append({"gene": grp["gene"].iloc[0], "peptide": pep, "sign": int(sign),
                        "n_datasets": len(dsets), "datasets": ",".join(dsets),
                        "ad_regions": sum(x in ad_regions for x in dsets),
                        "min_BH": grp["BH"].min()})
    rep = pd.DataFrame(agg).sort_values(["n_datasets", "ad_regions"], ascending=False)
    rep.to_csv(out / "REPRODUCIBLE_discordant_peptides.csv", index=False)
    print(f"{len(rep)} reproducible discordant peptides -> discordance/REPRODUCIBLE_discordant_peptides.csv")
    return rep
