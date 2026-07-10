"""
Cryptic-peptide differential analysis with an intensity-trend variance prior (limma-trend).

Runs a moderated linear model TWO ways for comparison and writes one CSV:
  - moderation="limma"            (global variance prior; == standard limma / inmoose)
  - moderation="intensity_trend"  (variance prior conditioned on intensity; limma-trend)

Also flags cryptic peptides and adds honest per-group DETECTION rates (from the
transition-level DetectionQValue), because for near-detection-limit cryptic peptides the
dense abundance matrix conflates "absent" (borrowed baseline) with "low abundance".
See README.md for the methodology story (why detection rate matters).

The limma code lives in the vendored `moderated_limma/` package (forked from
proteomics-toolkit v26.4.0). PRISM loading uses prism-diff-explorer (set PDE_REPO below).

USAGE
-----
1. conda activate lf01   (needs numpy/scipy/statsmodels/pandas/duckdb/pyarrow + inmoose)
2. Set the CONFIG block. On the FIRST run leave GROUP_COLUMN = None: the script prints the
   available metadata columns + values and exits, so you can choose the contrast. Then set
   GROUP_COLUMN / GROUP_A_VALUES / GROUP_B_VALUES and re-run.
3. python intensity_prior_limma.py
   (the intensity_trend fit takes a few minutes on ~40k peptides — that's expected.)

logFC is group B vs group A (positive = higher in B).
"""
import os
import sys

import numpy as np
import pandas as pd

# vendored limma / moderated-linear-model code (self-contained fork)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from moderated_limma import StatisticalConfig, run_moderated_linear_model

# PRISM loading / cryptic mapping comes from prism-diff-explorer (adjust if it moves)
PDE_REPO = r"C:\Users\field\repos\prism-diff-explorer"
sys.path.insert(0, PDE_REPO)
import prism_diff_explorer as pde

# ============================== CONFIG ================================================
OUTPUT_DIR = r"G:\Manuscripts\CrypticPeptide\YubinPlasmaAnalysis\output_dir_2026June18"
CLINICAL_CSV = r"G:\Manuscripts\CrypticPeptide\YubinPlasmaAnalysis\Replicates_formatted.csv"
CLINICAL_ID_COLUMN = None          # None = auto-detect the sample-id column by value;
                                   # set to the exact column name if auto-detect picks wrong.
LEVEL = "peptide"                  # cryptic peptides live at peptide level

# --- the contrast: fill these in after the first run prints the options ---------------
GROUP_COLUMN = None                # e.g. "Sample_Key"
GROUP_A_VALUES = []                # reference group, e.g. ["Healthy control"]  (logFC baseline)
GROUP_B_VALUES = []                # case group,      e.g. ["ALS"]              (logFC = B vs A)

CRYPTIC_ONLY = False               # True = write only cryptic peptides; False = all
OUT_CSV = os.path.join(os.path.dirname(OUTPUT_DIR) or ".", "cryptic_intensity_prior_results.csv")
# ======================================================================================


def bh(pvals):
    """Benjamini-Hochberg FDR (NaN-safe)."""
    return pde._bh(np.asarray(pvals, dtype=float))


def merged_columns(merged_parquet):
    """Auto-detect the peptide / detection-q / sample-id column names (spaced vs not)."""
    import duckdb
    cols = duckdb.connect().execute(
        f"DESCRIBE SELECT * FROM read_parquet('{merged_parquet}')").df()["column_name"].tolist()

    def pick(cands):
        for c in cands:
            if c in cols:
                return c
        return None
    return {
        "pep": pick(["PeptideModifiedSequenceUnimodIds", "Peptide Modified Sequence Unimod Ids"]),
        "prot": pick(["Protein"]),
        "detq": pick(["DetectionQValue", "Detection Q Value"]),
        "samp": pick(["Sample ID", "Sample Id", "Replicate Name"]),
    }


def main():
    print(f"Loading PRISM ({LEVEL}) from:\n  {OUTPUT_DIR}")
    d = pde.load_prism(OUTPUT_DIR, LEVEL)
    expr = d["expr_log2"]                       # peptides x samples, log2 (dense, linear->log2)
    meta, cinfo = pde.attach_clinical(d["meta"], CLINICAL_CSV, CLINICAL_ID_COLUMN)
    print(f"  matrix: {expr.shape[0]:,} peptides x {expr.shape[1]} samples")
    print(f"  clinical join: {cinfo['match_rate']:.0%} matched via "
          f"'{cinfo['key_column']}' ({'manual' if cinfo.get('manual') else 'auto'})")

    # --- if the contrast isn't set, show the options and stop -------------------------
    meta_cols = [c for c in meta.columns if c not in pde._NON_META_COLS and meta[c].notna().any()]
    if not GROUP_COLUMN or not GROUP_A_VALUES or not GROUP_B_VALUES:
        print("\n>>> Set GROUP_COLUMN / GROUP_A_VALUES / GROUP_B_VALUES, then re-run.")
        print(">>> Available metadata columns and their values (experimental samples):\n")
        exp = meta[meta["sample_type"] == "experimental"] if "sample_type" in meta.columns else meta
        for c in meta_cols:
            vc = exp[c].astype(str).value_counts()
            if 2 <= len(vc) <= 25:              # skip free-text / continuous columns
                print(f"  {c!r}:")
                for v, n in vc.items():
                    print(f"       {v!r}: {n}")
        sys.exit(0)

    # --- define the two groups --------------------------------------------------------
    col = meta[GROUP_COLUMN].astype(str)
    group_a = [s for s in meta.index[col.isin([str(v) for v in GROUP_A_VALUES])] if s in expr.columns]
    group_b = [s for s in meta.index[col.isin([str(v) for v in GROUP_B_VALUES])] if s in expr.columns]
    samples = group_a + group_b
    if len(group_a) < 2 or len(group_b) < 2:
        sys.exit(f"Need >=2 samples per group; got A={len(group_a)}, B={len(group_b)}.")
    print(f"\nContrast on {GROUP_COLUMN}:  A(ref)={GROUP_A_VALUES} n={len(group_a)}  |  "
          f"B={GROUP_B_VALUES} n={len(group_b)}   (logFC = B vs A)")

    # --- cryptic peptide map (accession/protein) --------------------------------------
    merged = os.path.join(OUTPUT_DIR, "merged_data.parquet")
    mc = merged_columns(merged)
    cmap = {}
    try:
        cmap = pde.cryptic_peptide_map(merged, "CRYPTIC", peptide_col=mc["pep"], protein_col=mc["prot"])
    except Exception as e:  # noqa: BLE001
        print(f"  (cryptic map skipped: {e})")
    crypt_ids = set(cmap)
    print(f"  cryptic peptides identified: {len(crypt_ids & set(expr.index))}")

    def acc(pep):
        s = str(cmap.get(pep, ""))
        return s.split("|")[1] if "|" in s else ""

    # --- per-group detection rates (genuine detection, not borrowed baseline) ---------
    det_a = det_b = None
    try:
        import duckdb
        det = duckdb.connect().execute(f"""
            SELECT "{mc['pep']}" AS pep, "{mc['samp']}" AS samp,
                   MAX(CASE WHEN "{mc['detq']}" IS NOT NULL AND "{mc['detq']}" < 0.01
                            THEN 1 ELSE 0 END) AS det
            FROM read_parquet('{merged}') GROUP BY pep, samp
        """).df().pivot_table(index="pep", columns="samp", values="det", fill_value=0)
        dcols = set(det.columns)
        keycol = max(["sample_id", "sample"], key=lambda c: len(set(meta[c].astype(str)) & dcols))
        id2val = {}
        for sid in samples:
            k = str(meta.loc[sid, keycol])
            if k in dcols:
                id2val[sid] = k
        a_cols = [id2val[s] for s in group_a if s in id2val]
        b_cols = [id2val[s] for s in group_b if s in id2val]
        det_a = det[a_cols].sum(axis=1) if a_cols else None
        det_b = det[b_cols].sum(axis=1) if b_cols else None
        print(f"  detection matrix matched via '{keycol}': A cols {len(a_cols)}, B cols {len(b_cols)}")
    except Exception as e:  # noqa: BLE001
        print(f"  (detection rates skipped: {e})")

    # --- run the moderated model both ways --------------------------------------------
    feat_log2 = expr[samples].copy()
    feat_log2.index.name = "Protein"
    feat_raw = 2 ** feat_log2                                   # raw linear, for the intensity trend
    metadata_df = pd.DataFrame({
        "Sample": samples,
        "GROUP": ["A"] * len(group_a) + ["B"] * len(group_b),
    })

    def run(moderation):
        cfg = StatisticalConfig()
        cfg.analysis_type = "unpaired"
        cfg.statistical_test_method = "moderated_linear_model"
        cfg.group_column = "GROUP"
        cfg.group_labels = ["A", "B"]                          # logFC = B vs A
        cfg.moderation = moderation
        cfg.robust = False
        if moderation == "intensity_trend":
            cfg._raw_feature_data = feat_raw                   # trend prior sees raw intensity
        import time
        t0 = time.perf_counter()
        res = run_moderated_linear_model(feat_log2, metadata_df, cfg).set_index("Protein")
        print(f"    {moderation}: {time.perf_counter() - t0:.0f}s, {len(res):,} features")
        return res

    print("\nFitting moderated linear models:")
    res_l = run("limma")
    res_it = run("intensity_trend")

    # --- assemble the output table ----------------------------------------------------
    feats = res_it.index
    out = pd.DataFrame(index=feats)
    out["peptide"] = feats
    out["is_cryptic"] = [f in crypt_ids for f in feats]
    out["cryptic_accession"] = [acc(f) for f in feats]
    out["cryptic_protein"] = [str(cmap.get(f, "")) for f in feats]
    out["logFC"] = res_it["logFC"]
    out["AveExpr"] = res_it["AveExpr"]
    out["n_A"] = len(group_a)
    out["n_B"] = len(group_b)

    # intensity-trend (limma-trend) columns
    out["P_intensity_trend"] = res_it["P.Value"]
    out["FDR_genomewide_intensity_trend"] = bh(res_it["P.Value"].to_numpy())
    # plain-limma columns for comparison
    out["P_limma"] = res_l["P.Value"].reindex(feats)
    out["FDR_genomewide_limma"] = bh(out["P_limma"].to_numpy())

    # cryptic-family FDR (BH within the cryptic set only) for each engine
    for src, dst in [("P_intensity_trend", "FDR_crypticfamily_intensity_trend"),
                     ("P_limma", "FDR_crypticfamily_limma")]:
        out[dst] = np.nan
        cmask = out["is_cryptic"].to_numpy()
        cp = out.loc[cmask, src]
        out.loc[cmask, dst] = bh(cp.to_numpy())

    # detection rates
    if det_a is not None and det_b is not None:
        out["det_A"] = det_a.reindex(feats).fillna(0).astype(int)
        out["det_B"] = det_b.reindex(feats).fillna(0).astype(int)
        out["detrate_A"] = (out["det_A"] / max(len(group_a), 1)).round(3)
        out["detrate_B"] = (out["det_B"] / max(len(group_b), 1)).round(3)

    if CRYPTIC_ONLY:
        out = out[out["is_cryptic"]]
    out = out.sort_values("P_intensity_trend")

    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(out):,} rows -> {OUT_CSV}")

    # --- quick summary ---------------------------------------------------------------
    cr = out[out["is_cryptic"]]
    print(f"\nCryptic peptides: {len(cr)}")
    print(f"  intensity_trend: raw p<0.05 = {(cr['P_intensity_trend']<0.05).sum()}, "
          f"family-FDR<0.05 = {(cr['FDR_crypticfamily_intensity_trend']<0.05).sum()}")
    print(f"  limma          : raw p<0.05 = {(cr['P_limma']<0.05).sum()}, "
          f"family-FDR<0.05 = {(cr['FDR_crypticfamily_limma']<0.05).sum()}")
    show = ["cryptic_accession", "logFC", "P_intensity_trend",
            "FDR_crypticfamily_intensity_trend"]
    if "detrate_A" in out.columns:
        show += ["detrate_A", "detrate_B"]
    print("\nTop 15 cryptic by intensity-trend p:")
    print(cr[show].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
