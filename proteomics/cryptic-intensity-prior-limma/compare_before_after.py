"""
Before/after check for the intensity_trend LOWESS speedup.

Runs the intensity_trend moderated fit TWO ways on the exact same data/contrast:
  - "before": config.lowess_delta_frac = 0.0   -> exact O(n^2) LOWESS (the original, ~4-5 min)
  - "after" : config.lowess_delta_frac = 0.01  -> delta-interpolated LOWESS (the fast one, ~seconds)

It reuses the CONFIG you already set in intensity_prior_limma.py (OUTPUT_DIR, CLINICAL_CSV,
GROUP_COLUMN, GROUP_A/B_VALUES), writes a side-by-side comparison CSV, and prints whether any
result actually moved (runtimes, cryptic significance counts, max |delta P|, logFC agreement,
and the number of p<0.05 / FDR<0.05 calls that flip).

Run:  conda activate lf01  &&  python compare_before_after.py
"""
import os
import sys
import time

import numpy as np
import pandas as pd

# Reuse the tool module: its top-level sets sys.path for moderated_limma + prism-diff-explorer,
# defines merged_columns()/bh(), and holds the CONFIG constants you edited.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import intensity_prior_limma as ipl  # noqa: E402
from moderated_limma import StatisticalConfig, run_moderated_linear_model  # noqa: E402

pde = ipl.pde
OUT_CSV = os.path.join(os.path.dirname(ipl.OUT_CSV) or ".", "compare_before_after.csv")

# --- load exactly what the tool would load -------------------------------------------------
d = pde.load_prism(ipl.OUTPUT_DIR, ipl.LEVEL)
expr = d["expr_log2"]
meta, cinfo = pde.attach_clinical(d["meta"], ipl.CLINICAL_CSV, ipl.CLINICAL_ID_COLUMN)
print(f"matrix: {expr.shape[0]:,} peptides x {expr.shape[1]} samples "
      f"(clinical join {cinfo['match_rate']:.0%} via '{cinfo['key_column']}')")

col = meta[ipl.GROUP_COLUMN].astype(str)
group_a = [s for s in meta.index[col.isin([str(v) for v in ipl.GROUP_A_VALUES])] if s in expr.columns]
group_b = [s for s in meta.index[col.isin([str(v) for v in ipl.GROUP_B_VALUES])] if s in expr.columns]
samples = group_a + group_b
print(f"contrast {ipl.GROUP_COLUMN}: A={ipl.GROUP_A_VALUES} n={len(group_a)} | "
      f"B={ipl.GROUP_B_VALUES} n={len(group_b)}")

feat_log2 = expr[samples].copy()
feat_log2.index.name = "Protein"
feat_raw = 2 ** feat_log2
metadata_df = pd.DataFrame({"Sample": samples, "GROUP": ["A"] * len(group_a) + ["B"] * len(group_b)})

# cryptic set (for the significance counts / family FDR)
merged = os.path.join(ipl.OUTPUT_DIR, "merged_data.parquet")
mc = ipl.merged_columns(merged)
cmap = pde.cryptic_peptide_map(merged, "CRYPTIC", peptide_col=mc["pep"], protein_col=mc["prot"])
crypt_ids = set(cmap)


def run_intensity_trend(delta_frac):
    cfg = StatisticalConfig()
    cfg.analysis_type = "unpaired"
    cfg.statistical_test_method = "moderated_linear_model"
    cfg.group_column = "GROUP"
    cfg.group_labels = ["A", "B"]
    cfg.moderation = "intensity_trend"
    cfg.robust = False
    cfg.lowess_delta_frac = delta_frac
    cfg._raw_feature_data = feat_raw
    t0 = time.perf_counter()
    res = run_moderated_linear_model(feat_log2, metadata_df, cfg).set_index("Protein")
    return res, time.perf_counter() - t0


print("\nRunning intensity_trend BEFORE (delta=0.0, exact)  -- this is the slow one ...")
before, t_before = run_intensity_trend(0.0)
print(f"  before: {t_before:.1f}s")
print("Running intensity_trend AFTER  (delta=0.01, fast) ...")
after, t_after = run_intensity_trend(0.01)
print(f"  after : {t_after:.1f}s   (speedup {t_before / max(t_after, 1e-9):.0f}x)")

# --- side-by-side table --------------------------------------------------------------------
feats = before.index
cmp = pd.DataFrame(index=feats)
cmp["peptide"] = feats
cmp["is_cryptic"] = [f in crypt_ids for f in feats]
cmp["logFC_before"] = before["logFC"]
cmp["logFC_after"] = after["logFC"].reindex(feats)
cmp["P_before"] = before["P.Value"]
cmp["P_after"] = after["P.Value"].reindex(feats)
cmp["dP"] = (cmp["P_after"] - cmp["P_before"]).abs()

# cryptic-family FDR each way
for tag, src in [("before", before), ("after", after)]:
    fam = pd.Series(np.nan, index=feats)
    cm = cmp["is_cryptic"].to_numpy()
    fam.loc[cm] = ipl.bh(src["P.Value"].reindex(feats)[cm].to_numpy())
    cmp[f"FDRfam_{tag}"] = fam

cmp = cmp.sort_values("P_after")
cmp.to_csv(OUT_CSV, index=False)

# --- equivalence report --------------------------------------------------------------------
cr = cmp[cmp["is_cryptic"]]


def counts(pcol, fcol):
    return int((cr[pcol] < 0.05).sum()), int((cr[fcol] < 0.05).sum())


rb, fb = counts("P_before", "FDRfam_before")
ra, fa = counts("P_after", "FDRfam_after")
flip_p = int(((cmp["P_before"] < 0.05) != (cmp["P_after"] < 0.05)).sum())
flip_famcr = int(((cr["FDRfam_before"] < 0.05) != (cr["FDRfam_after"] < 0.05)).sum())
max_dp = float(np.nanmax(cmp["dP"].to_numpy()))
max_dfc = float(np.nanmax((cmp["logFC_after"] - cmp["logFC_before"]).abs().to_numpy()))

print("\n================ BEFORE vs AFTER ================")
print(f"runtime            : {t_before:6.1f}s  ->  {t_after:6.1f}s   ({t_before / max(t_after, 1e-9):.0f}x faster)")
print(f"cryptic raw p<0.05 : {rb:6d}    ->  {ra:6d}")
print(f"cryptic famFDR<0.05: {fb:6d}    ->  {fa:6d}")
print(f"max |dP| (all pep) : {max_dp:.2e}")
print(f"max |dlogFC|       : {max_dfc:.2e}   (0 = identical)")
print(f"p<0.05 calls flipped (all peptides) : {flip_p}")
print(f"famFDR<0.05 calls flipped (cryptic) : {flip_famcr}")
verdict = "IDENTICAL results" if (rb == ra and fb == fa and flip_famcr == 0) else "RESULTS DIFFER -- inspect"
print(f"VERDICT            : {verdict}")
print(f"\nwrote side-by-side table -> {OUT_CSV}")
