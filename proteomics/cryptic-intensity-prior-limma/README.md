# cryptic-intensity-prior-limma

Differential abundance of **cryptic peptides** (ALS-vs-HC style contrasts) from a PRISM
output directory, run through a **moderated linear model with an intensity-trend variance
prior** (limma-trend) — plus honest **detection (on/off) context** so a "hit" isn't
mistaken for an abundance change when it's really a presence/absence event.

Produces one tidy **CSV** you can open in Python/Excel.

## Why this exists

For near-detection-limit peptides (like cryptic peptides), the dense PRISM abundance matrix
integrates a **borrowed RT window** when the peptide isn't really detected — so an *absent*
peptide still gets a low, noisy number instead of a blank. Standard abundance tests treat
those fake numbers as real, which:

- **conflates detection loss with abundance change** (an on/off event reads as a fold-change), and
- **contaminates the variance**, so even a good variance model (the intensity-trend prior)
  can be misled.

So this tool reports the intensity-prior result **and** the plain-limma result **and** the
per-group detection rates side by side, and computes FDR both genome-wide and *within the
cryptic family*, so you can see exactly what each "hit" really is.

## What's in here

| File | What |
|------|------|
| `intensity_prior_limma.py` | the analysis script (edit the CONFIG block, run it) |
| `moderated_limma/` | **vendored fork** of the limma / moderated-linear-model code |
| `README.md` | this file |

### Vendored limma code (provenance)

`moderated_limma/` is a verbatim fork of three modules from
**proteomics-toolkit v26.4.0** (commit `25a0133`, 2026-07-02):
`statistical_analysis.py`, `normalization.py`, `preprocessing.py`. Forked so this tool is
self-contained and won't break if the toolkit changes. Exposed API:
`run_moderated_linear_model`, `StatisticalConfig`. `moderation` selects the variance prior:
`"limma"` (global, == standard limma), `"deqms"` (per-protein peptide-count prior),
`"intensity_trend"` (variance conditioned on intensity; limma-trend — used here).

## Dependencies

- `conda activate lf01` — needs numpy, scipy, statsmodels, pandas, duckdb, pyarrow.
- **PRISM loading** uses [prism-diff-explorer](https://github.com/laurenfields/prism-diff-explorer)
  for `load_prism` / `attach_clinical` / `cryptic_peptide_map`. Set `PDE_REPO` at the top of
  the script to its checkout path (default `C:\Users\field\repos\prism-diff-explorer`).
  (Not vendored — it's an actively-maintained repo, not "limma code".)

## Usage

1. `conda activate lf01`
2. Open `intensity_prior_limma.py`, edit the **CONFIG** block: `OUTPUT_DIR`, `CLINICAL_CSV`,
   and the contrast. Leave `GROUP_COLUMN = None` on the **first** run — the script prints the
   available metadata columns/values and exits so you can pick the contrast. Then set
   `GROUP_COLUMN` / `GROUP_A_VALUES` / `GROUP_B_VALUES` and re-run.
3. `python intensity_prior_limma.py`

`logFC` is group **B vs A** (positive = higher in B). Runs in a few seconds on ~40k peptides.
The LOWESS variance-trend fit uses `config.lowess_delta_frac` (default `0.01`, a fraction of
the intensity range) to interpolate between closely-spaced points; set it to `0.0` for the
exact — but ~100× slower — O(n²) fit.

## Output CSV columns

`peptide`, `is_cryptic`, `cryptic_accession`, `cryptic_protein`, `logFC`, `AveExpr`,
`n_A`, `n_B`, and:

- **intensity-trend:** `P_intensity_trend`, `FDR_genomewide_intensity_trend`, `FDR_crypticfamily_intensity_trend`
- **plain limma (comparison):** `P_limma`, `FDR_genomewide_limma`, `FDR_crypticfamily_limma`
- **detection:** `det_A`/`det_B` (# genuinely detected, DetectionQValue<0.01) and `detrate_A`/`detrate_B`

**Reading it:** if `detrate_A` and `detrate_B` differ a lot → it's a **detection (on/off)**
effect and `logFC` is *not* a true fold-change; if they're equal → a genuine abundance
difference. Sort by `P_intensity_trend`, but judge each hit with the detection columns.

## Honest caveats (learned on the Novartis/Yubin plasma cryptic data)

- At small n (e.g. 7 HC vs 18 ALS), **~20 cryptic peptides are nominally different (p<0.05)
  but ~0–2 survive family-FDR** — and the survivors are boundary-fragile (a re-processing of
  the same cohort moved them from FDR 0.037 to 0.06). Report the ~20 as **candidates**, not
  as "FDR-significant."
- The intensity-trend prior is the *right* variance model in principle, but on the dense
  cryptic matrix its apparent gains come from shrinking a **borrowed-baseline-inflated**
  variance — don't over-trust the extra hits vs plain limma here.
- The cleaner, more defensible signal is **detection-based** (a coherent "ALS loses cryptic
  peptides" pattern). Use the `detrate_*` columns / a Fisher on/off test as the headline for
  cryptic peptides; use the abundance test on genuinely-detected cells only.

## License

MIT (this repo). Vendored `moderated_limma/` retains its upstream proteomics-toolkit license.
