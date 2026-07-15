# PR notes — LOWESS `delta` for the moderated-LM variance priors

Draft description for an upstream PR to `uw-maccosslab/proteomics-toolkit`
(target base: the tip `25a0133` "Vectorise dense per-feature OLS"). Not yet opened.

---

## Title
`statistical_analysis: ~100x faster intensity_trend/deqms variance priors via LOWESS delta`

## Summary
The `intensity_trend` (and `deqms`) variance priors spend nearly all of their
runtime inside a single `statsmodels` `lowess` call. With the default `delta=0`,
LOWESS runs a full local regression at **every** input point, which is O(n²) and
dominates on peptide-level data — a ~40k-peptide unpaired fit produces ~80k
(feature, group) anchor points and the trend fit alone takes **~4–5 minutes**.

This adds `StatisticalConfig.lowess_delta_frac` (default `0.01`) and threads it
into both variance-prior LOWESS fits as `delta = frac * range(x)`. With `delta>0`,
LOWESS fits only at anchor points spaced ≥ `delta` apart and linearly
interpolates between them.

## Why `0.01` is the right default
`delta = 0.01 * range(x)` is **statsmodels' own documented recommendation** for
large n (`smoothers_lowess.py`) and is **R `lowess()`'s built-in default**
(`delta = 0.01 * diff(range(x))`). So this moves the code *toward* the reference
implementations' default behavior, not away from it.

The approximation is near-lossless here because the variance-vs-intensity trend
is smooth and monotone and the span (`frac=0.5`) is far wider than a 1%-of-range
window — no fine structure is skipped.

## Evidence (real 41,588-peptide unpaired ALS-vs-HC matrix)
| metric | `delta=0.0` (exact) | `delta=0.01` (default) |
|---|---|---|
| intensity_trend fit | ~250 s | ~2.6 s (**~95× faster**) |
| logFC | — | **bit-identical** (contrast never uses the prior) |
| max \|ΔP\| over all peptides | — | **1.3e-4** |
| significance-call flips (p<0.05) | — | **0** |
| prior hyperparameters | `d0=5.43`, `s0²=0.903` | unchanged |

## Behavior-change note
The default `0.01` slightly changes numeric output vs the old exact fit on real
data (max \|ΔP\| ~1e-4; no significance-count changes observed). On the existing
synthetic test matrices the effect is nil, so all prior tests pass unchanged.
Set `lowess_delta_frac = 0.0` to force the exact O(n²) fit and reproduce the old
numbers to the last digit — the escape hatch is wired and tested.

## Changes
- `StatisticalConfig.lowess_delta_frac = 0.01` (new, documented).
- `_lowess_delta(x, config)` helper: `frac * (max(x) - min(x))`, returning `0.0`
  (exact) for non-finite/≤0 `frac` or empty/degenerate `x`.
- Threaded `delta=` into the two `lowess()` calls (`_fit_intensity_trend_prior`,
  `_fit_count_dependent_prior`); `_fit_count_dependent_prior` gains a
  `config=None` parameter (default reads as `0.01`).
- Tests: `tests/test_lowess_delta.py` — helper edge cases + an equivalence
  assertion that `delta=0.0` vs `0.01` agree to a tight tolerance with logFC
  identical and no significance flips.

## Backward compatibility
Purely additive API. No signature changes except the new optional `config`
parameter on the private `_fit_count_dependent_prior`. Existing callers that
never set `lowess_delta_frac` transparently get the new default.
