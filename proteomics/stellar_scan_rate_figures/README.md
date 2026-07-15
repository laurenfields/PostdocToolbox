# stellar_scan_rate_figures

Notebook that builds eight publication figures comparing Thermo **Stellar** DIA
scan rates (67 / 125 / 200 kDa/s) under GPF and Wide-window acquisition on a
24-min gradient.

> **Project-specific example, not a general utility.** Paths and run layout are
> hardcoded for the StellarDIA ScanRateTest dataset (mzML files + Skyline report
> parquets are *not* included in this repo). It lives here as a worked reference
> for two reusable patterns: (1) violin/box distributions from Skyline `ions` /
> *Points across the peak* reports, and (2) per-file acquisition-timing metrics
> pulled straight from mzML scan headers. For the standalone, generalized timing
> extractor see [`../mzml_ion_timing`](../mzml_ion_timing).

## Figures
| # | Figure | Source |
|---|--------|--------|
| 1 | Apex total ion count (n>0 annotated) | Skyline `ions` report |
| 2 | LC peak analyte ion count | Skyline `ions` report |
| 3 | Points across the peak | Skyline *Points across the peak* report |
| 4 | Actual injection time (axis zoomed 0.2–60 ms) | mzML headers |
| 5 | Actual acquisition rate (GPF colored by m/z window) | mzML headers |
| 6 | Acquisition rate vs. injection time | mzML headers |
| 7 | Overhead vs. injection time | mzML headers |
| 8 | Cycle time vs. injection time | mzML headers |

Each figure cell is followed by a markdown cell documenting the exact
calculation and source, and a **Method parameters** cell at the top lists the
instrument settings read from the mzML.

## Key calculations
- **Acquisition rate (Hz)** = total spectra in the mzML ÷ total acquisition time.
- **Overhead (ms)** = (1000 × total time ÷ n MS2) − actual mean injection time.
- **Cycle time (s)** = mean ΔRT between precursor-m/z resets.
- Injection time, precursor m/z, TIC and RT are read header-only (no peak
  decoding); RT is converted minute→second.

## Running
Point `STELLAR` (in the setup cell) at the mzML root and the `RUNS` / `POINTS`
tables at the Skyline report parquets, then run top to bottom. mzML metrics are
cached to `output/`; figures are written to `Figures/` as PNG (300 dpi) + SVG.
Colour = scan rate (Okabe–Ito); scheme = facet or marker shape.

Requires: numpy, pandas, pyarrow, seaborn, matplotlib.
