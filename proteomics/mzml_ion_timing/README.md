# mzml_ion_timing

Extract **actual ion injection time**, **acquisition rate**, and **cycle time**
straight from mzML scan headers — no pyteomics, no mzR, no vendor libraries.

## Why this exists

The nominal scan rate / injection time you set on the instrument is not
necessarily what you get. To report the *actual* acquisition timing you have to
read it back from the data. This is a dependency-light Python port of the UW
MacCoss lab [`manuscript-ion-counting`](https://github.com/uw-maccosslab/manuscript-ion-counting)
mzML metadata logic (`bin/HeLa_IsoWindow_*.Rmd`, originally R/`mzR`).

It streams the mzML and reads only the scan headers — peak/binary arrays are
never decoded — so it is fast (tens of thousands of scans/second) and low-memory
even on multi-hour gradients.

## Metrics (match the original R exactly)

Per MS2 scan it reads retention time, ion injection time (ms), precursor m/z and
total ion current (TIC). Per file it computes:

| column | definition |
|---|---|
| `ScanRateHz` | `n_MS2 / last_MS2_RT(s)` — actual acquisition rate |
| `MeanInjectionTimeMs` | mean ion injection time over MS2 scans |
| `MedianInjectionTimeMs` | median of the same |
| `MeanCycleTimeSec` | mean ΔRT between precursor-m/z resets (start of each DIA cycle) |
| `MedianCycleTimeSec` | median of the same |
| `NumCycles`, `ScansPerCycle` | cycle count and scans per cycle |
| `n_ms2_all`, `n_ms2_filtered`, `rt_span_sec` | scan counts and RT span |

`ScanRateHz` and the injection-time means **exclude the first MS2 after each
MS1** (the R `is_after_ms1` filter) so cycle-boundary scans don't skew the rate.

Two correctness details carried over / added:

- **RT units.** `mzR` returns retention time in seconds, but mzML stores "scan
  start time" in minutes on many instruments (Thermo). The unit is read from the
  cvParam and converted, so `ScanRateHz` is not silently 60× off.
- **Precursor m/z.** Prefers "selected ion m/z" (`MS:1000744`) and falls back to
  "isolation window target m/z" (`MS:1000827`) when the former is absent, for
  broader converter/vendor compatibility.

## Dependencies

- `numpy`, `pandas`
- `pyarrow` — only if you write per-scan output to `.parquet`

```
pip install numpy pandas pyarrow
```

## Usage

```bash
# one summary row per file (directory scanned recursively)
python mzml_ion_timing.py DATA_DIR --out summary.csv

# also dump every MS2 scan (RT, injection time, TIC, precursor m/z, Ions)
python mzml_ion_timing.py DATA_DIR --out summary.csv --per-scan-out scans.parquet

# explicit files instead of a directory
python mzml_ion_timing.py a.mzML b.mzML --out summary.csv

# skip a subfolder; single-threaded; shallow scan
python mzml_ion_timing.py DATA_DIR --out summary.csv --exclude raw_files --workers 1 --no-recursive
```

Per-scan output also includes `Ions = injectionTime * TIC / 1000`, the
ion-count proxy from the source manuscript.

### Pulling experimental factors out of file names

The tool is instrument- and naming-scheme-agnostic: it keys results by file
name and assumes nothing about them. To add metadata columns, pass a regex with
**named groups** via `--regex` — each named group becomes a summary column.

```bash
# e.g. Stellar files like 600to700_24m_run94_200SR_TR1.mzML
python mzml_ion_timing.py DATA_DIR --out summary.csv \
  --regex '(?P<window>\d+to\d+)_(?P<gradient>\d+)m_run\d+_(?P<rate>\d+)SR(?:_TR(?P<replicate>\d+))?'
```

adds `window`, `gradient`, `rate`, `replicate` columns. Groups that don't match
are left blank.

## Use as a library

```python
import mzml_ion_timing as mit

summary, per_scan = mit.collect("DATA_DIR", per_scan=True,
                                regex=r"(?P<rate>\d+)SR")
# summary: one row per file; per_scan: one row per MS2 scan (or None)
```

`parse_ms2_headers(path)`, `compute_metrics(ms2_df)`, and `process_file(path)`
are also exposed for finer-grained use.

## Notes / limitations

- Reads the standard PSI-MS accessions listed above; unusual converters that
  omit "ion injection time" (`MS:1000927`) will yield `NaN` injection times.
- Cycle detection assumes precursor m/z increases monotonically within a DIA
  cycle and resets at the next cycle (a negative m/z step marks a new cycle).
  This is the manuscript's heuristic; it fits GPF / stepped-window DIA. Randomized
  or non-monotonic isolation schedules would need a different cycle definition.
- MS1-only files or files with no MS2 scans produce `NaN` metrics.
