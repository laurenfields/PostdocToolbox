# dark-tic-screen

Screen a DIA run for the ion current that the search did **not** assign to identified
peptides — the "dark TIC" — and decide whether any of it is worth chasing as a
proteoform, PTM, or novel species.

A library/targeted search only ever assigns signal to peptides already in its space, so
it never tells you what the *rest* of the ion current is. This tool measures that
remainder and, crucially, distinguishes **high-intensity genuinely-unexplained** signal
(a real lead) from **spillover** (an identified peptide's own fragment showing up at
another retention time) and from a **diffuse low-abundance haze** (nothing to chase).

## The three tiers

1. **Coarse** — fraction of MS2 ion current in scans with no eluting identification in
   their isolation window. Needs only the ID list + raw scan headers, so it is fast and
   runs across every file. A cohort-wide "how much did we assign" flag.
2. **Fine** — after subtracting each eluting identified peptide's own fragment ladder
   (monoisotopic b/y at charge 1–2, ¹³C isotopes, H₂O/NH₃ neutral losses, and the
   surviving precursor), how concentrated is the residual dark? Reported as the share of
   the dark carried by the top-1 / top-500 features and the number of features needed to
   reach 50 %. A few big features vs a haze of tiny ones.
3. **Novel + delta-mass** — the highest-intensity dark features that are *not* an
   identified fragment at another RT, each tested for whether it equals an identified
   fragment plus a common modification mass (deamidation, methyl, oxidation, water loss,
   acetyl, phospho, …). Because a narrow window's fragment set is dense, the same test is
   run with matched **decoy deltas**; a real proteoform signal requires the real-mod match
   rate to beat the decoy rate.

## Usage

```bash
python dark_tic_screen.py --report ids.parquet --raws /path/to/raws --out screen/
```

`--report` is a Skyline or DIA-NN export (`.parquet` or `.csv`) of **confident**
identifications, one row per (replicate, precursor, charge) with its elution boundaries.
Column names are auto-detected from common Skyline/DIA-NN headers; override any that do
not match with `--col-seq / --col-pmz / --col-pch / --col-rep / --col-start / --col-end`.
By default a row counts as confident only if a q-value column is present and populated
(Skyline emits q ≤ 0.01); pass `--all-ids` if the report is already filtered.

Modified sequences are parsed in both `C(UniMod:4)` / `(unimod:4)` and bracket-mass
`M[+15.9949]` notations, with fixed and variable mods handled per residue and at the
termini.

Key options:

| flag | meaning |
|---|---|
| `--detail STEM` | run the fine + novel screen on this raw stem (default: the run with the most MS2 TIC) |
| `--all-fine` | run the fine screen on every file (slow: full peak decode per file) |
| `--rt-range lo,hi` | restrict to a retention-time window in minutes (e.g. `5,20`) |
| `--ppm` | match tolerance (default 10) |
| `--top` | number of top features inspected for novelty (default 200) |
| `--no-rt-gate` | report has no per-run RT peak boundaries: subtract each ID's fragments across all RT in its window, giving a **conservative lower-bound** dark (the coarse tier is not meaningful in this mode). Use when only a precursor list is available (e.g. a Panorama abundance table). |
| `--out DIR` | write `coarse_summary.parquet` (per-run assigned/dark fractions) |

Modified-sequence masses shown rounded by Skyline (`C[+57]`, `M[+16]`) are snapped to the
exact modification mass, so a report exported with integer mass display still fragments
correctly.

### Example (Skyline document)

```bash
python dark_tic_screen.py \
  --report Latimer_PlasmaEVs_Cryptic.parquet \
  --raws  /data/astral_raws \
  --rt-range 5,20 --out darkscreen/
```

## Reading the verdict

For the detailed run the tool prints a concentration summary and a one-line verdict:

- **NEGATIVE** — the dark is diffuse (no single unexplained feature reaches ~0.5 % of the
  dark) and/or the biggest features are spillover, and the real-mod delta rate does not
  beat the decoy rate. Nothing high-abundance or proteoform-like is hiding; do not spend
  effort on an open search here.
- **FLAG** — a high-intensity unexplained feature is present *and* real modification
  offsets match more than decoy deltas. Worth a targeted open / delta-mass search
  (e.g. MSFragger open search) to name the proteoform.

The coarse per-run table (`--out`) is the cohort screen: rank runs by `dark_tic_frac` and
drop to the fine screen only on the ones that flag high *and* have raw files available.

## What it needs, and does not

- **Needs the raw files** for the fine and novel tiers. Skyline's `.skyd` caches only the
  extracted chromatograms of targeted transitions, not the full unassigned MS2 peaks, so
  the dark peak content can only come from the raw data. The coarse tier needs only scan
  headers (TIC per scan), so it is cheap.
- **Accepts `.raw` (Thermo) or `.mzML`.** `.raw` is read via the `thermo_raw` reader
  (ProteoWizard DLLs); `.mzML` via a built-in stdlib reader (no extra dependency).
  **mzML must be centroided** (peak-picked on conversion) — the dark accounting assumes
  centroids, and the tool warns once if a spectrum looks profile. For **overlapping /
  staggered DIA** (e.g. `ovlp` acquisitions), prefer a **demultiplexed mzML**: demux
  yields clean windows that every MS2 scan maps to, whereas the raw's first-cycle window
  inference cannot resolve staggered windows (they land in the "unassigned" slice).
- Built for **narrow-window DIA** (isolation windows inferred per file from the first
  cycle via the `thermo_raw` reader). Wide-window DIA and DDA are out of scope.
- The delta-mass tier is a **screen, not an identification**: a fragment-level delta match
  is suggestive, and the decoy calibration is what keeps it honest. Confirm any flag with
  a proper open/targeted search.

## Running from PRISM / Skyline transition reports (e.g. `merged_data.parquet`)

Point `--report` straight at a **transition-level** Skyline/PRISM parquet (one row per
fragment, with a `FragmentIon`/`ProductMz` column, such as PRISM's `merged_data.parquet`).
Both tools auto-detect it and roll it up per-precursor with duckdb (streaming, so the
88M-row file never loads into memory): `min(StartTime)`/`max(EndTime)`, `PrecursorMz`, and
a confident-ID cutoff of `DetectionQValue <= --qvalue-max` (default 0.01) for the requested
run. This reproduces a Skyline peak-boundary report exactly, and it is the source your PRISM
notes point to for real detection.

```bash
python dark_tic_dashboard.py --report prism_output_dir/merged_data.parquet --raw run.raw --out dash.html
```

Do **not** use PRISM's `corrected_*.parquet` (the rolled-up abundance matrix): it has no
precursor m/z and no per-run RT boundaries, and its cells are imputed - it cannot drive the
assigned-vs-dark accounting or the RT map.

The run is matched to a `ReplicateName` by the raw filename, in either direction (the raw
filename often carries an extra date/instrument prefix that the Skyline replicate name
omits, or vice versa). If the filename does not disambiguate, pass the exact name with
`--replicate "<ReplicateName>"`; on a no-match the error lists the available names.

## dark_tic_dashboard.py - one-command visual dashboard

`dark_tic_dashboard.py --report <report> --raw <one .raw> --out dash.html` writes a single
self-contained HTML page (charts embedded, no server) for one Skyline+raw pair: where the
ion current goes (MS1/MS2, assigned vs dark), a retention-time x window darkness map, the
dark concentration curve, the top unexplained features with delta-mass flags, and a
FLAG/NEGATIVE verdict. It auto-opens the page unless `--no-open`. Same `--col-*`,
`--no-rt-gate`, `--qvalue-max`, and PRISM auto-detection as the screen.

## Dependencies

- Python ≥ 3.9, `numpy`, `pyarrow`, `duckdb` (transition-report rollup), `matplotlib` (dashboard).
- The sibling PostdocToolbox reader `proteomics/thermo_raw/thermo_raw_reader.py`
  (pythonnet + ProteoWizard's ThermoFisher CommonCore DLLs). Point it at your install
  with `THERMO_RAW_DLL_DIR` if needed. Thermo `.raw` only.

## Provenance

Generalized from the Latimer plasma-EV dark-TIC investigation (stages 01/02 coarse,
03 subtraction, 27 concentration, 28 delta-mass). Validated on external Panorama Public
data (MacCoss mouse plasma-EV Astral DIA, PXD080860): a clean NEGATIVE verdict (diffuse
dark, no high-abundance novel feature, no proteoform signal). MIT licensed.
