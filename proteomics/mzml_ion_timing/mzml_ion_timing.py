"""Extract actual injection time and acquisition rate from mzML scan headers.

A dependency-light port of the UW MacCoss lab ``manuscript-ion-counting`` mzML
metadata logic (bin/HeLa_IsoWindow_*.Rmd, originally R/mzR), for reading real
per-scan acquisition timing straight out of mzML without pyteomics or mzR.

For every MS2 spectrum it reads, from the scan header only (peak/binary arrays
are never decoded, so it is fast and low-memory): retention time, ion injection
time (ms), precursor m/z and total ion current (TIC). Per file it computes,
matching the original R:

  ScanRateHz          = n_MS2 / last_MS2_RT(s)           (actual acquisition rate)
  MeanInjectionTimeMs = mean ion injection time over MS2 scans
      Both exclude the first MS2 after each MS1 (the R ``is_after_ms1`` filter),
      so cycle-boundary scans don't skew the rate.
  MeanCycleTimeSec    = mean ΔRT between precursor-m/z resets (start of each cycle)
  Ions                = injectionTime * TIC / 1000        (per scan)

mzR reports retention time in seconds; mzML stores "scan start time" in minutes
on many instruments, so the unit is read and converted here.

This tool is instrument/vendor-agnostic: it keys results by file name and does
NOT assume any particular naming scheme. To pull experimental factors (scan
rate, gradient, window, replicate, ...) out of your file names, pass a regex
with named groups via ``--regex`` — each named group becomes a summary column.

CLI
---
    # one summary row per file
    python mzml_ion_timing.py DATA_DIR --out summary.csv

    # also dump every MS2 scan to parquet
    python mzml_ion_timing.py DATA_DIR --out summary.csv --per-scan-out scans.parquet

    # pull metadata columns out of the file names
    python mzml_ion_timing.py DATA_DIR --out summary.csv \
        --regex '(?P<window>\\d+to\\d+)_(?P<grad>\\d+)m_run\\d+_(?P<rate>\\d+)SR'

    # explicit files instead of a directory
    python mzml_ion_timing.py a.mzML b.mzML --out summary.csv

Dependencies: numpy, pandas (+ pyarrow only if writing parquet per-scan output).
"""

from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

# mzML controlled-vocabulary accessions read from each spectrum header.
ACC_MS_LEVEL = "MS:1000511"
ACC_SCAN_START = "MS:1000016"    # retention time
ACC_INJ_TIME = "MS:1000927"      # ion injection time (ms)
ACC_TIC = "MS:1000285"           # total ion current
ACC_SELECTED_MZ = "MS:1000744"   # selected ion m/z (preferred precursor value)
ACC_ISOWIN_TARGET = "MS:1000827"  # isolation window target m/z (fallback)
UNIT_MINUTE = "UO:0000031"

_SCAN_NUM_RE = re.compile(r"scan=(\d+)")

PER_SCAN_COLUMNS = [
    "ScanNumber",
    "RetentionTimeSec",
    "PrecursorMZ",
    "IonInjectionTimeMs",
    "Intensity",
    "after_ms1",
]


def _strip_ns(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_ms2_headers(path) -> pd.DataFrame:
    """Stream an mzML file; return one row per MS2 spectrum (no peak decoding).

    ``after_ms1`` flags MS2 scans whose immediately preceding spectrum was MS1.
    Precursor m/z prefers "selected ion m/z" and falls back to the isolation
    window target m/z when the former is absent.
    """
    rows = []
    prev_ms_level = None
    ms_level = rt_sec = inj_ms = tic = sel_mz = iso_mz = None
    in_spectrum = False
    scan_num = None

    context = ET.iterparse(str(path), events=("start", "end"))
    _, root = next(context)
    for event, elem in context:
        tag = _strip_ns(elem.tag)
        if event == "start":
            if tag == "spectrum":
                in_spectrum = True
                ms_level = rt_sec = inj_ms = tic = sel_mz = iso_mz = None
                m = _SCAN_NUM_RE.search(elem.get("id", ""))
                scan_num = int(m.group(1)) if m else elem.get("index")
            continue

        # event == "end"
        if not in_spectrum:
            if tag == "spectrum":
                root.clear()
            continue

        if tag == "cvParam":
            acc = elem.get("accession")
            if acc == ACC_MS_LEVEL:
                ms_level = int(elem.get("value"))
            elif acc == ACC_SCAN_START:
                val = float(elem.get("value"))
                rt_sec = val * 60.0 if elem.get("unitAccession") == UNIT_MINUTE else val
            elif acc == ACC_INJ_TIME:
                inj_ms = float(elem.get("value"))
            elif acc == ACC_TIC:
                tic = float(elem.get("value"))
            elif acc == ACC_SELECTED_MZ:
                sel_mz = float(elem.get("value"))
            elif acc == ACC_ISOWIN_TARGET:
                iso_mz = float(elem.get("value"))
        elif tag == "spectrum":
            if ms_level == 2:
                prec_mz = sel_mz if sel_mz is not None else iso_mz
                rows.append(
                    (scan_num, rt_sec, prec_mz, inj_ms, tic, prev_ms_level == 1)
                )
            prev_ms_level = ms_level
            in_spectrum = False
            elem.clear()
            root.clear()

    return pd.DataFrame(rows, columns=PER_SCAN_COLUMNS)


def compute_metrics(ms2: pd.DataFrame) -> dict:
    """File-level acquisition-rate / injection-time / cycle-time metrics.

    Mirrors the R ``process_mzML_metadata``.
    """
    filt = ms2.loc[~ms2["after_ms1"]]  # exclude first MS2 after each MS1
    n_filt = len(filt)
    last_rt = filt["RetentionTimeSec"].iloc[-1] if n_filt else np.nan
    scan_rate_hz = n_filt / last_rt if (n_filt and last_rt and last_rt > 0) else np.nan

    cyc = ms2.loc[ms2["PrecursorMZ"].notna()].sort_values("RetentionTimeSec")
    pmz = cyc["PrecursorMZ"].to_numpy()
    rt = cyc["RetentionTimeSec"].to_numpy()
    reset_points = np.where(np.diff(pmz) < 0)[0] + 1  # negative step = new cycle
    cycle_times = np.diff(rt[reset_points])
    n_cycles = len(cycle_times)

    return dict(
        n_ms2_all=len(ms2),
        n_ms2_filtered=n_filt,
        rt_span_sec=ms2["RetentionTimeSec"].max() if len(ms2) else np.nan,
        ScanRateHz=scan_rate_hz,
        MeanInjectionTimeMs=filt["IonInjectionTimeMs"].mean(),
        MedianInjectionTimeMs=filt["IonInjectionTimeMs"].median(),
        MeanCycleTimeSec=float(np.nanmean(cycle_times)) if n_cycles else np.nan,
        MedianCycleTimeSec=float(np.nanmedian(cycle_times)) if n_cycles else np.nan,
        NumCycles=n_cycles,
        ScansPerCycle=(n_filt / n_cycles) if n_cycles else np.nan,
    )


def process_file(path, per_scan: bool = False):
    """Return (summary_row, per_scan_df|None) for one mzML. Picklable for a pool."""
    path = Path(path)
    ms2 = parse_ms2_headers(path)
    metrics = compute_metrics(ms2)
    summary = dict(folder=path.parent.name, filename=path.name, **metrics)

    per_scan_df = None
    if per_scan:
        per_scan_df = ms2.copy()
        per_scan_df["Ions"] = (
            per_scan_df["IonInjectionTimeMs"] * per_scan_df["Intensity"] / 1000.0
        )
        per_scan_df.insert(0, "filename", path.name)
        per_scan_df["ScanRateHz"] = metrics["ScanRateHz"]
        per_scan_df["MeanInjectionTimeMs"] = metrics["MeanInjectionTimeMs"]
    return summary, per_scan_df


def _process_file_star(args):
    return process_file(*args)


def add_filename_metadata(df: pd.DataFrame, pattern: str) -> pd.DataFrame:
    """Add columns from a regex with named groups, applied to ``filename``."""
    rx = re.compile(pattern)

    def extract(name):
        m = rx.search(name)
        return m.groupdict() if m else {}

    meta = df["filename"].map(extract).apply(pd.Series)
    return pd.concat([df, meta], axis=1)


def find_mzml(paths, recursive: bool = True, exclude_dirs=()):
    """Expand a mix of directories and files into a sorted list of mzML paths."""
    exclude = set(exclude_dirs)
    out = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            it = p.rglob("*.mzML") if recursive else p.glob("*.mzML")
            out.extend(f for f in it if f.parent.name not in exclude)
        elif p.suffix.lower() == ".mzml":
            out.append(p)
    # de-dup, keep sorted
    return sorted(set(out))


def collect(paths, recursive=True, exclude_dirs=(), workers=None,
            per_scan=False, regex=None, verbose=True):
    """Extract metrics from many mzML files in parallel.

    Returns (summary_df, per_scan_df|None). ``paths`` may be dirs and/or files.
    """
    files = find_mzml(
        [paths] if isinstance(paths, (str, Path)) else paths,
        recursive=recursive, exclude_dirs=exclude_dirs,
    )
    if not files:
        raise FileNotFoundError("No .mzML files found")

    summaries, per_scans = [], []
    workers = workers or min(20, len(files))
    if workers == 1:
        for i, f in enumerate(files, 1):
            summ, ps = process_file(f, per_scan)
            summaries.append(summ)
            if ps is not None:
                per_scans.append(ps)
            if verbose:
                print(f"[{i}/{len(files)}] {f.name}: {summ['ScanRateHz']:.2f} Hz, "
                      f"IT {summ['MeanInjectionTimeMs']:.2f} ms", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_process_file_star, (str(f), per_scan)): f for f in files}
            for i, fut in enumerate(as_completed(futs), 1):
                f = futs[fut]
                summ, ps = fut.result()
                summaries.append(summ)
                if ps is not None:
                    per_scans.append(ps)
                if verbose:
                    print(f"[{i}/{len(files)}] {f.name}: {summ['ScanRateHz']:.2f} Hz, "
                          f"IT {summ['MeanInjectionTimeMs']:.2f} ms", flush=True)

    summary_df = pd.DataFrame(summaries).sort_values(["folder", "filename"]).reset_index(drop=True)
    if regex:
        summary_df = add_filename_metadata(summary_df, regex)
    per_scan_df = pd.concat(per_scans, ignore_index=True) if per_scans else None
    return summary_df, per_scan_df


def main():
    ap = argparse.ArgumentParser(
        description="Extract injection time / acquisition rate / cycle time from mzML.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("paths", nargs="+", help="mzML file(s) and/or directories")
    ap.add_argument("--out", required=True, help="Per-file summary CSV output path")
    ap.add_argument("--per-scan-out", default=None,
                    help="Optional per-MS2-scan output (.parquet or .csv)")
    ap.add_argument("--regex", default=None,
                    help="Regex with named groups to parse metadata from file names")
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="Folder names to skip when scanning directories")
    ap.add_argument("--no-recursive", action="store_true",
                    help="Do not descend into subdirectories")
    ap.add_argument("--workers", type=int, default=None,
                    help="Parallel worker processes (default: min(20, n_files))")
    args = ap.parse_args()

    summary_df, per_scan_df = collect(
        args.paths,
        recursive=not args.no_recursive,
        exclude_dirs=args.exclude,
        workers=args.workers,
        per_scan=args.per_scan_out is not None,
        regex=args.regex,
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(out, index=False)
    print(f"\nWrote {out} ({len(summary_df)} files)")

    if per_scan_df is not None:
        ps_out = Path(args.per_scan_out)
        ps_out.parent.mkdir(parents=True, exist_ok=True)
        if ps_out.suffix.lower() == ".parquet":
            per_scan_df.to_parquet(ps_out, index=False)
        else:
            per_scan_df.to_csv(ps_out, index=False)
        print(f"Wrote {ps_out} ({len(per_scan_df):,} MS2 scans)")


if __name__ == "__main__":
    main()
