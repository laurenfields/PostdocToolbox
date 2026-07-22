"""Profile DIA cycle time across the gradient, straight from Thermo .raw files.

Reads scan headers only (peak arrays are never decoded) via the ProteoWizard-bundled
ThermoFisher.CommonCore RawFileReader DLLs, so no mzML conversion step is needed and a
600-file sweep takes minutes rather than hours.

Cycle time is the retention-time delta between consecutive cycle starts. Two detection
modes, chosen automatically per file:

  ms1     ΔRT between consecutive MS1 survey scans (default when MS1 scans are present)
  mzreset ΔRT between precursor-m/z resets, for MS2-only DIA methods with no survey scan

Unlike a per-file summary, this emits EVERY cycle with the retention time at which it
started (``--per-cycle-out``), which is what makes cycle time resolvable against the
gradient. That matters because cycle time is often not constant within a run: when few
ions elute (gradient edges) the injection time saturates and the cycle stretches, while
mid-gradient the AGC target fills fast and the cycle collapses. A single median hides a
swing that can exceed 200%.

See also ``proteomics/mzml_ion_timing`` -- same family of metrics computed from mzML.
Use that one if you already have mzML, need TIC/ion counts, or are on non-Thermo data.
Use this one to skip conversion, or when you need RT-resolved per-cycle values.

Vendor/instrument-agnostic beyond requiring Thermo .raw: results are keyed by file name
and no naming scheme is assumed. To pull experimental factors out of file names, pass a
regex with named groups via ``--regex`` -- each named group becomes a column.

CLI
---
    # one summary row per file
    python raw_cycle_time_profile.py RAW_DIR --out summary.csv

    # plus every individual cycle, for RT-resolved analysis
    python raw_cycle_time_profile.py RAW_DIR --out summary.csv \
        --per-cycle-out cycles.parquet

    # parse factors out of file names
    python raw_cycle_time_profile.py RAW_DIR --out summary.csv \
        --regex '(?P<window>\\d+to\\d+)_(?P<grad>\\d+m)_run(?P<run>\\d+)'

Dependencies
------------
pythonnet, pandas, numpy (pyarrow only if writing .parquet), and a local ProteoWizard
install providing ThermoFisher.CommonCore.Data.dll + .RawFileReader.dll.
"""
import argparse
import os
import re
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

# pythonnet must target .NET Framework for these DLLs on Windows.
os.environ.setdefault("PYTHONNET_RUNTIME", "netfx")

SUMMARY_COLUMNS = [
    "File", "NumScans", "NumMS1", "NumMS2", "ScansPerCycle", "CycleDetection",
    "MedianCycleTimeSec", "MeanCycleTimeSec", "SdCycleTimeSec",
    "P05CycleTimeSec", "P95CycleTimeSec", "MinCycleTimeSec", "MaxCycleTimeSec",
    "MedianMS2ScanTimeMs", "MedianMS2InjectionTimeMs", "InjectionTimeFraction",
    "RTStartMin", "RTEndMin", "MS2ScanRange", "NumCycles", "Error",
]

_PWIZ_GLOBS = [
    r"C:\Users\*\AppData\Local\Apps\ProteoWizard*",
    r"C:\Program Files\ProteoWizard\ProteoWizard*",
    r"C:\Program Files (x86)\ProteoWizard\ProteoWizard*",
]
_REQUIRED_DLLS = ("ThermoFisher.CommonCore.Data.dll",
                  "ThermoFisher.CommonCore.RawFileReader.dll")

_raw_adapter = None      # set per worker process by _init_clr
_Device = None


def find_pwiz(explicit=None):
    """Locate a ProteoWizard install containing the Thermo reader DLLs."""
    candidates = []
    if explicit:
        candidates = [Path(explicit)]
    else:
        if os.environ.get("PWIZ_DIR"):
            candidates.append(Path(os.environ["PWIZ_DIR"]))
        import glob
        for pat in _PWIZ_GLOBS:
            candidates += [Path(p) for p in glob.glob(pat)]
    for c in sorted(candidates, reverse=True):     # prefer the newest build
        if all((c / d).exists() for d in _REQUIRED_DLLS):
            return c
    raise SystemExit(
        "Could not find ProteoWizard's Thermo reader DLLs (%s).\n"
        "Pass --pwiz <dir> or set PWIZ_DIR to the ProteoWizard install directory."
        % ", ".join(_REQUIRED_DLLS))


def _init_clr(pwiz_dir):
    global _raw_adapter, _Device
    import clr
    for dll in _REQUIRED_DLLS:
        clr.AddReference(str(Path(pwiz_dir) / dll))
    from ThermoFisher.CommonCore.Data.Business import Device
    from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
    _raw_adapter, _Device = RawFileReaderAdapter, Device


def _cycle_starts_from_mzresets(raw, first, last):
    """Cycle boundaries for MS2-only methods: where precursor m/z steps back down."""
    starts, prev = [], None
    for s in range(first, last + 1):
        try:
            mz = raw.GetScanEventForScanNumber(s).GetReaction(0).PrecursorMass
        except Exception:
            continue
        if prev is None or mz < prev:
            starts.append(s)
        prev = mz
    return starts


def analyze(path, want_cycles=False):
    """Return (summary dict, cycles DataFrame or None) for one raw file."""
    rec = {c: None for c in SUMMARY_COLUMNS}
    rec["File"] = os.path.basename(path)
    raw = None
    try:
        raw = _raw_adapter.FileFactory(str(path))
        if raw.IsError:
            rec["Error"] = "raw file reports error state"
            return rec, None
        raw.SelectInstrument(_Device.MS, 1)
        hdr = raw.RunHeaderEx
        first, last = hdr.FirstSpectrum, hdr.LastSpectrum
        rec["NumScans"] = last - first + 1
        rec["RTStartMin"] = round(hdr.StartTime, 4)
        rec["RTEndMin"] = round(hdr.EndTime, 4)

        ms1 = list(raw.GetFilteredScanEnumerator(raw.GetFilterFromString("ms")))
        rec["NumMS1"] = len(ms1)
        rec["NumMS2"] = rec["NumScans"] - len(ms1)

        if len(ms1) >= 3:
            starts, rec["CycleDetection"] = ms1, "ms1"
        else:
            starts, rec["CycleDetection"] = _cycle_starts_from_mzresets(raw, first, last), "mzreset"
        if len(starts) < 3:
            rec["Error"] = "fewer than 3 detectable cycles"
            return rec, None

        gaps = [starts[i + 1] - starts[i] - 1 for i in range(len(starts) - 1)]
        rec["ScansPerCycle"] = float(statistics.median(gaps))

        rts = [raw.RetentionTimeFromScanNumber(s) for s in starts]
        cyc = np.array([(rts[i + 1] - rts[i]) * 60.0 for i in range(len(rts) - 1)])
        rec["NumCycles"] = int(cyc.size)
        rec["MedianCycleTimeSec"] = round(float(np.median(cyc)), 4)
        rec["MeanCycleTimeSec"] = round(float(cyc.mean()), 4)
        rec["SdCycleTimeSec"] = round(float(cyc.std(ddof=1)), 4) if cyc.size > 1 else 0.0
        rec["P05CycleTimeSec"] = round(float(np.percentile(cyc, 5)), 4)
        rec["P95CycleTimeSec"] = round(float(np.percentile(cyc, 95)), 4)
        rec["MinCycleTimeSec"] = round(float(cyc.min()), 4)
        rec["MaxCycleTimeSec"] = round(float(cyc.max()), 4)
        if rec["ScansPerCycle"]:
            rec["MedianMS2ScanTimeMs"] = round(
                rec["MedianCycleTimeSec"] * 1000.0 / (rec["ScansPerCycle"] + 1), 3)

        # Sample MS2 scans for injection time + scan range. Sampled, not exhaustive:
        # per-scan trailer access is the expensive call in this reader.
        sample, step = [], max(1, (last - first) // 200)
        for s in range(first, last + 1, step):
            try:
                f = raw.GetFilterForScanNumber(s)
                if str(f.MSOrder) != "Ms2":
                    continue
                if rec["MS2ScanRange"] is None:
                    rec["MS2ScanRange"] = "%.0f-%.0f" % (
                        f.GetMassRange(0).Low, f.GetMassRange(0).High)
                tr = raw.GetTrailerExtraInformation(s)
                for lbl, val in zip(tr.Labels, tr.Values):
                    if "Injection Time" in lbl:
                        sample.append(float(val))
                        break
            except Exception:
                continue
        if sample:
            rec["MedianMS2InjectionTimeMs"] = round(statistics.median(sample), 3)
            if rec["MedianMS2ScanTimeMs"]:
                rec["InjectionTimeFraction"] = round(
                    rec["MedianMS2InjectionTimeMs"] / rec["MedianMS2ScanTimeMs"], 4)

        cycles = None
        if want_cycles:
            cycles = pd.DataFrame({
                "File": rec["File"],
                "CycleIndex": np.arange(cyc.size, dtype=np.int32),
                "RTMin": np.round(rts[:-1], 5),
                "CycleTimeSec": np.round(cyc, 5),
            })
        return rec, cycles
    except Exception as e:
        rec["Error"] = "%s: %s" % (type(e).__name__, str(e)[:200].replace("\n", " "))
        return rec, None
    finally:
        if raw is not None:
            try:
                raw.Dispose()
            except Exception:
                pass


def _work(args):
    return analyze(*args)


def collect_paths(paths, recursive=True):
    out = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            out += sorted(p.rglob("*.raw") if recursive else p.glob("*.raw"))
        elif p.suffix.lower() == ".raw":
            out.append(p)
    seen, uniq = set(), []
    for p in out:
        if p.name not in seen:
            seen.add(p.name)
            uniq.append(p)
    return uniq


def apply_regex(df, pattern):
    rx = re.compile(pattern)
    rows = []
    for name in df["File"]:
        m = rx.search(name)
        rows.append(m.groupdict() if m else {})
    meta = pd.DataFrame(rows, index=df.index)
    for c in meta.columns:
        # numeric where the whole column converts; leave as strings otherwise.
        # (pandas 3.0 removed errors="ignore", so convert-or-revert explicitly)
        try:
            meta[c] = pd.to_numeric(meta[c])
        except (ValueError, TypeError):
            pass
    unmatched = int(meta.isna().all(axis=1).sum()) if len(meta.columns) else len(df)
    if unmatched:
        print("  warning: --regex matched no groups in %d/%d file names"
              % (unmatched, len(df)), file=sys.stderr)
    return pd.concat([df, meta], axis=1) if len(meta.columns) else df


def write_table(df, path):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.suffix.lower() == ".parquet":
        df.to_parquet(p, index=False)
    else:
        df.to_csv(p, index=False)


def main():
    ap = argparse.ArgumentParser(
        description="Profile DIA cycle time across the gradient from Thermo .raw files.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help=".raw file(s) and/or directories")
    ap.add_argument("--out", required=True, help="Per-file summary output (.csv/.parquet)")
    ap.add_argument("--per-cycle-out", default=None,
                    help="Optional per-cycle output with RT (.parquet or .csv)")
    ap.add_argument("--regex", default=None,
                    help="Regex with named groups to parse metadata from file names")
    ap.add_argument("--workers", type=int, default=8, help="Parallel worker processes")
    ap.add_argument("--pwiz", default=None,
                    help="ProteoWizard install dir (else auto-detect / $PWIZ_DIR)")
    ap.add_argument("--no-recursive", action="store_true",
                    help="Do not descend into subdirectories")
    args = ap.parse_args()

    pwiz = find_pwiz(args.pwiz)
    files = collect_paths(args.paths, recursive=not args.no_recursive)
    if not files:
        raise SystemExit("No .raw files found in: %s" % ", ".join(args.paths))
    print("%d raw files | ProteoWizard: %s | %d workers"
          % (len(files), pwiz, args.workers), flush=True)

    want = args.per_cycle_out is not None
    summaries, cycle_parts, failed = [], [], 0
    with ProcessPoolExecutor(max_workers=args.workers,
                             initializer=_init_clr, initargs=(pwiz,)) as ex:
        for i, (rec, cycles) in enumerate(
                ex.map(_work, [(str(f), want) for f in files], chunksize=2), 1):
            summaries.append(rec)
            if cycles is not None:
                cycle_parts.append(cycles)
            if rec["Error"]:
                failed += 1
                print("  ! %s -> %s" % (rec["File"], rec["Error"]), flush=True)
            if i % 50 == 0 or i == len(files):
                print("  %d/%d" % (i, len(files)), flush=True)

    df = pd.DataFrame(summaries, columns=SUMMARY_COLUMNS)
    if args.regex:
        df = apply_regex(df, args.regex)
    write_table(df, args.out)
    print("wrote %s (%d rows, %d failed)" % (args.out, len(df), failed))

    if want:
        if not cycle_parts:
            print("no cycles extracted; skipping --per-cycle-out", file=sys.stderr)
        else:
            cd = pd.concat(cycle_parts, ignore_index=True)
            write_table(cd, args.per_cycle_out)
            print("wrote %s (%d cycles across %d files)"
                  % (args.per_cycle_out, len(cd), cd["File"].nunique()))

    ok = df[df["Error"].isna()]
    if len(ok):
        print("median cycle time across files: %.3f s (range %.3f-%.3f)"
              % (ok["MedianCycleTimeSec"].median(),
                 ok["MedianCycleTimeSec"].min(), ok["MedianCycleTimeSec"].max()))


if __name__ == "__main__":
    main()
