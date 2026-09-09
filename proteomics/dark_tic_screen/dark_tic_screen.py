"""dark-tic-screen - how much of a DIA run's ion current is NOT assigned to identified
peptides, is any of it high-intensity, and does it look like a proteoform?

A search only assigns signal to peptides in its space. This screens any narrow-window
DIA run (e.g. Orbitrap Astral) for the ion current the search leaves behind, as a
practical trigger for proteoform / PTM / novel-species follow-up:

  1. COARSE  - fraction of MS2 ion current in scans with no eluting identification
               (needs only the ID list + raw scan headers; fast).
  2. FINE    - after subtracting each eluting identified peptide's own fragment ladder
               (mono b/y z1-2 + 13C isotopes + H2O/NH3 losses + surviving precursor),
               how concentrated is the residual "dark" - a few big features or a haze?
  3. NOVEL   - the highest-intensity dark features that are NOT an identified peptide's
               fragment at another retention time (spillover), annotated with the delta
               mass to the nearest identified fragment for a panel of modifications, with
               decoy deltas so a hit must beat fragment-density coincidence.

A run that flags high concentration + genuinely-unexplained features at real mod offsets
is worth a proteoform / open search. A diffuse, spillover-dominated dark is not.

Inputs: a Skyline / DIA-NN report (parquet or csv) of confident identifications, and the
Thermo .raw files. Reuses the sibling PostdocToolbox `thermo_raw` reader.

    python dark_tic_screen.py --report ids.parquet --raws /path/to/raws --out screen/
    python dark_tic_screen.py --report ids.parquet --raws /raws --detail SAMPLE_STEM

MIT licensed. Lauren Fields / PostdocToolbox.
"""
from __future__ import annotations
import os, sys, glob, argparse, importlib
from collections import defaultdict
import numpy as np

# --- reuse the sibling thermo_raw reader -------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from thermo_raw.thermo_raw_reader import RawFile, assign_window_index  # noqa: E402

PROTON = 1.0072764669
WATER = 18.0105646863
NH3 = 17.0265491015
C13 = 1.0033548378
AA = {  # monoisotopic residue masses
    "G": 57.02146372, "A": 71.03711379, "S": 87.03202841, "P": 97.05276385,
    "V": 99.06841392, "T": 101.04767847, "C": 103.00918448, "L": 113.08406398,
    "I": 113.08406398, "N": 114.04292745, "D": 115.02694302, "Q": 128.05857751,
    "K": 128.09496302, "E": 129.04259309, "M": 131.04048509, "H": 137.05891186,
    "F": 147.06841392, "R": 156.10111102, "Y": 163.06332854, "W": 186.07931295,
    "U": 150.95363500, "O": 237.14772686,
}
UNIMOD = {1: 42.01056468, 4: 57.02146372, 5: 43.00581368, 7: 0.98401558,
          21: 79.96633052, 34: 14.01565006, 35: 15.99491462, 36: 28.03130013,
          28: -17.02654910, 27: -18.01056468, 385: -17.02654910}
# common modification neutral deltas for the NOVEL delta-mass panel
MODS = {"deamid +0.98": 0.98402, "methyl +14.02": 14.01565, "oxid +15.99": 15.99491,
        "water_loss -18.01": -18.01056, "dioxid +31.99": 31.98983,
        "acetyl +42.01": 42.01057, "trimethyl +42.05": 42.04695, "phospho +79.97": 79.96633}
DECOYS = [11.531, 23.273, 37.711, 6.191, 61.442, 88.131, 51.007, 70.884]
# snap a rounded bracket mass (Skyline display, e.g. C[+57]) to the exact mod mass
_SNAP = {57: 57.02146372, 16: 15.99491462, 80: 79.96633052, 42: 42.01056468,
         1: 0.98401558, 14: 14.01565006, 43: 43.00581368, 32: 31.98982924,
         -17: -17.02654910, -18: -18.01056468, 28: 28.03130013, 100: 100.01604399}


# --- modified-sequence parsing (DIA-NN/Skyline UnimodIds and bracket-mass) ----------
def parse_residue_masses(mod_seq: str):
    """Return (residue_masses, nterm_mod, cterm_mod). Supports 'C(UniMod:4)',
    'C(unimod:4)', bracket mass 'M[+15.9949]' / '[+42]-AA', and plain sequence."""
    s = mod_seq.strip()
    res, mods = [], []          # per-residue base mass, per-residue mod mass
    nterm = cterm = 0.0
    i, n = 0, len(s)

    def read_mod(j):
        """Parse a (…) or […] group starting at s[j]=='('/'['; return (mass, next_j)."""
        close = ")" if s[j] == "(" else "]"
        k = s.index(close, j)
        body = s[j + 1:k]
        low = body.lower()
        if low.startswith("unimod:"):
            m = UNIMOD.get(int(body.split(":")[1]), 0.0)
        else:
            try:
                m = float(body.replace("+", ""))
                r = round(m)                      # Skyline often shows rounded mass, e.g.
                if abs(m - r) < 0.5 and r in _SNAP:  # C[+57] -> snap to exact 57.02146
                    m = _SNAP[r]
            except ValueError:
                m = 0.0
        return m, k + 1

    # leading N-term mod, e.g. "(UniMod:1)" or "[+42]-"
    while i < n and s[i] in "([":
        m, i = read_mod(i)
        nterm += m
        if i < n and s[i] == "-":
            i += 1
    while i < n:
        c = s[i]
        if c in AA:
            res.append(AA[c]); mods.append(0.0); i += 1
            while i < n and s[i] in "([":
                m, i = read_mod(i)
                mods[-1] += m
        elif c in "-.":
            i += 1
        else:                    # unknown char: skip
            i += 1
    rm = np.array([r + m for r, m in zip(res, mods)], dtype=float)
    return rm, nterm, cterm


def explained_targets(mod_seq, precursor_mz=None, precursor_charge=None, max_z=2):
    """Sorted unique m/z an identified peptide legitimately owns: mono b/y (z1..max_z)
    + 13C(+1,+2) + H2O/NH3 losses + surviving precursor(+isotopes)."""
    rm, nterm, cterm = parse_residue_masses(mod_seq)
    n = rm.size
    if n < 2:
        return np.empty(0)
    pref = np.cumsum(rm)                       # b ladder neutral (before +proton)
    b_neu = pref[:-1] + nterm                  # b1..b(n-1)
    total = pref[-1] + nterm + cterm
    y_neu = (total + WATER) - b_neu            # complementary y
    out = []
    for z in range(1, max_z + 1):
        for neu in (b_neu, y_neu):
            base = (neu + z * PROTON) / z
            out.append(base)
            out.append(base + C13 / z); out.append(base + 2 * C13 / z)
            out.append((neu - WATER + z * PROTON) / z)
            out.append((neu - NH3 + z * PROTON) / z)
    if precursor_mz:
        pz = precursor_charge or 2
        out.append(np.array([precursor_mz, precursor_mz + C13 / pz,
                             precursor_mz + 2 * C13 / pz]))
    return np.unique(np.concatenate([np.atleast_1d(o) for o in out]))


def matched_mask(obs_mz, T, ppm):
    if T.size == 0 or obs_mz.size == 0:
        return np.zeros(obs_mz.size, bool)
    pos = np.searchsorted(T, obs_mz)
    left = np.clip(pos - 1, 0, T.size - 1); right = np.clip(pos, 0, T.size - 1)
    d = np.minimum(np.abs(T[left] - obs_mz), np.abs(T[right] - obs_mz))
    return d <= ppm * 1e-6 * obs_mz


# --- ID report loading --------------------------------------------------------------
def load_report(path, cols):
    """Return {replicate_stem: [ {mod_seq,pmz,pch,start,end}, ... ]}. Column names are
    auto-detected from common Skyline/DIA-NN headers; override with --col-*."""
    import pyarrow.parquet as pq
    if path.lower().endswith((".parquet", ".pq")):
        tbl = pq.read_table(path).to_pydict()
    else:
        import csv
        tbl = defaultdict(list)
        with open(path, newline="", encoding="utf-8-sig") as f:
            first = f.readline(); f.seek(0)
            delim = "\t" if first.count("\t") > first.count(",") else ","  # sniff CSV vs TSV
            for row in csv.DictReader(f, delimiter=delim):
                for k, v in row.items():
                    tbl[k].append(v)
    keys = list(tbl)

    def pick(cands, override):
        if override:
            return override
        for c in cands:
            for k in keys:
                if k.lower() == c.lower():
                    return k
        for c in cands:                        # substring fallback
            for k in keys:
                if c.lower() in k.lower():
                    return k
        raise SystemExit(f"[dark-tic] could not find a column among {cands}; "
                         f"available: {keys}\n  use the matching --col-* override.")
    c_seq = pick(["PeptideModifiedSequenceUnimodIds", "ModifiedPeptide",
                  "ModifiedSequence", "PeptideModifiedSequence", "Modified.Sequence"], cols["seq"])
    c_pmz = pick(["PrecursorMz", "Precursor.Mz", "PrecursorMzCalc"], cols["pmz"])
    c_pch = pick(["PrecursorCharge", "Charge", "Precursor.Charge"], cols["pch"])
    c_rep = pick(["ReplicateName", "Run", "FileName", "R.FileName"], cols["rep"])
    no_rt = cols.get("no_rt")
    c_st = c_en = None
    if not no_rt:
        c_st = pick(["StartTime", "MinStartTime", "RT.Start", "EG.StartRT"], cols["start"])
        c_en = pick(["EndTime", "MaxEndTime", "RT.Stop", "EG.EndRT"], cols["end"])
    c_q = None
    if not cols["all_ids"]:
        for c in ["DetectionQValue", "QValue", "Q.Value", "PEP"]:
            for k in keys:
                if k.lower() == c.lower():
                    c_q = k
    out = defaultdict(list)
    n = len(tbl[c_seq])
    for i in range(n):
        if c_q is not None:
            qv = tbl[c_q][i]
            if qv in (None, "", "NA", "#N/A") or (isinstance(qv, float) and np.isnan(qv)):
                continue                        # confident = q reported (Skyline emits q<=0.01)
        try:
            pmz = float(tbl[c_pmz][i]); pch = int(float(tbl[c_pch][i]))
            st = -1e9 if no_rt else float(tbl[c_st][i])
            en = 1e9 if no_rt else float(tbl[c_en][i])
        except (TypeError, ValueError):
            continue
        rep = str(tbl[c_rep][i])
        out[rep].append(dict(mod_seq=str(tbl[c_seq][i]), pmz=pmz, pch=pch, start=st, end=en))
    return out


def ids_for_raw(stem, report):
    """Match a raw stem to a replicate key (exact, then containment)."""
    if stem in report:
        return report[stem]
    for k in report:
        if k == stem or stem in k or k in stem:
            return report[k]
    return []


def build_window_ids(ids, wt):
    """Assign each ID to its DIA window; precompute its explained-fragment array."""
    lo = np.array([w["lo_mz"] for w in wt]); hi = np.array([w["hi_mz"] for w in wt])
    widx = np.array([w["win_idx"] for w in wt])
    by = defaultdict(list)
    for d in ids:
        sel = np.where((d["pmz"] >= lo) & (d["pmz"] < hi))[0]
        if not sel.size:
            continue
        frag = explained_targets(d["mod_seq"], d["pmz"], d["pch"])
        if frag.size:
            by[int(widx[sel[0]])].append(dict(start=d["start"], end=d["end"], frag=frag))
    return by


# --- the screen ---------------------------------------------------------------------
def screen_run(rawpath, ids, ppm, rt_lo, rt_hi, fine=False, top=200):
    with RawFile(rawpath) as rf:
        wt = rf.window_table()
        by = build_window_ids(ids, wt)
        win_full = {w: np.sort(np.unique(np.concatenate([e["frag"] for e in es])))
                    for w, es in by.items()} if fine else {}
        ms2_tic = dark_tic = 0.0
        n_ms2 = n_dark_scans = 0
        acc = defaultdict(lambda: defaultdict(lambda: defaultdict(float))) if fine else None
        dark_total = 0.0
        for sn in range(rf.first_scan, rf.last_scan + 1):
            h = rf.header(sn, with_inj=False)
            if h.ms_level != 2:
                continue
            if h.rt_min < rt_lo or h.rt_min > rt_hi:
                continue
            w = assign_window_index(h.center_mz, wt)
            es = by.get(w)
            active = [e for e in es if e["start"] <= h.rt_min <= e["end"]] if es else []
            ms2_tic += h.tic; n_ms2 += 1
            if not active:
                dark_tic += h.tic; n_dark_scans += 1
            if fine:
                mz, it = rf.peaks(sn)
                if mz.size == 0:
                    continue
                T = np.sort(np.unique(np.concatenate([e["frag"] for e in active]))) if active else np.empty(0)
                dark = ~matched_mask(mz, T, ppm) if T.size else np.ones(mz.size, bool)
                dmz, dit = mz[dark], it[dark]
                dark_total += float(dit.sum())
                rb = int(h.rt_min * 4)                    # 0.25-min bins
                for m, i in zip(dmz, dit):
                    acc[w][round(np.log(m) / np.log(1 + ppm * 1e-6))][rb] += float(i)
    res = dict(ms2_tic=ms2_tic, n_ms2=n_ms2, n_dark_scans=n_dark_scans,
               dark_scan_tic_frac=dark_tic / ms2_tic if ms2_tic else float("nan"),
               assigned_scan_tic_frac=1 - dark_tic / ms2_tic if ms2_tic else float("nan"))
    if fine:
        res["fine_dark_tic_frac"] = dark_total / ms2_tic if ms2_tic else float("nan")
        res.update(_fine(acc, win_full, dark_total, ppm, top))
    return res


def _fine(acc, win_full, dark_total, ppm, top):
    feats = []
    for w, wm in acc.items():
        for b, prof in wm.items():
            rbs = sorted(prof); run = []
            def flush(run):
                if run:
                    tot = sum(prof[r] for r in run)
                    feats.append((w, np.exp(b * np.log(1 + ppm * 1e-6)), tot,
                                  max(run, key=lambda r: prof[r]) / 4.0, len(run)))
            for r in rbs:
                if run and r != run[-1] + 1:
                    flush(run); run = []
                run.append(r)
            flush(run)
    feats.sort(key=lambda f: -f[2])
    inten = np.array([f[2] for f in feats]) if feats else np.array([0.0])
    cum = np.cumsum(inten) / inten.sum() if inten.sum() else inten
    half = int(np.searchsorted(cum, 0.5)) + 1 if feats else 0

    def infer_z(w, mz, rt):
        for z in (2, 3, 1):
            b = round(np.log(mz + 1.00335 / z) / np.log(1 + ppm * 1e-6))
            for bb in (b - 1, b, b + 1):
                p = acc[w].get(bb)
                if p and any(abs(r / 4.0 - rt) <= 0.5 for r in p):
                    return z
        return 0
    novel = []
    for w, mz, it, rt, nrt in feats[:top]:
        T = win_full.get(w, np.empty(0))
        if T.size and matched_mask(np.array([mz]), T, ppm)[0]:
            continue                                  # spillover of a real ID
        z = infer_z(w, mz, rt)
        rmh = [nm for nm, d in MODS.items()
               if z and (mz - d / z) > 0 and T.size and matched_mask(np.array([mz - d / z]), T, ppm)[0]]
        dch = sum(1 for d in DECOYS
                  if z and (mz - d / z) > 0 and T.size and matched_mask(np.array([mz - d / z]), T, ppm)[0])
        novel.append(dict(win=w, mz=mz, z=z, pct_dark=it / dark_total * 100 if dark_total else 0,
                          rt=rt, n_rt=nrt, mod_hits=rmh, n_decoy_hits=dch))
    return dict(n_features=len(feats), top1_pct=inten[0] / inten.sum() * 100 if inten.sum() else 0,
                top500_pct=cum[min(500, len(cum)) - 1] * 100 if feats else 0,
                feats_to_50pct=half, novel=novel)


def main():
    ap = argparse.ArgumentParser(description="Screen DIA runs for high-intensity unassigned (dark) ion current.")
    ap.add_argument("--report", required=True, help="Skyline/DIA-NN report (.parquet or .csv) of confident IDs")
    ap.add_argument("--raws", required=True, help="directory of .raw files (or a glob)")
    ap.add_argument("--out", default=None, help="output directory for the per-run summary")
    ap.add_argument("--detail", default=None, help="raw stem to run the FINE + novel screen on "
                    "(default: the run with the most MS2 ion current)")
    ap.add_argument("--all-fine", action="store_true", help="run the fine screen on every file (slow)")
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--rt-range", default="0,1e9", help="RT window in minutes, 'lo,hi'")
    ap.add_argument("--top", type=int, default=200, help="top-N features to inspect for novelty")
    ap.add_argument("--all-ids", action="store_true", help="treat every report row as confident (ignore q-value)")
    ap.add_argument("--no-rt-gate", action="store_true", help="report has no RT peak boundaries: "
                    "subtract each ID's fragments across ALL RT in its window (conservative lower-bound dark; "
                    "the coarse tier is not meaningful in this mode)")
    for c in ["seq", "pmz", "pch", "rep", "start", "end"]:
        ap.add_argument(f"--col-{c}", default=None)
    a = ap.parse_args()
    rt_lo, rt_hi = (float(x) for x in a.rt_range.split(","))
    cols = dict(seq=a.col_seq, pmz=a.col_pmz, pch=a.col_pch, rep=a.col_rep,
                start=a.col_start, end=a.col_end, all_ids=a.all_ids, no_rt=a.no_rt_gate)

    report = load_report(a.report, cols)
    raws = sorted(glob.glob(a.raws if any(c in a.raws for c in "*?[") else os.path.join(a.raws, "*.raw")))
    if not raws:
        raise SystemExit(f"[dark-tic] no .raw files under {a.raws}")
    print(f"[dark-tic] {len(report)} replicates in report; {len(raws)} raw files\n")

    print(f"[dark-tic] COARSE (assigned vs dark MS2 ion current per run):")
    print(f"  {'run':44} {'nID':>6} {'dark_TIC%':>9} {'dark_scans%':>11}")
    coarse = []
    for rp in raws:
        stem = os.path.splitext(os.path.basename(rp))[0]
        ids = ids_for_raw(stem, report)
        r = screen_run(rp, ids, a.ppm, rt_lo, rt_hi, fine=False)
        coarse.append((stem, len(ids), r))
        print(f"  {stem[:44]:44} {len(ids):>6} {r['dark_scan_tic_frac']*100:>8.1f} "
              f"{r['n_dark_scans']/max(r['n_ms2'],1)*100:>10.1f}")

    detail = a.detail or max(coarse, key=lambda c: c[2]["ms2_tic"])[0]
    targets = [c[0] for c in coarse] if a.all_fine else [detail]
    for stem in targets:
        rp = next(r for r in raws if os.path.splitext(os.path.basename(r))[0] == stem)
        ids = ids_for_raw(stem, report)
        print(f"\n[dark-tic] FINE + NOVEL screen: {stem}")
        r = screen_run(rp, ids, a.ppm, rt_lo, rt_hi, fine=True, top=a.top)
        print(f"  fine dark = {r['fine_dark_tic_frac']*100:.1f}% of MS2 TIC "
              f"(coarse scan-level {r['dark_scan_tic_frac']*100:.1f}%); {r['n_features']:,} dark features")
        print(f"  concentration: top1={r['top1_pct']:.2f}%  top500={r['top500_pct']:.1f}%  "
              f"features to 50% = {r['feats_to_50pct']:,}")
        nov = r["novel"]
        nmod = sum(1 for x in nov if x["mod_hits"])
        ndec = sum(1 for x in nov if x["n_decoy_hits"])
        print(f"  of top {a.top} features: {len(nov)} genuinely unexplained (not spillover)")
        print(f"  delta-mass: {nmod}/{len(nov)} match a real mod, {ndec}/{len(nov)} match a decoy delta "
              f"-> {'PROTEOFORM SIGNAL' if nmod > ndec * 1.5 else 'coincidence (no proteoform signal)'}")
        big = [x for x in nov if x["pct_dark"] >= 0.5]
        if big:
            print(f"  high-intensity unexplained features (>=0.5% of dark):")
            for x in big[:15]:
                print(f"    m/z {x['mz']:9.4f} z{x['z']} rt {x['rt']:5.1f}  {x['pct_dark']:5.2f}% "
                      f"mods:[{','.join(x['mod_hits']) or '-'}]")
        else:
            print("  no single unexplained feature reaches 0.5% of the dark -> diffuse, low-abundance.")
        flag = big and nmod > ndec * 1.5
        print(f"\n  VERDICT: {'FLAG - high-intensity unexplained proteoform-like signal; worth an open/delta-mass search.' if flag else 'NEGATIVE - dark is diffuse and/or spillover; no high-abundance novel species.'}")

    if a.out:
        import pyarrow as pa, pyarrow.parquet as pq
        os.makedirs(a.out, exist_ok=True)
        pq.write_table(pa.table({
            "run": [c[0] for c in coarse], "n_ids": [c[1] for c in coarse],
            "ms2_tic": [c[2]["ms2_tic"] for c in coarse],
            "dark_tic_frac": [c[2]["dark_scan_tic_frac"] for c in coarse],
            "dark_scan_frac": [c[2]["n_dark_scans"] / max(c[2]["n_ms2"], 1) for c in coarse]}),
            os.path.join(a.out, "coarse_summary.parquet"))
        print(f"\n[dark-tic] wrote {a.out}/coarse_summary.parquet")


if __name__ == "__main__":
    main()
