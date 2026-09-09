"""dark-tic-cohort - aggregate the dark-TIC screen across every run in a PRISM/Skyline
report to pinpoint cohort-wide trends and gaps in coverage.

Point it at ONE report (ideally a PRISM/Skyline transition report like merged_data.parquet,
which lists every replicate) and a directory of raw/mzML files. For each replicate that has
a matching file it runs the fast COARSE pass (scan headers only: which MS2 scans have an
eluting identification in their window, RT-gated), then aggregates:

  1. Per-run dark fraction (% of MS2 ion current in scans with no ID) - ranked, so QC
     outliers stand out.
  2. A COHORT-AVERAGED retention-time x isolation-window darkness map: cells dark across
     the whole cohort are systematic coverage gaps; cells dark in only some runs are
     sample-specific.
  3. A consistency map: in what fraction of runs each cell is majority-dark.

Writes a self-contained HTML dashboard + a per-run parquet. The coarse pass needs only
scan headers, so it is cheap and parallelizable; the fine/novel per-scan tier (peak decode)
is not run here (it is ~10 min/file) - use dark_tic_dashboard.py on individual runs a
cohort map flags as interesting.

    python dark_tic_cohort.py --report merged_data.parquet --raws /dir/of/raws --out cohort.html
    python dark_tic_cohort.py --report merged_data.parquet --raws "/dir/*.mzML" --workers 6

MIT licensed.
"""
from __future__ import annotations
import os, sys, glob, argparse, io, base64, importlib.util, time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("dts", os.path.join(_HERE, "dark_tic_screen.py"))
E = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(E)   # engine

RT_BIN = 0.5   # min


def window_intervals(ids, wt):
    """Per window: array of (start,end) elution intervals of its IDs (no fragment calc)."""
    lo = np.array([w["lo_mz"] for w in wt]); hi = np.array([w["hi_mz"] for w in wt])
    widx = np.array([w["win_idx"] for w in wt])
    by = defaultdict(list)
    for d in ids:
        sel = np.where((d["pmz"] >= lo) & (d["pmz"] < hi))[0]
        if sel.size:
            by[int(widx[sel[0]])].append((d["start"], d["end"]))
    return {w: np.array(v) for w, v in by.items()}


def coarse_map(rawpath, ids, rt_bin=RT_BIN):
    """Header-only pass: per-run dark fraction + a (rt_bin, win_idx) darkness grid."""
    from thermo_raw.thermo_raw_reader import assign_window_index
    sys.path.insert(0, os.path.join(_HERE, ".."))
    with E.open_run(rawpath) as rf:
        wt = rf.window_table()
        valid = {w["win_idx"] for w in wt}
        iv = window_intervals(ids, wt)
        ms2 = dark = unassigned = 0.0; n_ms2 = n_un = 0
        dmap = defaultdict(lambda: [0.0, 0.0])
        for sn in range(rf.first_scan, rf.last_scan + 1):
            h = rf.header(sn, with_inj=False)
            if h.ms_level != 2:
                continue
            w = assign_window_index(h.center_mz, wt)
            ms2 += h.tic; n_ms2 += 1
            if w not in valid:
                unassigned += h.tic; n_un += 1; continue
            ints = iv.get(w)
            lit = bool(ints is not None and np.any((ints[:, 0] <= h.rt_min) & (ints[:, 1] >= h.rt_min)))
            if not lit:
                dark += h.tic
            rb = int(h.rt_min / rt_bin); c = dmap[(rb, w)]; c[0] += h.tic
            if not lit:
                c[1] += h.tic
        centers = {w["win_idx"]: w["center_mz"] for w in wt}
    return dict(dark_frac=dark / ms2 if ms2 else float("nan"),
                un_frac=unassigned / ms2 if ms2 else float("nan"),
                n_ms2=n_ms2, ms2_tic=ms2, n_ids=len(ids), centers=centers,
                dmap={k: v for k, v in dmap.items()})


def _worker(report, raw, replicate, cols, rt_bin):
    try:
        if E.is_transition_level(report):
            c = dict(cols); c["replicate"] = replicate
            ids = E.report_from_transitions(report, replicate, c).get(replicate, [])
        else:
            ids = E.ids_for_raw(replicate, E.load_report(report, cols))
        if not ids:
            return replicate, {"error": "no IDs matched"}
        r = coarse_map(raw, ids, rt_bin)
        r["raw"] = os.path.basename(raw)
        return replicate, r
    except Exception as e:
        return replicate, {"error": repr(e)}


def match_files(reps, raws):
    """Pair each replicate name to a raw file (bidirectional containment)."""
    stems = {os.path.splitext(os.path.basename(p))[0]: p for p in raws}
    pairs = []
    for rep in reps:
        hit = None
        if rep in stems:
            hit = stems[rep]
        else:
            for st, p in stems.items():
                if rep in st or st in rep:
                    hit = p; break
        if hit:
            pairs.append((rep, hit))
    return pairs


def _png(fig):
    import matplotlib.pyplot as plt
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def figures(runs, cohort, centers):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ORANGE, INK, MUT = "#C65A1E", "#191723", "#64607A"
    out = {}
    stems = [r["rep"] for r in runs]; df = np.array([r["dark_frac"] * 100 for r in runs])

    # per-run dark fraction, ranked
    order = np.argsort(df); n = len(df)
    fig, ax = plt.subplots(figsize=(6.6, max(2.2, n * 0.12)))
    ax.barh(range(n), df[order], color=ORANGE, height=0.8)
    ax.axvline(np.median(df), color=MUT, ls="--", lw=1, label=f"median {np.median(df):.0f}%")
    ax.set_yticks([]); ax.set_xlabel("dark % of MS2 ion current"); ax.set_ylabel(f"{n} runs (ranked)")
    ax.legend(frameon=False, fontsize=8);
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.set_title("Per-run dark fraction", fontsize=11)
    out["perrun"] = _png(fig)

    # cohort-mean darkness map (rt_bin x window)
    wins = sorted(centers); wi = {w: k for k, w in enumerate(wins)}
    rbs = sorted({rb for rb, w in cohort})
    if rbs and wins:
        ri = {rb: k for k, rb in enumerate(rbs)}
        mean = np.full((len(wins), len(rbs)), np.nan)
        cons = np.full((len(wins), len(rbs)), np.nan)
        for (rb, w), (tic, dk, nr, ndarkrun) in cohort.items():
            if w in wi and tic > 0:
                mean[wi[w], ri[rb]] = dk / tic
                cons[wi[w], ri[rb]] = ndarkrun / nr if nr else np.nan
        ext = [min(rbs) * RT_BIN, max(rbs) * RT_BIN, centers[wins[0]], centers[wins[-1]]]
        for key, grid, lab in (("map", mean, "cohort-mean dark fraction"),
                               ("consistency", cons, "fraction of runs cell is majority-dark")):
            fig, ax = plt.subplots(figsize=(7.2, 3.4))
            im = ax.imshow(grid, aspect="auto", origin="lower", cmap="magma_r", vmin=0, vmax=1, extent=ext)
            ax.set_xlabel("retention time (min)"); ax.set_ylabel("isolation-window m/z")
            cb = fig.colorbar(im, ax=ax); cb.set_label(lab)
            out[key] = _png(fig)
    return out


def html(runs, figs, report, n_reps, n_matched):
    df = np.array([r["dark_frac"] * 100 for r in runs])
    un = np.array([r["un_frac"] * 100 for r in runs])
    ranked = sorted(runs, key=lambda r: -r["dark_frac"])
    def row(r):
        return (f"<tr><td>{r['rep'][:60]}</td><td>{r['dark_frac']*100:.1f}%</td>"
                f"<td>{r['un_frac']*100:.1f}%</td><td>{r['n_ids']:,}</td><td>{r['n_ms2']:,}</td></tr>")
    top = "".join(row(r) for r in ranked[:12])
    bot = "".join(row(r) for r in ranked[-6:]) if len(ranked) > 18 else ""
    tiles = [("runs aggregated", f"{len(runs)} / {n_reps}"),
             ("median dark", f"{np.median(df):.1f}%"),
             ("dark range", f"{df.min():.0f}-{df.max():.0f}%"),
             ("median unassigned", f"{np.median(un):.1f}%")]
    tilehtml = "".join(f'<div class=tile><div class=big>{v}</div><div class=lab>{k}</div></div>' for k, v in tiles)
    cons = f'<h2>3 · Darkness consistency</h2><p class=cap>Fraction of runs in which each RT x window cell is majority-dark. Cells near 1.0 (bright) are dark in <b>every</b> run - systematic coverage gaps, not sample noise.</p><img src="{figs["consistency"]}">' if "consistency" in figs else ""
    mapimg = f'<h2>2 · Cohort-mean darkness map</h2><p class=cap>Mean dark fraction of scan TIC across retention time (x) and isolation-window m/z (y), averaged over all runs. Bright bands are where the cohort is consistently unassigned.</p><img src="{figs["map"]}">' if "map" in figs else ""
    return f"""<!doctype html><meta charset=utf8><title>Dark-TIC cohort</title>
<style>
body{{font:15px/1.55 -apple-system,Segoe UI,sans-serif;color:#191723;background:#FAF9FC;margin:0}}
.wrap{{max-width:900px;margin:0 auto;padding:32px 24px 80px}}
h1{{font-size:24px;margin:0 0 2px}} .sub{{color:#64607A;margin:0 0 18px;font-family:monospace;font-size:12px;word-break:break-all}}
h2{{font-size:16px;margin:32px 0 4px}} .cap{{color:#64607A;font-size:13px;margin:0 0 10px}}
.tiles{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:14px 0}}
.tile{{background:#fff;border:1px solid #E6E3EE;border-radius:10px;padding:12px 14px}}
.big{{font-size:22px;font-weight:700}} .lab{{color:#64607A;font-size:11px;text-transform:uppercase;letter-spacing:.04em;margin-top:4px}}
img{{width:100%;border:1px solid #E6E3EE;border-radius:8px;background:#fff}}
table{{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}} th,td{{text-align:left;padding:5px 8px;border-bottom:1px solid #E6E3EE}}
th{{color:#64607A;font-size:11px;text-transform:uppercase}} td:nth-child(n+2){{font-variant-numeric:tabular-nums}}
.note{{color:#64607A;font-size:12px;margin-top:26px;border-top:1px solid #E6E3EE;padding-top:12px}}
</style>
<div class=wrap>
<h1>Dark-TIC cohort</h1><p class=sub>{os.path.basename(report)} &middot; {n_matched}/{n_reps} replicates matched to files</p>
<div class=tiles>{tilehtml}</div>
<h2>1 &middot; Per-run dark fraction</h2>
<p class=cap>Dark % of MS2 ion current per run (scans with no eluting ID in their window), ranked. Outliers at the high end are QC flags.</p>
<img src="{figs['perrun']}">
{mapimg}
{cons}
<h2>4 &middot; Darkest and cleanest runs</h2>
<table><tr><th>run</th><th>dark</th><th>unassigned</th><th>IDs</th><th>MS2 scans</th></tr>{top}{('<tr><td colspan=5 style=color:#B8B4C4>...</td></tr>'+bot) if bot else ''}</table>
<p class=note>dark-tic-cohort (PostdocToolbox). Coarse tier only: dark = MS2 scans with no confident ID (q&le;threshold) eluting in their isolation window, RT-gated to Skyline peak boundaries. Header-only, no peak decode. For per-scan intensity dark + unexplained features on a run of interest, use dark_tic_dashboard.py.</p>
</div>"""


def main():
    ap = argparse.ArgumentParser(description="Aggregate the dark-TIC coarse screen across a cohort.")
    ap.add_argument("--report", required=True, help="report listing all replicates (PRISM merged_data.parquet ideal)")
    ap.add_argument("--raws", required=True, help="directory of .raw/.mzML files (or a glob)")
    ap.add_argument("--out", default="dark_tic_cohort.html")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0, help="cap number of runs (0 = all; for a quick test)")
    ap.add_argument("--qvalue-max", type=float, default=0.01)
    ap.add_argument("--no-open", action="store_true")
    for c in ["seq", "pmz", "pch", "rep", "start", "end"]:
        ap.add_argument(f"--col-{c}", default=None)
    a = ap.parse_args()
    cols = dict(seq=a.col_seq, pmz=a.col_pmz, pch=a.col_pch, rep=a.col_rep,
                start=a.col_start, end=a.col_end, all_ids=False, no_rt=False,
                qmax=a.qvalue_max, replicate=None)
    raws = sorted(glob.glob(a.raws if any(c in a.raws for c in "*?[")
                  else os.path.join(a.raws, "*.raw")) + (
                  [] if any(c in a.raws for c in "*?[") else glob.glob(os.path.join(a.raws, "*.mzML"))))
    if E.is_transition_level(a.report):
        reps = E.distinct_replicates(a.report, limit=100000)
    else:
        reps = list(E.load_report(a.report, cols))
    pairs = match_files(reps, raws)
    if a.limit:
        pairs = pairs[:a.limit]
    print(f"[cohort] {len(reps)} replicates in report; {len(raws)} files; {len(pairs)} matched -> screening")
    if not pairs:
        raise SystemExit("[cohort] no replicate matched a file. Check --raws and names.")

    runs = []; cohort = defaultdict(lambda: [0.0, 0.0, 0, 0]); centers = {}
    t0 = time.time(); done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_worker, a.report, raw, rep, cols, RT_BIN): rep for rep, raw in pairs}
        for fut in as_completed(futs):
            rep, r = fut.result(); done += 1
            if "error" in r:
                print(f"  [skip] {rep}: {r['error']}"); continue
            r["rep"] = rep; runs.append(r)
            centers.update(r["centers"])
            for (rb, w), (tic, dk) in r["dmap"].items():          # merge into cohort accumulator
                cell = cohort[(rb, w)]
                cell[0] += tic; cell[1] += dk; cell[2] += 1
                if tic > 0 and dk / tic >= 0.5:
                    cell[3] += 1
            if done % 10 == 0 or done == len(pairs):
                print(f"  {done}/{len(pairs)} ({time.time()-t0:.0f}s)")
    if not runs:
        raise SystemExit("[cohort] no runs succeeded.")
    figs = figures(runs, cohort, centers)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(html(runs, figs, a.report, len(reps), len(pairs)))
    print(f"[cohort] {len(runs)} runs aggregated in {(time.time()-t0)/60:.1f} min -> {a.out}")
    try:
        import pyarrow as pa, pyarrow.parquet as pq
        pq.write_table(pa.table({"run": [r["rep"] for r in runs],
                                 "dark_frac": [r["dark_frac"] for r in runs],
                                 "unassigned_frac": [r["un_frac"] for r in runs],
                                 "n_ids": [r["n_ids"] for r in runs],
                                 "n_ms2": [r["n_ms2"] for r in runs]}),
                       os.path.splitext(a.out)[0] + "_perrun.parquet")
    except Exception:
        pass
    if not a.no_open:
        try:
            import webbrowser; webbrowser.open("file://" + os.path.abspath(a.out))
        except Exception:
            pass


if __name__ == "__main__":
    main()
