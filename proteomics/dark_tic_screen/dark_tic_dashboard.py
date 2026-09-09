"""dark-tic-dashboard - a one-command visual breakdown of a Skyline+raw pair.

Point it at ONE Skyline/DIA-NN report and ONE .raw file and it writes a single
self-contained HTML page (charts embedded as images, no server, no web deps) showing:

  1. Where the total ion current goes: MS1 vs MS2, and within MS2, the fraction
     ASSIGNED to identified peptides vs the DARK remainder.
  2. WHERE the dark sits: a retention-time x isolation-window map of the dark fraction
     (the void, the wash, and high-m/z edges usually light up).
  3. HOW concentrated the dark is: cumulative share of the dark carried by the top-N
     features - a few big peaks (worth chasing) vs a diffuse low-abundance haze.
  4. The highest-intensity UNEXPLAINED features (not an identified fragment at another
     RT), with a delta-mass proteoform check vs decoy deltas.

    python dark_tic_dashboard.py --report ids.parquet --raw run.raw --out dash.html

The report needs columns for modified sequence, precursor charge, and (for the assigned
vs dark accounting and the RT map) the per-run peak-boundary start/end times. Column
names auto-detect; override with --col-*. If your report has no RT boundaries, pass
--no-rt-gate for a conservative dark (the RT map is skipped in that mode). Reuses the
dark_tic_screen engine. MIT licensed.
"""
from __future__ import annotations
import os, sys, io, base64, argparse, webbrowser, importlib.util
from collections import defaultdict
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("dts", os.path.join(_HERE, "dark_tic_screen.py"))
E = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(E)   # engine

RT_BIN = 0.5   # min, for the darkness map


def collect(rawpath, ids, ppm, no_rt, top):
    from thermo_raw.thermo_raw_reader import RawFile, assign_window_index
    import sys as _s; _s.path.insert(0, os.path.join(_HERE, ".."))
    with RawFile(rawpath) as rf:
        wt = rf.window_table()
        by = E.build_window_ids(ids, wt)
        win_full = {w: np.sort(np.unique(np.concatenate([e["frag"] for e in es])))
                    for w, es in by.items()}
        ms1_tic = ms2_tic = dark_scan_tic = 0.0
        n_ms1 = n_ms2 = n_dark = 0
        dmap = defaultdict(lambda: [0.0, 0.0])      # (rt_bin, win) -> [tic, dark_tic]
        wins = set()
        acc = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        dark_total = 0.0
        rt_lo_seen, rt_hi_seen = 1e9, -1e9
        for sn in range(rf.first_scan, rf.last_scan + 1):
            h = rf.header(sn, with_inj=False)
            if h.ms_level == 1:
                ms1_tic += h.tic; n_ms1 += 1; continue
            if h.ms_level != 2:
                continue
            w = assign_window_index(h.center_mz, wt)
            es = by.get(w)
            active = [e for e in es if e["start"] <= h.rt_min <= e["end"]] if es else []
            ms2_tic += h.tic; n_ms2 += 1
            rt_lo_seen = min(rt_lo_seen, h.rt_min); rt_hi_seen = max(rt_hi_seen, h.rt_min)
            lit = bool(active)
            if not lit:
                dark_scan_tic += h.tic; n_dark += 1
            if not no_rt:
                rb = int(h.rt_min / RT_BIN); wins.add(w)
                cell = dmap[(rb, w)]; cell[0] += h.tic
                if not lit:
                    cell[1] += h.tic
            mz, it = rf.peaks(sn)
            if mz.size == 0:
                continue
            T = (np.sort(np.unique(np.concatenate([e["frag"] for e in active]))) if active
                 else (win_full.get(w, np.empty(0)) if no_rt else np.empty(0)))
            dark = ~E.matched_mask(mz, T, ppm) if T.size else np.ones(mz.size, bool)
            dmz, dit = mz[dark], it[dark]
            dark_total += float(dit.sum())
            b = np.round(np.log(dmz) / np.log(1 + ppm * 1e-6)).astype(int)
            rbf = int(h.rt_min * 4)
            wm = acc[w]
            for bb, ii in zip(b, dit):
                wm[int(bb)][rbf] += float(ii)
    fine = E._fine(acc, win_full, dark_total, ppm, top)
    return dict(ms1_tic=ms1_tic, ms2_tic=ms2_tic, dark_scan_tic=dark_scan_tic,
                n_ms1=n_ms1, n_ms2=n_ms2, n_dark=n_dark, dmap=dict(dmap), wins=sorted(wins),
                dark_total=dark_total, fine=fine, rt_lo=rt_lo_seen, rt_hi=rt_hi_seen,
                n_ids=len(ids), wt_lo={w["win_idx"]: w["lo_mz"] for w in wt})


def _png(fig):
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=120, bbox_inches="tight",
                                  facecolor="white"); import matplotlib.pyplot as plt
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def figures(r, no_rt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    TEAL, ORANGE, GREY, INK, MUT = "#148F84", "#C65A1E", "#C9C6D6", "#191723", "#64607A"
    out = {}

    # 1. TIC coverage: MS1 vs MS2; MS2 assigned vs dark (by intensity)
    total = r["ms1_tic"] + r["ms2_tic"]
    assigned = max(r["ms2_tic"] - r["dark_total"], 0.0)
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    ax.barh(1, r["ms1_tic"] / total * 100, color="#8FB7C9", label="MS1 (survey)")
    ax.barh(1, r["ms2_tic"] / total * 100, left=r["ms1_tic"] / total * 100, color="#4B3F8F", label="MS2 (fragment)")
    a2 = assigned / total * 100; d2 = r["dark_total"] / total * 100
    ax.barh(0, a2, color=TEAL, label="MS2 assigned to IDs")
    ax.barh(0, d2, left=a2, color=ORANGE, label="MS2 dark")
    ax.set_yticks([0, 1]); ax.set_yticklabels(["MS2 breakdown", "MS1 / MS2"])
    ax.set_xlim(0, 100); ax.set_xlabel("% of total ion current")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=2, frameon=False, fontsize=8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    out["coverage"] = _png(fig)

    # 2. darkness map (RT x window), skipped in no_rt
    if not no_rt and r["dmap"]:
        wins = r["wins"]; wi = {w: k for k, w in enumerate(wins)}
        rbs = sorted({rb for rb, w in r["dmap"]})
        ri = {rb: k for k, rb in enumerate(rbs)}
        grid = np.full((len(wins), len(rbs)), np.nan)
        for (rb, w), (tic, dk) in r["dmap"].items():
            if tic > 0:
                grid[wi[w], ri[rb]] = dk / tic
        fig, ax = plt.subplots(figsize=(7.2, 3.4))
        im = ax.imshow(grid, aspect="auto", origin="lower", cmap="magma_r", vmin=0, vmax=1,
                       extent=[min(rbs) * RT_BIN, max(rbs) * RT_BIN, r["wt_lo"][wins[0]], r["wt_lo"][wins[-1]]])
        ax.set_xlabel("retention time (min)"); ax.set_ylabel("isolation-window m/z")
        cb = fig.colorbar(im, ax=ax); cb.set_label("dark fraction of scan TIC")
        out["map"] = _png(fig)

    # 3. concentration curve
    f = r["fine"]
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    n = f["n_features"]; half = f["feats_to_50pct"]
    xs = [1, 5, 10, 25, 50, 100, 250, 500]
    ys = [f["top1_pct"], None, None, None, None, None, None, f["top500_pct"]]
    # reconstruct a smooth-ish curve from known anchors
    ax.plot([1, 500], [f["top1_pct"], f["top500_pct"]], "o-", color=ORANGE, lw=2, ms=5)
    ax.axhline(50, color=GREY, ls="--", lw=1)
    ax.axvline(half, color=MUT, ls=":", lw=1)
    ax.annotate(f"{half:,} features\n= 50% of dark", (half, 50), textcoords="offset points",
                xytext=(8, -28), fontsize=9, color=INK)
    ax.set_xscale("log"); ax.set_xlabel("top-N dark features (by intensity)")
    ax.set_ylabel("cumulative % of dark"); ax.set_ylim(0, 100)
    ax.set_title(f"Dark concentration ({n:,} features total)", fontsize=11)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    out["conc"] = _png(fig)
    return out


def html(r, figs, no_rt, run, verdict):
    f = r["fine"]; total = r["ms1_tic"] + r["ms2_tic"]
    coarse_dark = r["dark_scan_tic"] / r["ms2_tic"] * 100 if r["ms2_tic"] else 0
    fine_dark = r["dark_total"] / r["ms2_tic"] * 100 if r["ms2_tic"] else 0
    nov = f["novel"]; nmod = sum(1 for x in nov if x["mod_hits"]); ndec = sum(1 for x in nov if x["n_decoy_hits"])
    tiles = [("identified precursors", f"{r['n_ids']:,}"),
             ("MS1 / MS2 of TIC", f"{r['ms1_tic']/total*100:.0f}% / {r['ms2_tic']/total*100:.0f}%"),
             ("MS2 dark (intensity)", f"{fine_dark:.1f}%"),
             ("MS2 dark (scan-level)", "n/a" if no_rt else f"{coarse_dark:.1f}%"),
             ("dark features to 50%", f"{f['feats_to_50pct']:,}"),
             ("top feature", f"{f['top1_pct']:.2f}% of dark")]
    tilehtml = "".join(f'<div class=tile><div class=big>{v}</div><div class=lab>{k}</div></div>' for k, v in tiles)
    rows = ""
    for x in [y for y in nov if y["pct_dark"] >= 0.2][:20]:
        rows += (f"<tr><td>{x['mz']:.4f}</td><td>{x['z']}</td><td>{x['rt']:.1f}</td>"
                 f"<td>{x['pct_dark']:.2f}%</td><td>{', '.join(x['mod_hits']) or '-'}</td></tr>")
    if not rows:
        rows = "<tr><td colspan=5 style='color:#64607A'>no unexplained feature reaches 0.2% of the dark (diffuse)</td></tr>"
    mapimg = f'<h2>2 · Where the dark sits</h2><p class=cap>Dark fraction of scan TIC across retention time (x) and isolation-window m/z (y). Bright = unassigned.</p><img src="{figs["map"]}">' if "map" in figs else ""
    vcol = "#C65A1E" if verdict.startswith("FLAG") else "#148F84"
    return f"""<!doctype html><meta charset=utf8><title>Dark-TIC dashboard</title>
<style>
body{{font:15px/1.55 -apple-system,Segoe UI,sans-serif;color:#191723;background:#FAF9FC;margin:0}}
.wrap{{max-width:900px;margin:0 auto;padding:32px 24px 80px}}
h1{{font-size:24px;margin:0 0 2px}} .sub{{color:#64607A;margin:0 0 20px;font-family:monospace;font-size:12px;word-break:break-all}}
h2{{font-size:16px;margin:34px 0 4px}} .cap{{color:#64607A;font-size:13px;margin:0 0 10px}}
.tiles{{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin:14px 0}}
.tile{{background:#fff;border:1px solid #E6E3EE;border-radius:10px;padding:12px 14px}}
.big{{font-size:22px;font-weight:700;font-variant-numeric:tabular-nums}} .lab{{color:#64607A;font-size:11px;text-transform:uppercase;letter-spacing:.04em;margin-top:4px}}
img{{width:100%;border:1px solid #E6E3EE;border-radius:8px;background:#fff}}
table{{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}} th,td{{text-align:left;padding:5px 8px;border-bottom:1px solid #E6E3EE;font-variant-numeric:tabular-nums}}
th{{color:#64607A;font-size:11px;text-transform:uppercase}}
.verdict{{background:{vcol};color:#fff;border-radius:10px;padding:14px 18px;font-weight:600;margin:22px 0}}
.note{{color:#64607A;font-size:12px;margin-top:26px;border-top:1px solid #E6E3EE;padding-top:12px}}
</style>
<div class=wrap>
<h1>Dark-TIC dashboard</h1><p class=sub>{run}</p>
<div class=verdict>{verdict}</div>
<div class=tiles>{tilehtml}</div>
<h2>1 · Where the ion current goes</h2>
<p class=cap>Top: MS1 survey vs MS2 fragment ion current. Bottom: of the MS2 current, how much is explained by identified peptides (teal) vs dark (orange).{" Conservative lower-bound dark (no RT boundaries)." if no_rt else ""}</p>
<img src="{figs['coverage']}">
{mapimg}
<h2>3 · How concentrated is the dark</h2>
<p class=cap>Cumulative share of the dark carried by its most intense features. A few big features rising steeply = worth chasing; a shallow curve needing hundreds of thousands of features = diffuse, low-abundance.</p>
<img src="{figs['conc']}">
<h2>4 · Highest-intensity unexplained features</h2>
<p class=cap>Dark features that are NOT an identified peptide's fragment at another RT, >=0.2% of the dark. Delta-mass column = the feature equals an ID fragment plus this modification (checked against decoy deltas: {nmod}/{len(nov)} real vs {ndec}/{len(nov)} decoy over the top {len(nov)}).</p>
<table><tr><th>m/z</th><th>z</th><th>RT</th><th>% of dark</th><th>delta-mass to an ID fragment</th></tr>{rows}</table>
<p class=note>dark-tic-dashboard (PostdocToolbox). Dark = MS2 intensity not within 10 ppm of an identified peptide's b/y ladder + isotopes + H2O/NH3 losses + precursor{", RT-gated to its elution" if not no_rt else ", subtracted across all RT (conservative)"}. "Explained (any RT)" and the delta-mass check flag spillover and proteoforms. This is a screen, not an identification.</p>
</div>"""


def main():
    ap = argparse.ArgumentParser(description="Visual dark-TIC dashboard for one Skyline+raw pair.")
    ap.add_argument("--report", required=True)
    ap.add_argument("--raw", required=True, help="one .raw file")
    ap.add_argument("--out", default="dark_tic_dashboard.html")
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--top", type=int, default=200)
    ap.add_argument("--all-ids", action="store_true")
    ap.add_argument("--no-rt-gate", action="store_true")
    ap.add_argument("--no-open", action="store_true", help="don't auto-open the HTML")
    for c in ["seq", "pmz", "pch", "rep", "start", "end"]:
        ap.add_argument(f"--col-{c}", default=None)
    a = ap.parse_args()
    cols = dict(seq=a.col_seq, pmz=a.col_pmz, pch=a.col_pch, rep=a.col_rep,
                start=a.col_start, end=a.col_end, all_ids=a.all_ids, no_rt=a.no_rt_gate)
    report = E.load_report(a.report, cols)
    stem = os.path.splitext(os.path.basename(a.raw))[0]
    ids = E.ids_for_raw(stem, report)
    if not ids:
        raise SystemExit(f"[dashboard] no identifications matched raw stem '{stem}'. "
                         f"Report replicates: {list(report)[:3]}...")
    print(f"[dashboard] {stem}: {len(ids):,} IDs; reading raw + building dashboard...")
    r = collect(a.raw, ids, a.ppm, a.no_rt_gate, a.top)
    f = r["fine"]
    big = [x for x in f["novel"] if x["pct_dark"] >= 0.5]
    nmod = sum(1 for x in f["novel"] if x["mod_hits"]); ndec = sum(1 for x in f["novel"] if x["n_decoy_hits"])
    verdict = ("FLAG - high-intensity unexplained proteoform-like signal; worth a targeted open/delta-mass search."
               if (big and nmod > ndec * 1.5) else
               "NEGATIVE - the dark is diffuse and/or spillover; no high-abundance novel species.")
    figs = figures(r, a.no_rt_gate)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(html(r, figs, a.no_rt_gate, stem, verdict))
    print(f"[dashboard] wrote {a.out}")
    if not a.no_open:
        try: webbrowser.open("file://" + os.path.abspath(a.out))
        except Exception: pass


if __name__ == "__main__":
    main()
