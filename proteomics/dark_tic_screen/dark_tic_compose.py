"""dark-tic-compose - what the dark is MADE OF, aggregated over a subset of runs.

The cohort tool says WHERE the dark is; this says WHAT it is. It runs the fine tier
(per-scan peak decode + ID subtraction) on a small subset of runs and aggregates:

  1. Mass-defect band - a 2D density of the dark peaks (m/z vs fractional mass), the
     signature of composition: peptidic ions form a tight diagonal band, chemical /
     polymer noise sits off it. The single most telling "is the dark real peptide?" view.
  2. Dark concentration - a few big features vs a diffuse haze (per-run, summarized).
  3. Unexplained features that RECUR across the subset runs - reproducible dark species
     (not an identified fragment at another RT), with their delta-mass to a known
     modification. Recurrence across independent runs is the strongest single-tool signal.

This is the expensive tier (~10 min/run of peak decode), so keep the subset small (default
3 patient runs; override with --runs). Reuses the dark_tic_screen engine.

    python dark_tic_compose.py --report merged_data.parquet --raws <dir> --subset 3 --out compose.html
    python dark_tic_compose.py --report merged_data.parquet --raws <dir> --runs stemA,stemB --out compose.html
"""
from __future__ import annotations
import os, sys, glob, argparse, io, base64, importlib.util, time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("dts", os.path.join(_HERE, "dark_tic_screen.py"))
E = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(E)

MZLO, MZHI, NB, ND = 400.0, 905.0, 202, 120   # mass-defect histogram grid


def fine_pass(rawpath, ids, ppm, rt_lo, rt_hi, top):
    from thermo_raw.thermo_raw_reader import assign_window_index
    sys.path.insert(0, os.path.join(_HERE, ".."))
    with E.open_run(rawpath) as rf:
        wt = rf.window_table()
        by = E.build_window_ids(ids, wt)
        win_full = {w: np.sort(np.unique(np.concatenate([e["frag"] for e in es])))
                    for w, es in by.items()}
        valid = {w["win_idx"] for w in wt}
        acc = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        dark_total = 0.0
        H = np.zeros((NB, ND))
        for sn in range(rf.first_scan, rf.last_scan + 1):
            h = rf.header(sn, with_inj=False)
            if h.ms_level != 2 or h.rt_min < rt_lo or h.rt_min > rt_hi:
                continue
            w = assign_window_index(h.center_mz, wt)
            if w not in valid:
                continue
            es = by.get(w)
            active = [e for e in es if e["start"] <= h.rt_min <= e["end"]] if es else []
            mz, it = rf.peaks(sn)
            if mz.size == 0:
                continue
            T = (np.sort(np.unique(np.concatenate([e["frag"] for e in active]))) if active
                 else np.empty(0))
            dark = ~E.matched_mask(mz, T, ppm) if T.size else np.ones(mz.size, bool)
            dmz, dit = mz[dark], it[dark]
            if dmz.size == 0:
                continue
            dark_total += float(dit.sum())
            sel = (dmz >= MZLO) & (dmz < MZHI)
            m, iw = dmz[sel], dit[sel]
            if m.size:
                defect = m - np.floor(m)
                xi = ((m - MZLO) / (MZHI - MZLO) * NB).astype(int).clip(0, NB - 1)
                yi = (defect * ND).astype(int).clip(0, ND - 1)
                np.add.at(H, (xi, yi), iw)
            b = np.round(np.log(dmz) / np.log(1 + ppm * 1e-6)).astype(int)
            rbf = int(h.rt_min * 4); wm = acc[w]
            for bb, ii in zip(b, dit):
                wm[int(bb)][rbf] += float(ii)
    fine = E._fine(acc, win_full, dark_total, ppm, top)
    return H, fine


def _worker(report, raw, replicate, cols, rt_lo, rt_hi, top):
    try:
        c = dict(cols); c["replicate"] = replicate
        ids = (E.report_from_transitions(report, replicate, c).get(replicate, [])
               if E.is_transition_level(report) else E.ids_for_raw(replicate, E.load_report(report, cols)))
        if not ids:
            return replicate, None, {"error": "no IDs"}
        H, fine = fine_pass(raw, ids, 10.0, rt_lo, rt_hi, top)
        return replicate, H, fine
    except Exception as e:
        return replicate, None, {"error": repr(e)}


def is_patient(s):
    return not any(k in s for k in ("Pool", "GPF", "Tot", "-MM-", "CATLAT", "ADPD"))


def _png(fig):
    import matplotlib.pyplot as plt
    b = io.BytesIO(); fig.savefig(b, format="png", dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


def figures(H, fines):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ORANGE, INK, MUT = "#C65A1E", "#191723", "#64607A"
    out = {}
    # 1. mass-defect band
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    with np.errstate(divide="ignore"):
        img = np.log10(np.where(H.T > 0, H.T, np.nan))
    im = ax.imshow(img, origin="lower", aspect="auto", cmap="magma",
                   extent=[MZLO, MZHI, 0, 1])
    xs = np.linspace(MZLO, MZHI, 50)
    ax.plot(xs, (xs * 4.95e-4) % 1.0, color="#37C2B3", lw=1.2, ls="--", label="peptide band (~averagine)")
    ax.set_xlabel("dark peak m/z"); ax.set_ylabel("fractional mass (defect)")
    ax.set_title("Mass-defect band of the dark peaks", fontsize=11)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    cb = fig.colorbar(im, ax=ax); cb.set_label("log10 dark intensity")
    out["defect"] = _png(fig)
    # 2. concentration (per-run anchors)
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    for fn in fines:
        ax.plot([1, 500], [fn["top1_pct"], fn["top500_pct"]], "o-", color=ORANGE, alpha=0.6, ms=4)
    ax.set_xscale("log"); ax.set_xlabel("top-N dark features"); ax.set_ylabel("cumulative % of dark")
    ax.set_ylim(0, 100); ax.set_title("Dark concentration (per run)", fontsize=11)
    med = int(np.median([f["feats_to_50pct"] for f in fines]))
    ax.annotate(f"~{med:,} features = 50% of dark", (2, 55), fontsize=9, color=INK)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    out["conc"] = _png(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description="What the dark is made of, over a subset of runs (fine tier).")
    ap.add_argument("--report", required=True)
    ap.add_argument("--raws", required=True)
    ap.add_argument("--out", default="dark_tic_compose.html")
    ap.add_argument("--runs", default=None, help="comma-separated replicate stems (overrides --subset)")
    ap.add_argument("--subset", type=int, default=3, help="how many patient runs to use if --runs not given")
    ap.add_argument("--rt-range", default="2,23")
    ap.add_argument("--top", type=int, default=200)
    ap.add_argument("--qvalue-max", type=float, default=0.01)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--no-open", action="store_true")
    for c in ["seq", "pmz", "pch", "rep", "start", "end"]:
        ap.add_argument(f"--col-{c}", default=None)
    a = ap.parse_args()
    rt_lo, rt_hi = (float(x) for x in a.rt_range.split(","))
    cols = dict(seq=a.col_seq, pmz=a.col_pmz, pch=a.col_pch, rep=a.col_rep,
                start=a.col_start, end=a.col_end, all_ids=False, no_rt=False,
                qmax=a.qvalue_max, replicate=None)
    raws = sorted(glob.glob(a.raws if any(c in a.raws for c in "*?[") else os.path.join(a.raws, "*.raw"))
                  + ([] if any(c in a.raws for c in "*?[") else glob.glob(os.path.join(a.raws, "*.mzML"))))
    stems = {os.path.splitext(os.path.basename(p))[0]: p for p in raws}
    reps = (E.distinct_replicates(a.report, 100000) if E.is_transition_level(a.report)
            else list(E.load_report(a.report, cols)))

    def match(rep):
        if rep in stems: return stems[rep]
        for st, p in stems.items():
            if rep in st or st in rep: return p
        return None
    if a.runs:
        chosen = [r.strip() for r in a.runs.split(",")]
    else:
        chosen = [r for r in reps if is_patient(r) and match(r)][:a.subset]
    pairs = [(r, match(r)) for r in chosen if match(r)]
    if not pairs:
        raise SystemExit("[compose] no runs matched. Use --runs with exact replicate stems.")
    print(f"[compose] fine-tier composition on {len(pairs)} runs: {[p[0][:40] for p in pairs]}")

    Hsum = np.zeros((NB, ND)); fines = []; feat_runs = defaultdict(list)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_worker, a.report, raw, rep, cols, rt_lo, rt_hi, a.top): rep for rep, raw in pairs}
        for fut in as_completed(futs):
            rep, H, fine = fut.result()
            if H is None:
                print(f"  [skip] {rep}: {fine.get('error')}"); continue
            Hsum += H; fines.append(fine)
            for x in fine["novel"]:
                feat_runs[round(x["mz"], 2)].append((rep, x))
            print(f"  done {rep[:44]} ({time.time()-t0:.0f}s)")
    if not fines:
        raise SystemExit("[compose] no runs succeeded.")
    # features recurring across >=2 runs
    recurring = []
    for mzk, hits in feat_runs.items():
        runs_seen = {h[0] for h in hits}
        if len(runs_seen) >= 2:
            x = hits[0][1]
            recurring.append((mzk, len(runs_seen), np.mean([h[1]["pct_dark"] for h in hits]),
                              x["z"], ",".join(x["mod_hits"]) or "-"))
    recurring.sort(key=lambda r: (-r[1], -r[2]))
    figs = figures(Hsum, fines)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(html(figs, fines, recurring, len(pairs), a.rt_range))
    print(f"[compose] {len(fines)} runs in {(time.time()-t0)/60:.1f} min -> {a.out}")
    if not a.no_open:
        try:
            import webbrowser; webbrowser.open("file://" + os.path.abspath(a.out))
        except Exception:
            pass


def html(figs, fines, recurring, n, rt_range):
    import numpy as np
    md = int(np.median([f["feats_to_50pct"] for f in fines]))
    t1 = np.median([f["top1_pct"] for f in fines])
    rows = "".join(f"<tr><td>{mz:.2f}</td><td>{z}</td><td>{nr}</td><td>{pd:.2f}%</td><td>{mod}</td></tr>"
                   for mz, nr, pd, z, mod in recurring[:20]) or \
        "<tr><td colspan=5 style=color:#64607A>no unexplained feature recurred across runs above threshold</td></tr>"
    return f"""<!doctype html><meta charset=utf8><title>Dark composition</title>
<style>
body{{font:15px/1.55 -apple-system,Segoe UI,sans-serif;color:#191723;background:#FAF9FC;margin:0}}
.wrap{{max-width:860px;margin:0 auto;padding:32px 24px 80px}}
h1{{font-size:24px;margin:0 0 4px}} .sub{{color:#64607A;font-size:13px;margin:0 0 18px}}
h2{{font-size:16px;margin:30px 0 4px}} .cap{{color:#64607A;font-size:13px;margin:0 0 10px}}
img{{width:100%;border:1px solid #E6E3EE;border-radius:8px;background:#fff}}
table{{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}} th,td{{text-align:left;padding:5px 8px;border-bottom:1px solid #E6E3EE}}
th{{color:#64607A;font-size:11px;text-transform:uppercase}} td{{font-variant-numeric:tabular-nums}}
.note{{color:#64607A;font-size:12px;margin-top:26px;border-top:1px solid #E6E3EE;padding-top:12px}}
</style>
<div class=wrap>
<h1>What the dark is made of</h1><p class=sub>fine tier over {n} runs &middot; RT {rt_range} min &middot; median top-feature {t1:.2f}% of dark</p>
<h2>1 &middot; Mass-defect band</h2>
<p class=cap>Density of the dark peaks by m/z (x) and fractional mass / defect (y), summed over the runs. Peptidic ions fall on a tight diagonal band (dashed = approximate peptide/averagine line); chemical or polymer signal sits off it. This is the "is the dark real peptide" test.</p>
<img src="{figs['defect']}">
<h2>2 &middot; Concentration</h2>
<p class=cap>Cumulative share of the dark carried by its most intense features, per run. A shallow curve needing ~{md:,} features to reach 50% means the dark is a diffuse low-abundance haze, not a few big species.</p>
<img src="{figs['conc']}">
<h2>3 &middot; Unexplained features recurring across runs</h2>
<p class=cap>Dark features that are not an identified fragment at another RT and appear in &ge;2 of the {n} runs - reproducible dark species. Delta-mass = the feature equals an identified fragment plus this modification.</p>
<table><tr><th>m/z</th><th>z</th><th>runs</th><th>mean % of dark</th><th>delta-mass</th></tr>{rows}</table>
<p class=note>dark-tic-compose (PostdocToolbox). Fine tier: dark = MS2 intensity not within 10 ppm of an identified peptide's b/y ladder + isotopes + losses + precursor, RT-gated. Mass defect is intensity-weighted over all dark peaks 400-905 m/z. A screen, not an identification.</p>
</div>"""


if __name__ == "__main__":
    main()
