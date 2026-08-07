#!/usr/bin/env python
"""Generate LEADERBOARD.md from the cryptic re-mine summaries.

Ranks every gene in the working DB by the cryptic-consistent signature
(canonical peptide DECREASE in the TDP-43 knockdown model), then lists the
disease/plasma-stratum "signal" genes and subpopulation hits, plus a tally.

Usage:  python make_leaderboard.py [--db-dir DIR] [--out FILE]
"""
import argparse, glob, json, os, re
from collections import defaultdict

KD_DATASET = "ineuron_sirna_tdp43kd"
KD_LEVEL = "TDP43_KD"


def is_confound(cond):
    return bool(cond) and bool(re.search(r"(?i)plate|batch|prep|acq", str(cond)))


def real_bio(x):
    """A genuine, detected condition-linked canonical decrease (not a confound/noise)."""
    if not x.get("n_down") or is_confound(x.get("top_condition")):
        return False
    dt = x.get("detected_frac_in_level")
    return not (isinstance(dt, (int, float)) and dt < 0.3)


def load(db_dir):
    by = defaultdict(dict)
    for f in glob.glob(os.path.join(db_dir, "remine_runs", "summaries", "*.json")):
        try:
            d = json.load(open(f, encoding="utf-8"))
            by[d["gene"]][d["name"]] = d
        except Exception:
            pass
    return by


def build(by, generated):
    L = []
    L.append("# Cryptic re-mine leaderboard")
    L.append("")
    L.append(f"_Generated {generated} from {len(by)} genes in the working database._")
    L.append("")
    L.append("**Cryptic-consistent** = a condition-linked DECREASE in a protein's "
             "canonical peptides — the stoichiometry footprint expected when a cryptic "
             "splicing event diverts transcript away from the canonical isoform.")
    L.append("")

    # --- KD-direction leaderboard ---
    kd = []
    for g, ds in by.items():
        x = ds.get(KD_DATASET)
        if x and x.get("status") == "ok" and x.get("top_level") == KD_LEVEL and real_bio(x):
            kd.append((x["n_down"], x["median_effect_log2"], g))
    kd.sort(reverse=True)

    L.append("## 1. KD-direction leaderboard — canonical drop in TDP-43 knockdown")
    L.append("")
    L.append("The primary cryptic-consistent ranking (iNeuron `siRNA=TDP43_KD`).")
    L.append("")
    L.append("| Rank | Gene | Peptides ↓ | median log2 |")
    L.append("|---:|---|---:|---:|")
    for i, (n, e, g) in enumerate(kd, 1):
        L.append(f"| {i} | {g} | {n} | {e:+.2f} |")
    L.append("")
    L.append(f"_{len(kd)} genes show the cryptic-consistent KD-direction decrease._")
    L.append("")

    # --- non-KD signal genes (disease / plasma strata) ---
    rows = []
    for g, ds in by.items():
        hits = []
        for name, x in ds.items():
            if name == KD_DATASET:
                continue
            if x.get("status") == "ok" and real_bio(x):
                hits.append((name, x.get("top_condition"), x.get("top_level"),
                             x["n_down"], x["median_effect_log2"]))
        if hits:
            rows.append((g, hits))
    rows.sort(key=lambda r: -max(h[3] for h in r[1]))

    L.append("## 2. Disease / plasma-stratum signal (non-KD datasets)")
    L.append("")
    L.append("Detected canonical decreases tied to a biological condition outside the KD model "
             "(AD/PD brain strata, plasma-EV cohort, iNeuron controls). Technical axes excluded.")
    L.append("")
    L.append("| Gene | Dataset | Condition = level | Peptides ↓ | log2 |")
    L.append("|---|---|---|---:|---:|")
    for g, hits in rows:
        for name, cond, lvl, n, e in sorted(hits, key=lambda h: -h[3]):
            L.append(f"| {g} | {name} | {cond} = {lvl} | {n} | {e:+.2f} |")
    L.append("")

    # --- subpopulation hits ---
    subs = []
    for g, ds in by.items():
        for name, x in ds.items():
            if x.get("status") == "ok" and x.get("n_subpop_specific"):
                subs.append((g, name, x["n_subpop_specific"]))
    subs.sort(key=lambda s: -s[2])
    L.append("## 3. Subpopulation hits (specificity-controlled)")
    L.append("")
    L.append("| Gene | Dataset | Specific samples |")
    L.append("|---|---|---:|")
    for g, name, n in subs:
        L.append(f"| {g} | {name} | {n} |")
    L.append("")

    # --- tally ---
    tally = defaultdict(int)
    for g, ds in by.items():
        det = [x for x in ds.values() if x.get("status") == "ok"]
        strict = [x for x in det if real_bio(x)]
        sub = [x for x in det if x.get("n_subpop_specific")]
        if not det:
            v = "not measurable"
        elif strict:
            v = "signal"
        elif sub:
            v = "subpopulation"
        elif any(x.get("n_down") for x in det):
            v = "technical only"
        else:
            v = "flat"
        tally[v] += 1
    L.append("## 4. Verdict tally")
    L.append("")
    L.append("| Verdict | Genes |")
    L.append("|---|---:|")
    for k in ["signal", "subpopulation", "technical only", "flat", "not measurable"]:
        L.append(f"| {k} | {tally[k]} |")
    L.append(f"| **total** | **{len(by)}** |")
    L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db-dir", default=r"G:\Manuscripts\LatimerCryptic\cryptic_db")
    ap.add_argument("--out", default=None)
    ap.add_argument("--date", default="unknown date")
    a = ap.parse_args()
    out = a.out or os.path.join(a.db_dir, "LEADERBOARD.md")
    by = load(a.db_dir)
    open(out, "w", encoding="utf-8").write(build(by, a.date))
    print(f"wrote {out}  ({len(by)} genes)")


if __name__ == "__main__":
    main()
