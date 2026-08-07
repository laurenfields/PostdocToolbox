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

# Disease-direction classification: the miner reports the level where a protein is
# LOWEST. A genuine disease-DECREASE target is lowest in a DISEASE/high-pathology
# group; "lowest in a control group" means the protein is actually UP in disease.
CONTROL_LEVELS = {"no dementia", "hcf", "hcn", "pdcn", "hc", "not-als", "control",
                  "normal", "cognitively normal", "healthy", "nonad", "lowad", "non", "low"}
DISEASE_LEVELS = {"dementia", "add", "autosomal dominant add", "sporadic add", "pdd",
                  "als", "pd", "hd", "ftd", "ftld", "ftld-tdp", "lbd", "dlb",
                  "highad", "high"}


def disease_direction(level):
    """DOWN_in_disease (target) / UP_in_disease / ambiguous, from the lowest level."""
    l = str(level).strip().lower()
    if l in DISEASE_LEVELS:
        return "down"      # lowest in disease  -> real decrease target
    if l in CONTROL_LEVELS:
        return "up"        # lowest in control  -> protein rises in disease
    return "ambiguous"     # numeric stage / genotype / other


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

    # --- disease-DECREASE targets, direction-classified ---
    # For each gene take its strongest non-KD, non-confound decrease; classify by
    # whether the lowest level is a disease group (real decrease) or a control group.
    dd = {"down": [], "up": [], "ambiguous": []}
    for g, ds in by.items():
        best = None
        for name, x in ds.items():
            if name == KD_DATASET or x.get("status") != "ok" or not real_bio(x):
                continue
            if best is None or x["n_down"] > best["n_down"]:
                best = {"name": name, "cond": x.get("top_condition"),
                        "lvl": x.get("top_level"), "n_down": x["n_down"],
                        "eff": x["median_effect_log2"]}
        if best:
            dd[disease_direction(best["lvl"])].append((best["n_down"], best["eff"], g, best))

    L.append("## 2. Disease-DECREASE targets (direction-classified)")
    L.append("")
    L.append("Strongest non-KD condition-linked decrease per gene, split by whether the "
             "protein is lowest in a **disease** group (a genuine decrease-in-disease "
             "target) vs lowest in a **control** group (protein actually rises in disease). "
             "Any neurodegenerative stratum counts — AD/PD/FTD/DLB/ALS, brain/CSF/plasma.")
    L.append("")
    L.append("### 2a. DOWN in disease — decrease targets")
    L.append("")
    L.append("| Gene | Peptides ↓ | log2 | Dataset / condition = level |")
    L.append("|---|---:|---:|---|")
    for n, e, g, b in sorted(dd["down"], reverse=True):
        L.append(f"| {g} | {n} | {e:+.2f} | {b['name']} / {b['cond']} = {b['lvl']} |")
    L.append("")
    L.append(f"_{len(dd['down'])} decrease-in-disease targets; "
             f"{len(dd['ambiguous'])} ambiguous (numeric stage/genotype); "
             f"{len(dd['up'])} lowest-in-control (up in disease — see §3)._")
    L.append("")
    if dd["ambiguous"]:
        L.append("### 2b. Ambiguous (stage/genotype — check direction)")
        L.append("")
        L.append("| Gene | Peptides ↓ | log2 | Dataset / condition = level |")
        L.append("|---|---:|---:|---|")
        for n, e, g, b in sorted(dd["ambiguous"], reverse=True):
            L.append(f"| {g} | {n} | {e:+.2f} | {b['name']} / {b['cond']} = {b['lvl']} |")
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

    L.append("## 3. Disease / plasma-stratum signal — full detail (non-KD)")
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
    L.append("## 4. Subpopulation hits (specificity-controlled)")
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
    L.append("## 5. Verdict tally")
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
