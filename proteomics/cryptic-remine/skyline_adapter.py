#!/usr/bin/env python
"""Turn a Skyline transition-level report (CSV/TSV/parquet) into a minimal
PRISM-shaped directory the cryptic re-miner can consume unchanged.

Why: Panorama / collaborators publish Skyline documents, not PRISM output. Export
(or pull) a Skyline report once and this builds:

  <out_dir>/merged_data.parquet    slim transition-level table (gene, peptide,
                                   replicate, area, detq, truncated) -> detection
                                   + gene->peptide fallback
  <out_dir>/peptides.parquet       peptide x replicate matrix (transitions summed,
                                   PRISM-style linear area) -> abundance / limma
  <out_dir>/protein_groups.csv     LeadingGeneName / LeadingUniProtID / AllPeptides
  <out_dir>/sample_metadata.csv    one row per replicate (sample, sample_type)

Then register it in datasets.yaml with:
  peptide_matrix: peptides.parquet
  merged_data: merged_data.parquet
and add name_conditions / a clinical join for the biological contrast.

The report's own column names are used as-is; the miner's COL_ALIASES already
accepts the modern Skyline names (PeptideModifiedSequenceUnimodIds, ProteinGene,
DetectionQValue, Truncated, Area, ReplicateName).
"""
import argparse, os
from pathlib import Path
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.dataset as pads
import pyarrow.parquet as pq
import pandas as pd

# canonical Skyline report columns (with a couple of fallbacks)
C_GENE = ["ProteinGene", "Protein Gene"]
C_ACC  = ["ProteinAccession", "Protein Accession"]
C_PROT = ["Protein", "Protein Name", "ProteinName"]
C_PEP  = ["PeptideModifiedSequenceUnimodIds", "Peptide Modified Sequence"]
C_REP  = ["ReplicateName", "Replicate Name", "Replicate"]
C_AREA = ["Area", "Total Area", "Total Area Fragment"]
C_DETQ = ["DetectionQValue", "Detection QValue"]
C_TRUNC = ["Truncated"]


def _pick(schema_names, options, required=True, what=""):
    have = set(schema_names)
    for o in options:
        if o in have:
            return o
    if required:
        raise SystemExit(f"Skyline report missing a {what} column (looked for {options}); "
                         f"have e.g. {list(schema_names)[:12]}")
    return None


def import_skyline_report(report_path, out_dir, fmt=None):
    report_path = str(report_path)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    if fmt is None:
        fmt = "parquet" if report_path.lower().endswith(".parquet") else "csv"
    if fmt == "csv":
        popt = pads.CsvFileFormat(parse_options=pacsv.ParseOptions(
            delimiter="\t" if report_path.lower().endswith((".tsv", ".txt")) else ","))
        dataset = pads.dataset(report_path, format=popt)
    else:
        dataset = pads.dataset(report_path, format="parquet")

    names = dataset.schema.names
    pep = _pick(names, C_PEP, what="peptide")
    rep = _pick(names, C_REP, what="replicate"); area = _pick(names, C_AREA, what="area")
    acc = _pick(names, C_ACC, required=False); detq = _pick(names, C_DETQ, required=False)
    trunc = _pick(names, C_TRUNC, required=False)
    # gene column, falling back to the protein name when the doc has no gene annotation
    prot = _pick(names, C_PROT, required=False)
    gene = _pick(names, C_GENE, required=False) or prot
    if gene is None:
        raise SystemExit("Skyline report needs a gene or protein column "
                         f"(looked for {C_GENE + C_PROT}); have e.g. {list(names)[:12]}")
    acc = acc or prot

    # de-dupe: with no gene/accession columns, gene==acc==prot all point at the
    # same Protein column; selecting it twice yields a 2-D (duplicate-named) frame.
    keep = list(dict.fromkeys(c for c in [gene, acc, pep, rep, area, detq, trunc] if c))
    tbl = dataset.to_table(columns=keep)
    n_tx = tbl.num_rows

    # 1) slim merged_data.parquet (transition level, names preserved for COL_ALIASES)
    pq.write_table(tbl, out / "merged_data.parquet")

    df = tbl.to_pandas()
    df[area] = pd.to_numeric(df[area], errors="coerce")

    # 2) peptide x replicate matrix: SUM transition area (PRISM sum rollup), linear
    roll = (df.groupby([pep, rep], observed=True)[area].sum().reset_index())
    n_tx_per = df.groupby(pep, observed=True)[area].size().rename("n_transitions")
    mat = roll.pivot(index=pep, columns=rep, values=area)
    mat.insert(0, "n_transitions",
               df.groupby([pep, rep], observed=True).size().groupby(pep).mean().round().astype(int))
    mat = mat.reset_index()               # peptide key column kept with its Skyline name
    mat.columns.name = None
    mat.to_parquet(out / "peptides.parquet", index=False)
    n_pep = mat.shape[0]; n_rep = mat.shape[1] - 2

    # 3) protein_groups.csv (one group per gene; parsimony not attempted)
    pg_cols = list(dict.fromkeys([gene, pep] + ([acc] if acc else [])))
    pg = df[pg_cols].drop_duplicates()
    rows = []
    for g, sub in pg.groupby(gene, observed=True):
        rows.append({
            "LeadingGeneName": g,
            "LeadingUniProtID": (sub[acc].dropna().iloc[0] if acc and sub[acc].notna().any() else ""),
            "AllPeptides": ";".join(sorted(sub[pep].dropna().unique())),
        })
    pd.DataFrame(rows).to_csv(out / "protein_groups.csv", index=False)

    # 4) sample_metadata.csv
    reps = pd.DataFrame({"sample": sorted(df[rep].dropna().unique())})
    reps["sample_type"] = "experimental"
    reps.to_csv(out / "sample_metadata.csv", index=False)

    return {"transitions": n_tx, "peptides": n_pep, "replicates": n_rep,
            "genes": len(rows), "out_dir": str(out)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, help="Skyline report (csv/tsv/parquet)")
    ap.add_argument("--out-dir", required=True, help="destination PRISM-shaped dir")
    ap.add_argument("--format", default=None, choices=[None, "csv", "parquet"])
    a = ap.parse_args()
    r = import_skyline_report(a.report, a.out_dir, fmt=a.format)
    print(f"imported {r['transitions']:,} transitions -> {r['peptides']:,} peptides x "
          f"{r['replicates']} replicates, {r['genes']} genes")
    print(f"wrote: merged_data.parquet, peptides.parquet, protein_groups.csv, "
          f"sample_metadata.csv  in {r['out_dir']}")
    print("Register in datasets.yaml with  peptide_matrix: peptides.parquet  "
          "merged_data: merged_data.parquet")


if __name__ == "__main__":
    main()
