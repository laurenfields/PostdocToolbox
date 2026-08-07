#!/usr/bin/env python
"""
cryptic_remine.py -- log, database, and canonical re-mining for cryptic-peptide analyses.

Premise
-------
Producing a cryptic peptide changes the stoichiometry of the parent protein: a
cryptic/skiptic splicing event (a TDP-43 loss-of-function phenomenon) alters the
transcript, so the *canonical* protein -- or the canonical peptides spanning/adjacent
to the affected region -- DECREASES in the affected condition. We always observe the
canonical tryptic peptides in a PRISM search, so that decrease is already sitting in
completed results. When a paper reports that protein X makes a cryptic peptide, we add
that context to a database and re-mine existing PRISM runs for a condition-linked drop
in X's canonical peptides that we previously overlooked.

The reported cryptic peptide does NOT need to be in our search FASTA -- the database
entry is context that tells us which canonical protein to zoom into. (If the cryptic
peptide happens to be in the results, `mine` also runs a two-sided sanity check.)

This is a SCREEN (triage), not a claim. Promotion to a reported finding uses the
moderated limma/DEqMS tool (cryptic-intensity-prior-limma), with a batch covariate.

Subcommands
-----------
  seed-db   Build the cryptic-peptide database from the search FASTA custom headers.
  add       Register a protein/cryptic peptide (from a paper) as cryptic-relevant.
  mine      Re-mine registered PRISM datasets for a gene's canonical-peptide changes.
  log       Append a run summary to the human-readable analysis log.

Self-contained: reads PRISM parquet directly (pandas/pyarrow); no cross-repo imports.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

# --------------------------------------------------------------------------- #
# Configuration / paths
# --------------------------------------------------------------------------- #
DEFAULT_DB_DIR = Path(
    os.environ.get("CRYPTIC_DB_DIR", r"G:\Manuscripts\LatimerCryptic\cryptic_db")
)
DB_CSV_NAME = "cryptic_peptides.csv"
REGISTRY_NAME = "datasets.yaml"
RUNS_SUBDIR = "remine_runs"
LOG_NAME = "cryptic_analysis_log.md"  # written to db_dir.parent

DB_COLUMNS = [
    "entry_id", "gene", "canonical_uniprot", "cryptic_peptide_sequence",
    "cryptic_type", "source", "doi_pmid", "date_added", "in_fasta",
    "mined_datasets", "notes",
]

# Screen thresholds (documented in README; overridable via CLI where noted)
MIN_LEVEL_N = 3          # skip a condition level with fewer than this many samples
EFFECT_DOWN_LOG2 = -1.0  # a peptide is "down" in a level if median drops >= 1 log2 (2x)
STABLE_ABS_LOG2 = 0.5    # a peptide is "stable" if |effect| < this
ADJ_P_CUTOFF = 0.10      # BH-adjusted KW p to call a peptide's change notable
DETQ_CUTOFF = 0.01       # DetectionQValue threshold for a "truly detected" transition

# Condition-level labels that are technical (QC/reference/pool injections), never a
# biological condition -- dropped from every scan as a backstop.
TECHNICAL_LEVEL_RE = re.compile(r"(?:batchqc|\bqc\b|pool|totalref|\bref\b)", re.IGNORECASE)

# A condition with more distinct levels than this is treated as continuous / an ID
# (e.g. Donor Age) and skipped -- categorical scans over continuous fields manufacture
# spurious tiny-group "hits".
MAX_CONDITION_LEVELS = 12

# "Dive deeper" (liberal) tier -- surface anything worth a closer look, not just the
# strict concordant decrease. Deliberately generous; these are leads, not findings.
SUGGESTIVE_EFFECT_LOG2 = -0.6   # whole-group trend: median down >= this ...
SUGGESTIVE_P = 0.05             # ... at nominal (uncorrected) KW p
SUBSET_Z = -2.5                 # robust z below which one sample is a low-outlier
SUBSET_MIN_OUTLIERS = 2         # >= this many low-outliers in a level -> a subset flag

# Metadata-AGNOSTIC discovery: find samples that drop MANY of a protein's canonical
# peptides at once (a coherent subpopulation with reduced canonical protein), using no
# labels at all. Metadata, when present, only annotates what was found -- never drives it.
UNSUP_MIN_PEPTIDES = 3          # a sample must be low across >= this many canonical peptides
UNSUP_MIN_FRAC = 0.34           # ... and across >= this fraction of its testable peptides
# Specificity control: a subpopulation sample is only "protein-specific" if its low-rate
# for THIS protein exceeds its low-rate across a random proteome background by this margin
# (else it's just a globally low-input sample -- a loading artifact, not biology).
SPECIFICITY_MARGIN = 0.25
BACKGROUND_PEPTIDES = 3000      # random background size for the global low-rate


def today() -> str:
    return _dt.date.today().isoformat()


# --------------------------------------------------------------------------- #
# Schema tolerance -- PRISM output column names differ across pipeline versions
# (modern: `PeptideModifiedSequenceUnimodIds`, `ProteinGene`, `DetectionQValue`;
#  older iNeuron-era: `Peptide Modified Sequence`, `Protein Gene`, `Detection QValue`).
# --------------------------------------------------------------------------- #
COL_ALIASES = {
    "peptide": ["PeptideModifiedSequenceUnimodIds", "Peptide Modified Sequence"],
    "gene":    ["ProteinGene", "Protein Gene"],
    "protein": ["Protein", "Protein Name"],
    "detq":    ["DetectionQValue", "Detection QValue"],
    "trunc":   ["Truncated"],
    "area":    ["Area"],
    "sampleid": ["Sample ID", "Replicate Name", "ReplicateName"],
}


def resolve_cols(schema_names, needed) -> dict:
    """Map each canonical name in `needed` to the actual column present, or None."""
    have = set(schema_names)
    out = {}
    for canon in needed:
        out[canon] = next((a for a in COL_ALIASES[canon] if a in have), None)
    return out


def peptide_key(schema_names) -> str:
    k = next((a for a in COL_ALIASES["peptide"] if a in set(schema_names)), None)
    if k is None:
        raise KeyError("No peptide-sequence column found in matrix "
                       f"(looked for {COL_ALIASES['peptide']}).")
    return k


def matrix_sample_columns(schema_names) -> list[str]:
    """Sample columns = everything except the peptide key and per-peptide annotations."""
    non_sample = set(COL_ALIASES["peptide"]) | {"n_transitions", "mean_rt"}
    return [c for c in schema_names if c not in non_sample]


# --------------------------------------------------------------------------- #
# FASTA cryptic-header parsing  ->  database seed
# --------------------------------------------------------------------------- #
CRYPTIC_CATEGORIES = ("CRYPTIC", "CRYPTIC_UNIQUE", "CRYPTIC_EXON", "ALT_EXON")


def _kv(tokens):
    """Split 'key:value' tokens into a dict; return (kv_dict, positional_list)."""
    kv, pos = {}, []
    for t in tokens:
        if ":" in t and not t.lower().startswith("chr"):
            k, v = t.split(":", 1)
            kv[k.strip().lower()] = v.strip()
        else:
            pos.append(t)
    return kv, pos


def parse_cryptic_header(header: str) -> dict | None:
    """
    Parse one cryptic FASTA header (without the leading '>') into a DB row dict.
    Tolerant of the four observed schemas:
      CRYPTIC|GENE|CR_id|type|priority:...|source:...
      CRYPTIC|GENE|GENE_greek|TDP43_knockdown|source:...
      CRYPTIC_UNIQUE|GENE|sp|<coords>|unique:a-b|tryptic:c-d|source:...
      CRYPTIC_EXON|GENE|CE_id|chrX:.. (hg38)|tryptic_pep_N|desc|ref:...
      ALT_EXON|GENE|isoform_unknown|alt_region:a-b|alt_len:..|pep_len:..|C-term|source:UniProt
    """
    toks = header.split("|")
    if not toks or toks[0] not in CRYPTIC_CATEGORIES:
        return None
    category = toks[0]
    gene = toks[1] if len(toks) > 1 else ""
    rest = [t for t in toks[2:] if t != "sp"]  # drop the literal 'sp' filler token
    kv, pos = _kv(rest)

    source = kv.get("source") or kv.get("ref") or ""
    priority = kv.get("priority", "")
    # genomic coordinates: prefer an explicit chr..., else transcript/coords positional,
    # else alt_region.
    coords = ""
    for t in rest:
        if t.lower().startswith("chr") or "ENST" in t:
            coords = t
            break
    if not coords and "alt_region" in kv:
        coords = "alt_region:" + kv["alt_region"]

    # entry id + type per category
    if category == "CRYPTIC":
        entry_core = pos[0] if pos else gene
        ctype = pos[1] if len(pos) > 1 else "cryptic"
    elif category == "CRYPTIC_UNIQUE":
        entry_core = coords or (pos[0] if pos else gene)
        rng = kv.get("unique", "")
        entry_core = f"{entry_core}:{rng}" if rng else entry_core
        ctype = "unique"
    elif category == "CRYPTIC_EXON":
        ce_id = pos[0] if pos else "CE"
        pep = next((p for p in pos if p.startswith("tryptic_pep")), "")
        entry_core = f"{ce_id}_{pep}" if pep else ce_id
        ctype = "cryptic_exon"
    else:  # ALT_EXON
        entry_core = "alt_" + kv.get("alt_region", (pos[0] if pos else ""))
        ctype = "alt_exon"

    entry_id = f"{gene}:{entry_core}"
    return {
        "entry_id": entry_id,
        "gene": gene,
        "canonical_uniprot": "",
        "cryptic_peptide_sequence": "",  # filled from the sequence line
        "cryptic_type": ctype,
        "source": source,
        "doi_pmid": "",
        "date_added": today(),
        "in_fasta": True,
        "mined_datasets": "",
        "notes": f"category={category}" + (f"; coords={coords}" if coords else "")
                 + (f"; priority={priority}" if priority else ""),
    }


def parse_fasta_cryptics(fasta_path: Path) -> pd.DataFrame:
    """Scan a FASTA and return ALL cryptic-family entries (incl. CRYPTIC_UNIQUE)."""
    rows, cur, seq = [], None, []
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    cur["cryptic_peptide_sequence"] = "".join(seq)
                    rows.append(cur)
                seq = []
                cur = parse_cryptic_header(line[1:].strip())
            elif cur is not None:
                seq.append(line.strip())
        if cur is not None:
            cur["cryptic_peptide_sequence"] = "".join(seq)
            rows.append(cur)
    df = pd.DataFrame(rows, columns=DB_COLUMNS)
    return df.drop_duplicates(subset="entry_id", keep="first").reset_index(drop=True)


def seed_db_from_fasta(fasta_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return (curated_db, unique_reference).

    The main working DB is the CURATED cryptic catalog: CRYPTIC (Novartis + Seddighi
    isoform) + CRYPTIC_EXON + ALT_EXON -- entries with real gene symbols and literature
    provenance. The `CRYPTIC_UNIQUE` family (Seddighi DataS2) is an exhaustive
    machine-generated search-space catalog keyed by transcript/accession (tens of
    thousands of rows); it is returned separately as a reference sidecar, not mixed
    into the curated DB.
    """
    allc = parse_fasta_cryptics(fasta_path)
    is_uniq = allc["cryptic_type"] == "unique"
    return allc[~is_uniq].reset_index(drop=True), allc[is_uniq].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Database load / save / add
# --------------------------------------------------------------------------- #
def db_path(db_dir: Path) -> Path:
    return db_dir / DB_CSV_NAME


def load_db(db_dir: Path) -> pd.DataFrame:
    p = db_path(db_dir)
    if p.exists():
        return pd.read_csv(p, dtype=str).fillna("")
    return pd.DataFrame(columns=DB_COLUMNS)


def save_db(df: pd.DataFrame, db_dir: Path) -> None:
    db_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(db_path(db_dir), index=False)


def add_entry(db_dir: Path, gene: str, sequence: str = "", source: str = "",
              doi: str = "", ctype: str = "reported", uniprot: str = "",
              notes: str = "") -> pd.DataFrame:
    """Insert/refresh a cryptic-context row. Gene-only entries are allowed."""
    df = load_db(db_dir)
    seq = (sequence or "").strip().upper()
    entry_id = f"{gene}:{seq[:12]}" if seq else f"{gene}:reported_{today()}"
    fasta_hit = False
    if not df.empty:
        fasta_hit = bool(
            ((df["gene"] == gene) & (df["in_fasta"].astype(str) == "True")).any()
        )
    row = {
        "entry_id": entry_id, "gene": gene, "canonical_uniprot": uniprot,
        "cryptic_peptide_sequence": seq, "cryptic_type": ctype, "source": source,
        "doi_pmid": doi, "date_added": today(),
        "in_fasta": fasta_hit,  # True if this gene already has a FASTA-derived entry
        "mined_datasets": "",
        "notes": notes or "manually added from literature",
    }
    df = df[df["entry_id"] != entry_id]  # replace if same id
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    save_db(df, db_dir)
    return df


# --------------------------------------------------------------------------- #
# Dataset registry
# --------------------------------------------------------------------------- #
def load_registry(db_dir: Path) -> dict:
    p = db_dir / REGISTRY_NAME
    if not p.exists():
        raise FileNotFoundError(f"No dataset registry at {p}. Create it (see README).")
    if yaml is None:
        raise ImportError("pyyaml is required to read the dataset registry.")
    with open(p, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def resolve_datasets(reg: dict, name: str | None) -> list[dict]:
    datasets = reg.get("datasets", [])
    if name:
        datasets = [d for d in datasets if d.get("name") == name]
        if not datasets:
            raise ValueError(f"Dataset '{name}' not found in registry.")
    return datasets


# --------------------------------------------------------------------------- #
# Peptide -> gene resolution & abundance loading
# --------------------------------------------------------------------------- #
def canonical_peptides_from_merged(prism_dir: Path, ds: dict, gene: str) -> list[str]:
    """
    Fallback peptide->gene resolution straight from transition-level merged_data
    (`ProteinGene`). Catches genes that never formed a protein group (e.g. single-peptide
    or low-confidence proteins) -- important for low-abundance cryptic targets in plasma.
    """
    import pyarrow.dataset as pads
    import pyarrow.compute as pc
    md = ds.get("merged_data", "merged_data.parquet")
    p = prism_dir / md
    if not p.exists():
        return []
    dataset = pads.dataset(str(p), format="parquet")
    r = resolve_cols(dataset.schema.names, ["gene", "peptide"])
    gcol, key = r["gene"], r["peptide"]
    if gcol is None or key is None:
        return []
    # ProteinGene can be a semicolon list of genes; match gene as a token.
    filt = pc.match_substring(pc.field(gcol), gene)
    tbl = dataset.to_table(columns=[gcol, key], filter=filt)
    df = tbl.to_pandas()
    if df.empty:
        return []
    tok = df[gcol].astype(str).str.split(r"[;,/]").apply(
        lambda xs: gene in [x.strip() for x in xs])
    peps = df.loc[tok, key].dropna().unique().tolist()
    return sorted(set(peps))


def canonical_peptides_for_gene(prism_dir: Path, gene: str,
                                ds: dict | None = None) -> tuple[list[str], str, str]:
    """
    Return (peptide_modseq_list, leading_uniprot, source_tag).
    Prefer protein_groups.csv (fast, parsimony); fall back to merged_data ProteinGene
    when the gene formed no group or the file is absent (older PRISM output).
    """
    peps: list[str] = []
    uniprot = ""
    pg_path = prism_dir / "protein_groups.csv"
    if pg_path.exists():
        pg = pd.read_csv(pg_path)
        # LeadingGeneName may be a multi-gene group ("STMN2 / STMN3", "A;B") -> token match
        gtok = pg["LeadingGeneName"].astype(str).apply(
            lambda c: gene in [t.strip() for t in re.split(r"[;,/]", c)])
        hits = pg[gtok]
        for _, r in hits.iterrows():
            uniprot = uniprot or str(r.get("LeadingUniProtID", ""))
            allp = str(r.get("AllPeptides", ""))
            if allp and allp.lower() != "nan":
                peps.extend(allp.split(";"))
    peps = sorted(set(p for p in peps if p))
    if peps:
        return peps, uniprot, "protein_groups"
    if ds is not None:  # fallback
        peps = canonical_peptides_from_merged(prism_dir, ds, gene)
        if peps:
            return peps, "", "merged_data"
    return [], "", "none"


def load_abundance(prism_dir: Path, peptides: list[str], sample_cols: list[str],
                   matrix: str = "corrected_peptides.parquet") -> pd.DataFrame:
    """Return log2 abundance, index=peptide, columns=sample_cols (linear->log2)."""
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(prism_dir / matrix)
    key = peptide_key(pf.schema_arrow.names)
    cols = [key] + [c for c in sample_cols if c in pf.schema_arrow.names]
    df = pf.read(columns=cols).to_pandas()
    df = df[df[key].isin(peptides)].set_index(key)
    lin = df[cols[1:]].apply(pd.to_numeric, errors="coerce")
    log2 = np.log2(lin.where(lin > 0))
    return log2


# --------------------------------------------------------------------------- #
# Clinical join
# --------------------------------------------------------------------------- #
def _matrix_colmap(prism_dir: Path, ds: dict) -> tuple[dict, list]:
    """Map both full matrix column names and their base (pre-`__@__`) to the real column."""
    import pyarrow.parquet as pq
    matrix = ds.get("peptide_matrix", "corrected_peptides.parquet")
    sids = matrix_sample_columns(pq.ParquetFile(prism_dir / matrix).schema_arrow.names)
    m = {}
    for c in sids:
        m[c] = c
        m.setdefault(c.split("__@__")[0], c)
    return m, sids


def _sample_frame(prism_dir: Path, ds: dict) -> pd.DataFrame:
    """
    Base sample table indexed by the ACTUAL matrix column name (`sample_id`), with a
    `sample` (base name) and, when available, `batch`. Sample ids from
    sample_metadata.csv are reconciled against the real matrix columns (handles the
    `__@__` suffix and metadata files that lack a `sample_id` column); if there's no
    sample_metadata.csv, the sample list is taken straight from the matrix columns.
    """
    colmap, sids = _matrix_colmap(prism_dir, ds)
    smp = prism_dir / ds.get("sample_metadata", "sample_metadata.csv")
    if smp.exists():
        sm = pd.read_csv(smp)
        if "sample_type" in sm.columns:
            sm = sm[sm["sample_type"] == "experimental"].copy()
        base = sm["sample"] if "sample" in sm.columns else sm.get("sample_id")
        sid = sm["sample_id"] if "sample_id" in sm.columns else base
        resolved = [colmap.get(str(s)) or colmap.get(str(b))
                    for s, b in zip(sid, base)]
        out = pd.DataFrame({"__mc": resolved, "sample": list(base)})
        if "batch" in sm.columns:
            out["batch"] = list(sm["batch"])
        out = out[out["__mc"].notna()].set_index("__mc")
        out.index.name = "sample_id"
        return out
    base = [s.split("__@__")[0] for s in sids]
    return pd.DataFrame({"sample": base}, index=pd.Index(sids, name="sample_id"))


def _apply_name_conditions(sample: pd.Series, rules_by_cond: dict) -> pd.DataFrame:
    """Derive condition columns from the base sample name via ordered regex rules."""
    out = {}
    for cond_name, rules in rules_by_cond.items():
        def label(s, rules=rules):
            for rule in rules:
                if re.search(rule["pattern"], str(s), re.IGNORECASE):
                    return rule["label"]
            return np.nan
        out[cond_name] = sample.map(label)
    return pd.DataFrame(out, index=sample.index)


def build_sample_conditions(prism_dir: Path, ds: dict) -> pd.DataFrame:
    """
    Return a frame indexed by the matrix sample-column name, one column per condition
    to scan. Sources, all optional and combinable:
      - `batch` from sample_metadata (scanned as a confound-check unless scan_batch: false)
      - clinical columns joined from `clinical_metadata` on `clinical_join_key`
      - `name_conditions`: columns derived from the sample name by regex
    Technical samples are dropped via `exclude_sample_regex` and/or `patient_filter_*`.
    """
    frame = _sample_frame(prism_dir, ds)

    # global drop of technical injections by sample-name regex (REF/POOL/HeLa/QC...)
    excl_re = ds.get("exclude_sample_regex")
    if excl_re:
        frame = frame[~frame["sample"].astype(str).str.contains(excl_re, case=False,
                                                                 na=False, regex=True)]

    cond_cols: list[str] = []

    # batch as a confound-check condition
    if ds.get("scan_batch", True) and "batch" in frame.columns \
            and frame["batch"].nunique(dropna=True) > 1:
        cond_cols.append("batch")

    # clinical join (preserve the sample_id index across the merge)
    conds = ds.get("condition_columns", [])
    if ds.get("clinical_metadata") and conds:
        clinical = pd.read_csv(ds["clinical_metadata"], dtype=str)
        join_key = ds.get("clinical_join_key", "MS_FileName")
        filt_col = ds.get("patient_filter_column")
        keep = [c for c in dict.fromkeys(
            [join_key] + conds + ([filt_col] if filt_col else []))
            if c in clinical.columns]
        clinical = clinical[keep].drop_duplicates(subset=join_key)
        frame = (frame.reset_index()
                      .merge(clinical, left_on="sample", right_on=join_key, how="left")
                      .set_index("sample_id"))
        cond_cols += [c for c in conds if c in frame.columns]
        excl = ds.get("patient_filter_exclude", [])
        if filt_col and filt_col in frame.columns and excl:
            tech = frame[filt_col].astype(str).str.contains("|".join(excl),
                                                            case=False, na=False)
            frame = frame[~tech]

    # name-derived conditions
    name_conditions = ds.get("name_conditions")
    if name_conditions:
        derived = _apply_name_conditions(frame["sample"], name_conditions)
        for c in derived.columns:
            frame[c] = derived[c].values
            cond_cols.append(c)

    cond_cols = list(dict.fromkeys(c for c in cond_cols if c in frame.columns))
    return frame[cond_cols]


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment; NaN-safe."""
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    order = np.argsort(q)
    ranked = q[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def _clean_cond(cond: pd.Series) -> pd.Series:
    """Drop NaN/blank and technical (QC/ref/pool) labels."""
    cond = cond.dropna()
    cond = cond[cond.astype(str).str.strip().ne("")]
    return cond[~cond.astype(str).str.contains(TECHNICAL_LEVEL_RE)]


def _scannable_levels(cond: pd.Series) -> list | None:
    """Return usable levels, or None if the column is continuous/ID-like or too sparse."""
    if cond.nunique() > MAX_CONDITION_LEVELS:      # continuous field (e.g. Donor Age)
        return None
    levels = [lv for lv, n in cond.value_counts().items() if n >= MIN_LEVEL_N]
    return levels if len(levels) >= 2 else None


def peptide_condition_scan(log2: pd.DataFrame, cond: pd.Series) -> pd.DataFrame:
    """
    For one condition (a Series of level labels indexed by sample col), test each
    peptide across levels. Returns per-peptide: affected_level (lowest median),
    effect_log2 (affected level median - pooled other levels), kw_p, n_levels,
    plus `notable` (strict) and `suggestive` (liberal) flags.
    """
    from scipy.stats import kruskal
    cond = _clean_cond(cond)
    levels = _scannable_levels(cond)
    recs = []
    if levels is None:
        return pd.DataFrame(recs)
    samples_by_level = {lv: cond.index[cond == lv] for lv in levels}
    for pep, row in log2.iterrows():
        groups, medians = [], {}
        for lv in levels:
            vals = row.reindex(samples_by_level[lv]).dropna().values
            if len(vals) >= MIN_LEVEL_N:
                groups.append(vals)
                medians[lv] = float(np.median(vals))
        if len(groups) < 2:
            continue
        try:
            _, p = kruskal(*groups)
        except ValueError:
            p = np.nan
        aff = min(medians, key=medians.get)
        aff_n = int(np.isfinite(row.reindex(samples_by_level[aff]).values).sum())
        others = np.concatenate(
            [row.reindex(samples_by_level[lv]).dropna().values
             for lv in medians if lv != aff]
        )
        effect = medians[aff] - float(np.median(others)) if others.size else np.nan
        recs.append({"peptide": pep, "affected_level": aff, "affected_level_n": aff_n,
                     "effect_log2": effect,
                     "direction": "down" if effect < 0 else "up",
                     "kw_p": p, "n_levels": len(groups)})
    out = pd.DataFrame(recs)
    if not out.empty:
        out["kw_p_bh"] = bh_adjust(out["kw_p"].values)
        out["notable"] = (out["effect_log2"] <= EFFECT_DOWN_LOG2) & \
                         (out["kw_p_bh"] <= ADJ_P_CUTOFF)
        # liberal tier: a whole-group downward trend at a looser, uncorrected bar
        out["suggestive"] = (out["effect_log2"] <= SUGGESTIVE_EFFECT_LOG2) & \
                            (out["kw_p"] <= SUGGESTIVE_P) & (~out["notable"])
    return out


def subset_outlier_scan(log2: pd.DataFrame, cond: pd.Series) -> pd.DataFrame:
    """
    Liberal "some patients" detector. Independent of group-median shifts: for each
    peptide it finds INDIVIDUAL low-outlier samples (robust z <= SUBSET_Z vs the
    peptide's own median/MAD) and reports how they distribute across condition levels.
    A subset of one group (e.g. some ALS but not PD) dropping is exactly the signal the
    strict median test averages away. Returns one row per (peptide, level) with >=1
    low-outlier.
    """
    cond = _clean_cond(cond)
    levels = _scannable_levels(cond)
    if levels is None:
        return pd.DataFrame()
    keep = cond[cond.isin(levels)]
    recs = []
    for pep, row in log2.iterrows():
        vals = row.reindex(keep.index)
        v = vals.dropna()
        if len(v) < 2 * MIN_LEVEL_N:
            continue
        med = float(v.median())
        mad = float((v - med).abs().median()) * 1.4826
        if not np.isfinite(mad) or mad <= 0:
            continue
        low = ((vals - med) / mad) <= SUBSET_Z
        for lv in levels:
            lv_idx = keep.index[keep == lv]
            lo = low.reindex(lv_idx).fillna(False)
            n_out = int(lo.sum())
            if n_out:
                recs.append({
                    "peptide": pep, "level": lv, "n_outliers": n_out,
                    "n_level": int(vals.reindex(lv_idx).notna().sum()),
                    "outlier_samples": ";".join(map(str, lv_idx[lo.values])),
                })
    df = pd.DataFrame(recs)
    if not df.empty:
        df["outlier_frac"] = df["n_outliers"] / df["n_level"].clip(lower=1)
    return df


def unsupervised_discovery(log2: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fully metadata-agnostic. For each canonical peptide, compute a robust z per sample
    (vs the peptide's own median/MAD). A single low peptide is noise; a sample that is
    low across MANY of the protein's peptides at once is a coherent signal -- a candidate
    subpopulation with reduced canonical protein, discovered without any labels.

    Returns (samples, peptides):
      samples  -- one row per sample: n_low (peptides it's a low-outlier for), n_testable,
                  concordance = n_low/n_testable, mean_low_z; sorted, strongest first.
      peptides -- per peptide: how many samples are low (peptide-level discordance).
    Downstream code attaches whatever metadata exists as *description* only.
    """
    if log2.shape[0] == 0:
        return pd.DataFrame(), pd.DataFrame()
    z = pd.DataFrame(np.nan, index=log2.index, columns=log2.columns)
    for pep, row in log2.iterrows():
        v = row.dropna()
        if len(v) < 2 * MIN_LEVEL_N:          # too few samples to define an outlier
            continue
        med = float(v.median())
        mad = float((v - med).abs().median()) * 1.4826
        if not np.isfinite(mad) or mad <= 0:
            continue
        z.loc[pep] = (row - med) / mad
    z_defined = z.notna()
    low = (z <= SUBSET_Z) & z_defined
    n_low = low.sum(axis=0)
    n_testable = z_defined.sum(axis=0)
    conc = n_low / n_testable.clip(lower=1)
    samples = pd.DataFrame({
        "n_low": n_low.astype(int), "n_testable": n_testable.astype(int),
        "n_detected": log2.notna().sum(axis=0).astype(int),
        "concordance": conc.round(3),
        "mean_low_z": z.where(low).mean(axis=0).round(2),
    })
    samples = samples[(samples["n_low"] >= UNSUP_MIN_PEPTIDES)
                      & (samples["concordance"] >= UNSUP_MIN_FRAC)]
    samples = samples.sort_values(["n_low", "concordance"], ascending=False)
    peptides = pd.DataFrame({
        "n_low_samples": low.sum(axis=1).astype(int),
        "n_testable_samples": z_defined.sum(axis=1).astype(int),
    })
    peptides = peptides[peptides["n_low_samples"] > 0].sort_values(
        "n_low_samples", ascending=False)
    return samples, peptides


def global_low_rate(prism_dir: Path, ds: dict, sample_cols: list[str]) -> pd.Series:
    """
    Per-sample fraction of low-outlier calls across a proteome-wide background (not this
    protein). A sample with a high global low-rate is globally low-input; comparing a
    protein's low-rate to this separates a real protein-specific subpopulation from a
    loading artifact. Uses the same robust-z / SUBSET_Z definition, vectorized.
    """
    import pyarrow.parquet as pq
    matrix = ds.get("peptide_matrix", "corrected_peptides.parquet")
    pf = pq.ParquetFile(prism_dir / matrix)
    key = peptide_key(pf.schema_arrow.names)
    cols = [key] + [c for c in sample_cols if c in pf.schema_arrow.names]
    df = pf.read(columns=cols).to_pandas().set_index(key)
    if len(df) > BACKGROUND_PEPTIDES:                    # deterministic systematic subsample
        df = df.iloc[:: max(1, len(df) // BACKGROUND_PEPTIDES)]
    lin = df.apply(pd.to_numeric, errors="coerce")
    log2 = np.log2(lin.where(lin > 0))
    med = log2.median(axis=1)
    mad = (log2.sub(med, axis=0)).abs().median(axis=1) * 1.4826
    z = log2.sub(med, axis=0).div(mad.where(mad > 0), axis=0)
    zdef = z.notna()
    low = (z <= SUBSET_Z) & zdef
    return (low.sum(axis=0) / zdef.sum(axis=0).clip(lower=1))


def score_condition(scan: pd.DataFrame) -> dict:
    """Aggregate a per-peptide scan into a condition-level score for ranking."""
    if scan.empty:
        return {"n_down": 0, "n_stable": 0, "modal_level": "", "median_effect": np.nan,
                "discordant": False, "down_peptides": []}
    down = scan[(scan["effect_log2"] <= EFFECT_DOWN_LOG2)
                & (scan["kw_p_bh"] <= ADJ_P_CUTOFF)]
    modal = down["affected_level"].mode().iloc[0] if not down.empty else ""
    down_in_modal = down[down["affected_level"] == modal] if modal else down
    stable = scan[scan["effect_log2"].abs() < STABLE_ABS_LOG2]
    return {
        "n_down": int(len(down_in_modal)),
        "n_stable": int(len(stable)),
        "modal_level": modal,
        "median_effect": float(down_in_modal["effect_log2"].median())
        if len(down_in_modal) else np.nan,
        "discordant": bool(len(down_in_modal) > 0 and len(stable) > 0),
        "down_peptides": down_in_modal["peptide"].tolist(),
    }


# --------------------------------------------------------------------------- #
# Detection recovery (optional, reads transition-level merged_data)
# --------------------------------------------------------------------------- #
def detection_by_level(prism_dir: Path, ds: dict, peptides: list[str],
                       cond: pd.Series, level: str) -> float | None:
    """Fraction of (peptide, sample) pairs in `level` with >=1 truly-detected transition."""
    import pyarrow.dataset as pads
    import pyarrow.compute as pc
    md = ds.get("merged_data", "merged_data.parquet")
    dataset = pads.dataset(str(prism_dir / md), format="parquet")
    r = resolve_cols(dataset.schema.names, ["peptide", "sampleid", "detq", "trunc", "area"])
    key, sid = r["peptide"], r["sampleid"]
    if key is None or sid is None or r["area"] is None:
        return None
    cols = [c for c in [key, sid, r["detq"], r["trunc"], r["area"]] if c]
    tbl = dataset.to_table(columns=cols, filter=pc.field(key).isin(peptides))
    df = tbl.to_pandas()
    if df.empty:
        return None
    area = pd.to_numeric(df[r["area"]], errors="coerce")
    det = area > 0
    if r["trunc"]:
        det &= ~df[r["trunc"]].astype(str).str.lower().isin(["true", "1", "1.0"])
    if r["detq"]:
        det &= pd.to_numeric(df[r["detq"]], errors="coerce") <= DETQ_CUTOFF
    df["_det"] = det
    # merged_data sample id matches cond.index (the matrix column name)
    level_samples = set(cond.index[cond == level])
    sub = df[df[sid].astype(str).isin(level_samples)]
    if sub.empty:
        return None
    return float(sub.groupby([key, sid])["_det"].any().mean())


# --------------------------------------------------------------------------- #
# Mining driver
# --------------------------------------------------------------------------- #
def mine_gene(db_dir: Path, gene: str, dataset_name: str | None,
              detection: bool = True) -> dict:
    reg = load_registry(db_dir)
    datasets = resolve_datasets(reg, dataset_name)
    runs_dir = db_dir / RUNS_SUBDIR
    runs_dir.mkdir(parents=True, exist_ok=True)
    summary = {"gene": gene, "date": today(), "datasets": []}

    for ds in datasets:
        prism_dir = Path(ds["prism_dir"])
        peps, uniprot, pep_source = canonical_peptides_for_gene(prism_dir, gene, ds)
        if not peps:
            summary["datasets"].append(
                {"name": ds["name"], "status": "no canonical peptides found"})
            continue
        conds = build_sample_conditions(prism_dir, ds)
        sample_cols = list(conds.index)
        log2 = load_abundance(prism_dir, peps, sample_cols,
                              ds.get("peptide_matrix", "corrected_peptides.parquet"))
        if log2.empty:
            summary["datasets"].append(
                {"name": ds["name"], "status": "peptides not in matrix"})
            continue

        rows, ranked, subset_rows = [], [], []
        n_suggestive = 0
        for c in conds.columns:
            scan = peptide_condition_scan(log2, conds[c])
            if not scan.empty:
                scan.insert(0, "condition", c)
                rows.append(scan)
                n_suggestive += int(scan["suggestive"].sum())
                sc = score_condition(scan)
                sc["condition"] = c
                ranked.append(sc)
            ss = subset_outlier_scan(log2, conds[c])
            if not ss.empty:
                ss.insert(0, "condition", c)
                subset_rows.append(ss)

        detail = (pd.concat(rows, ignore_index=True) if rows
                  else pd.DataFrame(columns=["condition", "peptide"]))
        ranked_df = pd.DataFrame(ranked).sort_values(
            ["n_down", "median_effect"], ascending=[False, True]
        ) if ranked else pd.DataFrame()

        # subset "dive deeper" findings (some-patients signals)
        subset = (pd.concat(subset_rows, ignore_index=True) if subset_rows
                  else pd.DataFrame())
        subset_flags = (subset[subset["n_outliers"] >= SUBSET_MIN_OUTLIERS]
                        .sort_values("n_outliers", ascending=False)
                        if not subset.empty else pd.DataFrame())
        top_subset = None
        if not subset_flags.empty:
            r0 = subset_flags.iloc[0]
            top_subset = {"condition": str(r0["condition"]), "level": str(r0["level"]),
                          "peptide": str(r0["peptide"]), "n_outliers": int(r0["n_outliers"]),
                          "n_level": int(r0["n_level"]), "samples": str(r0["outlier_samples"])}

        # metadata-AGNOSTIC discovery: samples that coordinately drop many canonical
        # peptides at once (a subpopulation), found with no labels; annotate post-hoc.
        unsup, unsup_peps = unsupervised_discovery(log2)
        if not unsup.empty:
            # specificity control: is the drop protein-specific or globally low-input?
            try:
                grate = global_low_rate(prism_dir, ds, list(log2.columns))
                unsup["global_low_frac"] = grate.reindex(unsup.index).round(3)
                unsup["specificity"] = (unsup["concordance"]
                                        - unsup["global_low_frac"]).round(3)
                unsup["specific"] = unsup["specificity"] >= SPECIFICITY_MARGIN
                unsup = unsup.sort_values(["specific", "specificity", "n_low"],
                                          ascending=[False, False, False])
            except Exception:  # specificity is best-effort; keep raw discovery
                unsup["global_low_frac"] = np.nan
                unsup["specific"] = pd.NA
            annot = conds.reindex(unsup.index)
            annot.insert(0, "sample_name", [s.split("__@__")[0] for s in unsup.index])
            unsup = unsup.join(annot)
        n_specific = int(unsup["specific"].sum()) if "specific" in unsup else 0
        top_unsup = None
        if not unsup.empty:
            # headline the strongest protein-SPECIFIC sample if any, else the strongest overall
            spec = unsup[unsup.get("specific") == True] if "specific" in unsup else unsup
            pick = spec.iloc[0] if not spec.empty else unsup.iloc[0]
            desc = {c: str(pick[c]) for c in conds.columns if pd.notna(pick.get(c))}
            top_unsup = {"sample": str(pick.get("sample_name")),
                         "n_low": int(pick["n_low"]), "n_testable": int(pick["n_testable"]),
                         "concordance": float(pick["concordance"]),
                         "global_low_frac": (float(pick["global_low_frac"])
                                             if pd.notna(pick.get("global_low_frac"))
                                             else None),
                         "specific": bool(pick.get("specific"))
                         if pd.notna(pick.get("specific")) else None,
                         "annot": desc}

        # detection recovery for the single top-ranked (condition, level)
        det_frac = None
        top = None
        if not ranked_df.empty and ranked_df.iloc[0]["n_down"] > 0:
            top = ranked_df.iloc[0]
            if detection:
                try:
                    det_frac = detection_by_level(
                        prism_dir, ds, top["down_peptides"],
                        conds[top["condition"]], top["modal_level"])
                except Exception as e:  # detection is best-effort
                    det_frac = None
                    detail.attrs["detection_error"] = str(e)

        out_csv = runs_dir / f"{today()}_{gene}_{ds['name']}.csv"
        detail.to_csv(out_csv, index=False)
        dd_csv = None
        if not subset.empty:
            dd_csv = runs_dir / f"{today()}_{gene}_{ds['name']}_divedeeper.csv"
            subset.sort_values("n_outliers", ascending=False).to_csv(dd_csv, index=False)
        subpop_csv = None
        if not unsup.empty:
            subpop_csv = runs_dir / f"{today()}_{gene}_{ds['name']}_subpops.csv"
            unsup.to_csv(subpop_csv, index_label="sample_id")

        ds_summary = {
            "name": ds["name"], "status": "ok", "n_canonical_peptides": len(peps),
            "n_peptides_in_matrix": int(log2.shape[0]),
            "leading_uniprot": uniprot, "peptide_source": pep_source,
            "csv": str(out_csv),
            "top_condition": None if top is None else str(top["condition"]),
            "top_level": None if top is None else str(top["modal_level"]),
            "n_down": 0 if top is None else int(top["n_down"]),
            "median_effect_log2": None if top is None else top["median_effect"],
            "discordant": None if top is None else bool(top["discordant"]),
            "detected_frac_in_level": det_frac,
            "n_suggestive": int(n_suggestive),
            "n_subset_flags": int(len(subset_flags)),
            "top_subset": top_subset,
            "divedeeper_csv": None if dd_csv is None else str(dd_csv),
            "n_subpop_samples": int(len(unsup)),
            "n_subpop_specific": n_specific,
            "top_unsup": top_unsup,
            "subpop_csv": None if subpop_csv is None else str(subpop_csv),
        }
        summary["datasets"].append(ds_summary)
    return summary


# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #
LOG_HEADER = """# Cryptic-peptide analysis log

A running record of cryptic-peptide re-mining runs. Each block is one `mine` invocation.
Screens are triage, not claims: promote a hit with the moderated limma/DEqMS tool
(`cryptic-intensity-prior-limma`) before reporting.

| Date | Gene | Dataset | Top condition | Level | #down | median log2 | discordant | detect frac | dive deeper | CSV |
|------|------|---------|---------------|-------|-------|-------------|------------|-------------|-------------|-----|
"""


def log_run(db_dir: Path, summary: dict, source: str = "") -> Path:
    log_file = db_dir.parent / LOG_NAME
    if not log_file.exists():
        log_file.write_text(LOG_HEADER, encoding="utf-8")
    lines = []
    for ds in summary["datasets"]:
        if ds.get("status") != "ok":
            lines.append(f"| {summary['date']} | {summary['gene']} | "
                         f"{ds['name']} | _{ds.get('status')}_ |  |  |  |  |  |  |")
            continue
        eff = ds["median_effect_log2"]
        eff = f"{eff:.2f}" if isinstance(eff, (int, float)) and eff == eff else ""
        det = ds["detected_frac_in_level"]
        det = f"{det:.2f}" if isinstance(det, (int, float)) and det == det else ""
        csv_name = Path(ds["csv"]).name if ds.get("csv") else ""
        ts = ds.get("top_subset")
        parts = []
        if ds.get("n_subpop_samples"):
            parts.append(f"{ds['n_subpop_samples']} subpop samples")
        if ts:
            s = f"{ts['condition']}='{ts['level']}' {ts['n_outliers']}/{ts['n_level']}"
            if ds.get("n_subset_flags", 0) > 1:
                s += f" (+{ds['n_subset_flags'] - 1})"
            parts.append(s)
        elif ds.get("n_suggestive"):
            parts.append(f"{ds['n_suggestive']} suggestive")
        dig = "; ".join(parts)
        lines.append(
            f"| {summary['date']} | {summary['gene']} | {ds['name']} | "
            f"{ds.get('top_condition') or ''} | {ds.get('top_level') or ''} | "
            f"{ds.get('n_down', '')} | {eff} | {ds.get('discordant', '')} | "
            f"{det} | {dig} | {csv_name} |"
        )
    with open(log_file, "a", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return log_file


# --------------------------------------------------------------------------- #
# Run summaries -> plain-language findings.md
# --------------------------------------------------------------------------- #
def _json_safe(o):
    import math
    if isinstance(o, dict):
        return {k: _json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_safe(x) for x in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return None if math.isnan(f) else f
    if o is pd.NA:
        return None
    return o


def write_run_summaries(db_dir: Path, summary: dict) -> None:
    """Persist each dataset's structured result so findings.md can be regenerated."""
    import json
    sdir = db_dir / RUNS_SUBDIR / "summaries"
    sdir.mkdir(parents=True, exist_ok=True)
    for ds in summary["datasets"]:
        rec = {"gene": summary["gene"], "date": summary["date"], **ds}
        (sdir / f"{summary['gene']}__{ds['name']}.json").write_text(
            json.dumps(_json_safe(rec), indent=2), encoding="utf-8")


LOW_DETECT = 0.3  # a strict hit with detection below this is baseline noise, not signal


def _is_confound(cond) -> bool:
    """A condition axis that is technical (batch/plate/prep/acquisition), not biology."""
    return bool(cond) and bool(re.search(r"(?i)plate|batch|prep|\bacq", str(cond)))


def _real_bio(d: dict) -> bool:
    """A strict hit is a genuine biological lead only if not a confound and truly detected."""
    if not d.get("n_down") or _is_confound(d.get("top_condition")):
        return False
    det = d.get("detected_frac_in_level")
    if isinstance(det, (int, float)) and det < LOW_DETECT:
        return False
    return True


def _finding_sentence(d: dict) -> str:
    if d.get("status") != "ok":
        return "not detected here (no canonical peptides in this dataset)."
    parts = [f"{d['n_peptides_in_matrix']} canonical peptides detected."]
    if d.get("n_down"):
        cond, lvl = d.get("top_condition"), d.get("top_level")
        eff = d.get("median_effect_log2")
        eff_s = f"{eff:.2f} log2" if isinstance(eff, (int, float)) else "lower"
        conf = " (this is a technical axis - likely artifact, not disease)" \
            if _is_confound(cond) else ""
        det = d.get("detected_frac_in_level")
        det_s = (f" but only {det:.0%} truly detected, so mostly baseline noise"
                 if isinstance(det, (int, float)) and det < 0.3 else "")
        parts.append(f"STRICT: {d['n_down']} peptide(s) drop in {cond}='{lvl}' "
                     f"({eff_s}){conf}{det_s}.")
    else:
        parts.append("No group-level canonical decrease.")
    if d.get("n_subpop_specific"):
        tu = d.get("top_unsup") or {}
        parts.append(f"LABEL-FREE: {d['n_subpop_specific']} sample(s) specifically drop "
                     f"this protein (e.g. {tu.get('sample', '?')} low in "
                     f"{tu.get('n_low', '?')}/{tu.get('n_testable', '?')} peptides).")
    elif d.get("n_subpop_samples"):
        parts.append(f"{d['n_subpop_samples']} subpopulation candidate(s) but none passed "
                     f"the specificity control (likely just low-input samples).")
    return " ".join(parts)


def generate_findings(db_dir: Path) -> Path:
    """Aggregate all run summaries into a plain-language findings.md, grouped by target."""
    import json
    sdir = db_dir / RUNS_SUBDIR / "summaries"
    db = load_db(db_dir)
    src_by_gene = {}
    for _, r in db.iterrows():
        g = r["gene"]
        if g and g not in src_by_gene:
            s = (r.get("source", "") or "")
            if r.get("doi_pmid"):
                s = (s + f" ({r['doi_pmid']})").strip()
            if s:
                src_by_gene[g] = s
    recs = []
    if sdir.exists():
        for f in sorted(sdir.glob("*.json")):
            try:
                recs.append(json.loads(f.read_text(encoding="utf-8")))
            except Exception:
                pass
    by_gene: dict = {}
    for r in recs:
        by_gene.setdefault(r["gene"], []).append(r)

    # dataset order + tissue from the registry (fallback to observed order)
    ds_order, tissue = [], {}
    try:
        for d in load_registry(db_dir).get("datasets", []):
            ds_order.append(d["name"])
            tissue[d["name"]] = d.get("tissue", "")
    except Exception:
        pass
    seen = [d["name"] for r in recs for d in [r]]
    for r in recs:
        if r["name"] not in ds_order:
            ds_order.append(r["name"])

    def short_verdict(gene):
        dl = by_gene[gene]
        det = [d for d in dl if d.get("status") == "ok"]
        strict = [d for d in det if d.get("n_down")]
        real = [d for d in strict if _real_bio(d)]
        sub = [d for d in det if d.get("n_subpop_specific")]
        if not det:
            return "not measurable"
        if real:
            return f"signal ({len(real)} ds)"
        if sub:
            return f"subpopulation ({len(sub)} ds)"
        if strict:
            return "technical only"
        return "flat"

    def best_signal(gene):
        cands = [d for d in by_gene[gene] if d.get("status") == "ok" and _real_bio(d)]
        if cands:
            b = max(cands, key=lambda d: d["n_down"])
            return (f"{b['name']}: {b['top_condition']}='{b['top_level']}' "
                    f"({b['n_down']} pep, {b.get('median_effect_log2'):.2f})")
        subs = [d for d in by_gene[gene] if d.get("n_subpop_specific")]
        if subs:
            b = max(subs, key=lambda d: d["n_subpop_specific"])
            return f"{b['name']}: {b['n_subpop_specific']} subpop samples"
        return "-"

    out = [
        "# Cryptic re-mining findings",
        "",
        f"_Auto-generated {today()} from {len(recs)} target x dataset runs._",
        "",
        "Screens (leads), not claims - confirm with limma/DEqMS + batch covariate. "
        "STRICT = concordant canonical decrease; LABEL-FREE = specificity-checked "
        "subpopulation (no metadata).",
        "",
        "## Target verdicts",
        "",
        "| Target | Verdict | Detected | Best biological signal |",
        "|--------|---------|----------|------------------------|",
    ]
    for gene in sorted(by_gene):
        det = sum(1 for d in by_gene[gene] if d.get("status") == "ok")
        out.append(f"| {gene} | {short_verdict(gene)} | {det}/{len(by_gene[gene])} "
                   f"| {best_signal(gene)} |")

    # detection matrix: canonical peptides per (target x dataset)
    out += ["", "## Detection (canonical peptides per dataset)", "",
            "Cells = # canonical peptides in the matrix; `-` = not detected here.", ""]
    header = "| Dataset | Tissue | " + " | ".join(sorted(by_gene)) + " |"
    out.append(header)
    out.append("|" + "---|" * (2 + len(by_gene)))
    per_ds_detected = {}
    for dsn in ds_order:
        cells = []
        det_here = 0
        for gene in sorted(by_gene):
            rec = next((d for d in by_gene[gene] if d["name"] == dsn), None)
            if rec and rec.get("status") == "ok":
                cells.append(str(rec.get("n_peptides_in_matrix", "?")))
                det_here += 1
            else:
                cells.append("-")
        per_ds_detected[dsn] = det_here
        out.append(f"| {dsn} | {tissue.get(dsn, '')} | " + " | ".join(cells) + " |")
    out += ["",
            "_Detection rate per dataset (targets detected / queried):_ "
            + ", ".join(f"{d} {per_ds_detected[d]}/{len(by_gene)}" for d in ds_order
                        if d in per_ds_detected),
            ""]

    for gene in sorted(by_gene):
        ds_list = sorted(by_gene[gene], key=lambda d: d["name"])
        detected = [d for d in ds_list if d.get("status") == "ok"]
        strict = [d for d in detected if d.get("n_down")]
        real_strict = [d for d in strict if _real_bio(d)]
        subpop = [d for d in detected if d.get("n_subpop_specific")]
        if not detected:
            verdict = "Not measurable - no canonical peptides in any registered dataset."
        elif real_strict:
            verdict = (f"Possible signal - a canonical decrease tied to a biological "
                       f"condition in {len(real_strict)} dataset(s). Worth confirming.")
        elif subpop:
            verdict = (f"Label-free subpopulation in {len(subpop)} dataset(s), "
                       f"specificity-controlled - worth a look (may be sample-type/prep).")
        elif strict:
            verdict = "Only technical (plate/batch/prep) hits - no biological signal."
        else:
            verdict = "Detected but flat - no canonical decrease and no subpopulation."
        out.append(f"## {gene}")
        if gene in src_by_gene:
            out.append(f"*Reported by:* {src_by_gene[gene]}")
        out.append(f"**Verdict:** {verdict}")
        out.append(f"**Detected in** {len(detected)}/{len(ds_list)} datasets.")
        out.append("")
        for d in ds_list:
            out.append(f"- **{d['name']}** - {_finding_sentence(d)}")
        out.append("")
    path = db_dir / "findings.md"
    path.write_text("\n".join(out), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# FASTA generation -- append DB cryptic targets onto the latest search FASTA
# --------------------------------------------------------------------------- #
def _read_fasta_sequences(path: Path) -> set:
    """Return the set of (uppercased) sequences already present in a FASTA."""
    seqs, cur = set(), []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.add("".join(cur).upper())
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.add("".join(cur).upper())
    return seqs


def generate_fasta(db_dir: Path, base_fasta: Path, out_dir: Path | None = None) -> dict:
    """
    Write a new FASTA = the latest base search FASTA + every DB cryptic-peptide target
    whose sequence is not already present. Header uses the CRYPTIC| scheme so downstream
    tools (and this tool's seed-db) recognize the appended entries.
    """
    db = load_db(db_dir)
    have_seq = db[db["cryptic_peptide_sequence"].astype(str).str.strip().ne("")]
    n_noseq = int(len(db) - len(have_seq))
    base_seqs = _read_fasta_sequences(base_fasta)
    base_text = Path(base_fasta).read_text(encoding="utf-8", errors="replace")
    appended, skipped_present = [], 0
    for _, r in have_seq.iterrows():
        s = str(r["cryptic_peptide_sequence"]).upper()
        if s in base_seqs:
            skipped_present += 1
            continue
        src = r.get("source") or "cryptic-remine"
        hdr = f">CRYPTIC|{r['gene']}|{r['entry_id']}|{r['cryptic_type']}|source:{src}"
        appended.append((hdr, str(r["cryptic_peptide_sequence"])))
        base_seqs.add(s)
    out_dir = Path(out_dir) if out_dir else (db_dir / "fastas")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{Path(base_fasta).stem}_cryptic-remine_{today()}.fasta"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(base_text if base_text.endswith("\n") else base_text + "\n")
        for hdr, seq in appended:
            fh.write(f"{hdr}\n{seq}\n")
    return {"path": out, "n_appended": len(appended),
            "n_already_present": skipped_present, "n_no_sequence": n_noseq}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _print_summary(summary: dict) -> None:
    print(f"\n=== re-mine: {summary['gene']}  ({summary['date']}) ===")
    for ds in summary["datasets"]:
        if ds.get("status") != "ok":
            print(f"  [{ds['name']}] {ds.get('status')}")
            continue
        print(f"  [{ds['name']}] {ds['n_peptides_in_matrix']}/{ds['n_canonical_peptides']}"
              f" canonical peptides in matrix (leading {ds['leading_uniprot']}"
              f", via {ds['peptide_source']})")
        if ds["n_down"]:
            print(f"      STRICT: {ds['top_condition']} = '{ds['top_level']}' -> "
                  f"{ds['n_down']} peptides down (median {ds['median_effect_log2']:.2f} "
                  f"log2), discordant={ds['discordant']}, "
                  f"detect_frac={ds['detected_frac_in_level']}")
        else:
            print("      STRICT: no condition-linked canonical decrease passed")
        # metadata-agnostic discovery (headline "dive deeper" layer)
        tu = ds.get("top_unsup")
        if tu:
            desc = ", ".join(f"{k}={v}" for k, v in tu["annot"].items()) or "(no metadata)"
            spec = ds.get("n_subpop_specific", 0)
            print(f"      >> SUBPOPULATION (no labels used): {ds['n_subpop_samples']} samples "
                  f"drop many peptides at once; {spec} are protein-SPECIFIC")
            g = tu.get("global_low_frac")
            gtxt = (f" vs global {g:.2f} -> "
                    f"{'SPECIFIC' if tu.get('specific') else 'GLOBALLY-LOW (likely artifact)'}"
                    if g is not None else "")
            print(f"         top: {tu['sample']} low in {tu['n_low']}/{tu['n_testable']} "
                  f"peptides (conc {tu['concordance']}{gtxt}); metadata: {desc}")
        # liberal condition-based "dive deeper" tier
        ts = ds.get("top_subset")
        if ts:
            print(f"      >> DIVE DEEPER: {ts['n_outliers']}/{ts['n_level']} "
                  f"samples in {ts['condition']}='{ts['level']}' are low-outliers for "
                  f"peptide {ts['peptide']}")
            print(f"                   samples: {ts['samples']}")
        if ds.get("n_subset_flags"):
            print(f"                   ({ds['n_subset_flags']} subset flags, "
                  f"{ds['n_suggestive']} suggestive-trend peptides; see divedeeper CSV)")
        elif ds.get("n_suggestive"):
            print(f"      >> DIVE DEEPER: {ds['n_suggestive']} peptides show a "
                  f"suggestive downward trend (looser bar)")
        print(f"      CSV: {ds['csv']}")
        if ds.get("subpop_csv"):
            print(f"      subpopulation CSV: {ds['subpop_csv']}")
        if ds.get("divedeeper_csv"):
            print(f"      dive-deeper CSV: {ds['divedeeper_csv']}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR,
                    help=f"database directory (default {DEFAULT_DB_DIR})")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("seed-db", help="build curated DB from a search FASTA")
    sp.add_argument("--fasta", type=Path, required=True)
    sp.add_argument("--unique-reference", action="store_true",
                    help="also write the bulk CRYPTIC_UNIQUE catalog to "
                         "cryptic_unique_reference.csv (sidecar, not the main DB)")

    ap_add = sub.add_parser("add", help="register a cryptic-relevant protein")
    ap_add.add_argument("--gene", required=True)
    ap_add.add_argument("--sequence", default="")
    ap_add.add_argument("--source", default="")
    ap_add.add_argument("--doi", default="")
    ap_add.add_argument("--uniprot", default="")
    ap_add.add_argument("--notes", default="")

    ap_mine = sub.add_parser("mine", help="re-mine PRISM datasets for a gene")
    ap_mine.add_argument("--gene", required=True)
    ap_mine.add_argument("--dataset", default=None, help="registry name; omit = all")
    ap_mine.add_argument("--no-detection", action="store_true")
    ap_mine.add_argument("--no-log", action="store_true",
                         help="do not append to the analysis log")
    ap_mine.add_argument("--source", default="", help="paper/source, recorded in log")

    ap_conf = sub.add_parser("confirm",
                             help="limma-style moderated-t confirmation of a canonical "
                                  "decrease for a gene in one dataset (batch-adjusted)")
    ap_conf.add_argument("--gene", required=True)
    ap_conf.add_argument("--dataset", required=True, help="registry name")
    ap_conf.add_argument("--condition", required=True,
                         help="metadata column to contrast (e.g. Condition, CERAD, siRNA)")
    ap_conf.add_argument("--target", required=True,
                         help="level of --condition to test vs all other levels")
    ap_conf.add_argument("--covariate", default=None,
                         help="batch/plate column to adjust for (e.g. batch, Plate)")
    ap_conf.add_argument("--matrix", default="peptides_log2_internal.parquet",
                         help="peptide matrix (default: normalized log2, non-ComBat)")

    ap_sky = sub.add_parser("skyline-import",
                            help="build a mineable PRISM-shaped dir from a Skyline report "
                                 "(e.g. exported/pulled from Panorama)")
    ap_sky.add_argument("--report", required=True, help="Skyline report csv/tsv/parquet")
    ap_sky.add_argument("--out-dir", required=True, help="destination dir to register")
    ap_sky.add_argument("--format", default=None, choices=["csv", "parquet"])

    ap_find = sub.add_parser("findings", help="(re)generate the plain-language findings.md")

    ap_fa = sub.add_parser("make-fasta",
                           help="write latest FASTA + DB cryptic targets appended")
    ap_fa.add_argument("--base-fasta", type=Path, default=None,
                       help="base search FASTA (default: `fasta:` in datasets.yaml)")
    ap_fa.add_argument("--out-dir", type=Path, default=None)

    args = ap.parse_args(argv)

    if args.cmd == "seed-db":
        curated, uniq = seed_db_from_fasta(args.fasta)
        save_db(curated, args.db_dir)
        print(f"Seeded {len(curated)} curated cryptic entries -> {db_path(args.db_dir)}")
        print(curated["cryptic_type"].value_counts().to_string())
        print(f"\nGenes: {sorted(curated['gene'].unique())}")
        if args.unique_reference:
            ref = args.db_dir / "cryptic_unique_reference.csv"
            uniq.to_csv(ref, index=False)
            print(f"\nWrote {len(uniq)} CRYPTIC_UNIQUE reference rows -> {ref}")
        return 0

    if args.cmd == "add":
        df = add_entry(args.db_dir, args.gene, args.sequence, args.source,
                       args.doi, uniprot=args.uniprot, notes=args.notes)
        print(f"DB now has {len(df)} entries; {args.gene} registered.")
        _regenerate_fasta(args.db_dir)  # keep the target FASTA in sync with the DB
        return 0

    if args.cmd == "mine":
        summary = mine_gene(args.db_dir, args.gene, args.dataset,
                            detection=not args.no_detection)
        _print_summary(summary)
        write_run_summaries(args.db_dir, summary)
        if not args.no_log:
            lf = log_run(args.db_dir, summary, source=args.source)
            print(f"\nLogged to {lf}")
        fpath = generate_findings(args.db_dir)  # refresh plain-language findings
        print(f"Findings updated: {fpath}")
        _regenerate_leaderboard(args.db_dir)  # keep LEADERBOARD.md in sync
        return 0

    if args.cmd == "skyline-import":
        from skyline_adapter import import_skyline_report
        r = import_skyline_report(args.report, args.out_dir, fmt=args.format)
        print(f"imported {r['transitions']:,} transitions -> {r['peptides']:,} peptides x "
              f"{r['replicates']} replicates, {r['genes']} genes  ({r['out_dir']})")
        print("Register in datasets.yaml:  prism_dir: <out-dir>  "
              "peptide_matrix: peptides.parquet  merged_data: merged_data.parquet")
        return 0

    if args.cmd == "confirm":
        from limma_confirm import run_confirm
        run_confirm(args.db_dir, args.dataset, args.gene, args.condition, args.target,
                    covariate=args.covariate, matrix=args.matrix,
                    outdir=str(args.db_dir / "remine_runs" / "limma"))
        return 0

    if args.cmd == "findings":
        print(f"Findings written: {generate_findings(args.db_dir)}")
        _regenerate_leaderboard(args.db_dir)
        return 0

    if args.cmd == "make-fasta":
        base = args.base_fasta
        if base is None:
            base = load_registry(args.db_dir).get("fasta")
            if not base:
                print("No --base-fasta given and no `fasta:` in datasets.yaml.")
                return 1
        res = generate_fasta(args.db_dir, Path(base), args.out_dir)
        print(f"Wrote {res['path']}")
        print(f"  appended {res['n_appended']} DB targets; "
              f"{res['n_already_present']} already in base; "
              f"{res['n_no_sequence']} DB entries have no sequence (skipped).")
        return 0

    return 1


def _regenerate_leaderboard(db_dir: Path) -> None:
    """Best-effort LEADERBOARD.md refresh (mirrors the findings.md auto-refresh)."""
    try:
        import importlib.util
        mk = Path(__file__).with_name("make_leaderboard.py")
        spec = importlib.util.spec_from_file_location("make_leaderboard", mk)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        by = mod.load(str(db_dir))
        out = Path(db_dir) / "LEADERBOARD.md"
        out.write_text(mod.build(by, _dt.date.today().isoformat()), encoding="utf-8")
        print(f"Leaderboard updated: {out}")
    except Exception as e:
        print(f"(Leaderboard not refreshed: {e})")


def _regenerate_fasta(db_dir: Path) -> None:
    """Best-effort FASTA refresh after a DB change (uses registry `fasta:` as base)."""
    try:
        base = load_registry(db_dir).get("fasta")
        if not base:
            return
        res = generate_fasta(db_dir, Path(base))
        print(f"Target FASTA refreshed: {res['path']} "
              f"(+{res['n_appended']} appended, {res['n_no_sequence']} lack a sequence)")
    except Exception as e:
        print(f"(FASTA not refreshed: {e})")


if __name__ == "__main__":
    sys.exit(main())
