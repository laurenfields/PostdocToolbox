# cryptic-remine

Log, database, and **canonical re-mining** for cryptic-peptide analyses of PRISM
output. When a paper reports that protein *X* produces a cryptic peptide, register
that context and re-mine already-completed PRISM searches for a **condition-linked
decrease in X's canonical peptides** that was previously overlooked.

## The idea

Producing a cryptic peptide changes the **stoichiometry** of the parent protein. A
cryptic/skiptic splicing event (a TDP-43 loss-of-function phenomenon) alters the
transcript, so the **canonical** protein — or specifically the canonical peptides
spanning/adjacent to the affected region — **decreases** in the affected condition. We
always observe canonical tryptic peptides in a PRISM search, so that decrease is already
sitting in completed results; without the context that *X* is cryptic-relevant it looks
like noise.

Key consequence: **the reported cryptic peptide does not need to be in your search
FASTA.** The database entry is *context* that points the miner at the canonical protein.
(If the cryptic peptide happens to be in your results, that enables an optional two-sided
check — cryptic up while canonical goes down.)

**This is a screen (triage), not a claim.** Promote a hit with the moderated
limma/DEqMS tool (`../cryptic-intensity-prior-limma`), modeling batch as a covariate,
before reporting anything.

## Install / requirements

Python ≥ 3.11 with `pandas`, `pyarrow`, `scipy`, `pyyaml`. Self-contained — reads PRISM
parquet directly, no cross-repo imports.

```bash
python cryptic_remine.py --help
```

The database directory defaults to `G:\Manuscripts\LatimerCryptic\cryptic_db` and can be
overridden with `--db-dir` or the `CRYPTIC_DB_DIR` environment variable.

## Files it manages

| File | What it is |
|------|------------|
| `cryptic_db/cryptic_peptides.csv` | Curated cryptic-peptide database (context registry). |
| `cryptic_db/cryptic_unique_reference.csv` | Optional bulk Seddighi DataS2 unique-peptide catalog (sidecar). |
| `cryptic_db/datasets.yaml` | Registry of completed PRISM runs to re-mine. |
| `cryptic_db/remine_runs/<date>_<gene>_<dataset>.csv` | Per-run detailed output (+ `_divedeeper.csv`, `_subpops.csv`). |
| `cryptic_db/remine_runs/summaries/<gene>__<dataset>.json` | Structured result per run; the source for `findings.md`. |
| `cryptic_db/findings.md` | **Plain-language** summary of every target, one section per gene. Regenerated after each `mine`. |
| `cryptic_db/LEADERBOARD.md` | **Ranked** cross-target view: KD-direction (cryptic-consistent) leaderboard, disease/plasma-stratum signal, subpopulation hits, verdict tally. Regenerated after each `mine` (or `make_leaderboard.py`). |
| `cryptic_db/fastas/<base>_cryptic-remine_<date>.fasta` | Latest search FASTA + DB cryptic targets appended. Regenerated after each `add`. |
| `cryptic_analysis_log.md` | Human-readable running log (one row per run). |

## Usage

### 1. Seed the database from a search FASTA (once)
Parses the custom cryptic headers (`CRYPTIC|`, `CRYPTIC_EXON|`, `ALT_EXON|`) into the
curated DB. The huge machine-generated `CRYPTIC_UNIQUE|` set (Seddighi DataS2) is written
to a separate sidecar only with `--unique-reference`.
```bash
python cryptic_remine.py seed-db --fasta path/to/search.fasta --unique-reference
```

### 2. Register a newly-reported cryptic protein
```bash
python cryptic_remine.py add --gene HDGFL2 --sequence PEPTIDESEQ --source "Author et al. 2026 Journal" --doi 10.xxxx/yyyy
```
Gene-only is allowed (`--gene X`) — it just records X as cryptic-relevant.

### 2b. Plain-language findings + an updated search FASTA

Two artifacts regenerate automatically, and can be rebuilt on demand:

```bash
python cryptic_remine.py findings     # rebuild cryptic_db/findings.md from all runs
python cryptic_remine.py make-fasta   # base search FASTA + DB cryptic targets appended
```

- **`findings.md`** — one plain-English section per target: a verdict, how many datasets
  detected it, and a sentence per dataset (strict hit / label-free subpopulation / negative,
  with confounds called out). Rebuilt after every `mine`.
- **`make-fasta`** — writes a new dated FASTA = the latest base search FASTA (`fasta:` in
  `datasets.yaml`, or `--base-fasta`) with every DB cryptic peptide **that has a sequence**
  appended under a `CRYPTIC|` header. Rebuilt after every `add`. To have a newly-dropped
  target appear in this FASTA (so the *next* DIA-NN/Skyline search can detect the cryptic
  peptide itself), give its sequence when adding: `add --gene X --sequence PEPTIDER ...`.
  Gene-only entries are catalogued but can't be appended (nothing to search for).

### 3. Re-mine completed PRISM datasets
```bash
python cryptic_remine.py mine --gene HDGFL2               # all registered datasets
python cryptic_remine.py mine --gene HDGFL2 --dataset plasma_ev_prism_20260629
python cryptic_remine.py mine --gene HDGFL2 --no-detection  # skip the merged_data pass (faster)
```
For each dataset the miner:
1. Resolves the gene's **canonical** peptides (via `protein_groups.csv`, falling back to
   transition-level `merged_data.ProteinGene` for genes that never formed a group).
2. Scans **every** registered clinical stratum (`condition_columns`) for a
   condition-associated decrease, using a robust Kruskal–Wallis test + a log2 effect
   (affected level vs the rest), BH-corrected across the protein's peptides. Strata are
   **ranked** so the strongest overlooked decrease surfaces first.
3. Flags **within-protein discordance** — a subset of canonical peptides dropping while
   others hold steady, the stoichiometry fingerprint of a region-specific cryptic event.
4. **Recovers true detection** from `merged_data` (`Truncated=false`,
   `DetectionQValue ≤ 0.01`, `Area > 0`) so a "decrease" into integrated baseline noise
   is distinguished from a genuine on/off. (Dense PRISM matrices carry no missingness
   signal, so abundance alone can mislead for absent peptides.)
5. Writes a detail CSV and appends a summary to the analysis log (unless `--no-log`).

### 4. Confirm a hit with limma (batch-adjusted moderated-t)
```bash
python cryptic_remine.py confirm --gene KALRN --dataset ineuron_sirna_tdp43kd --condition siRNA --target TDP43_KD
python cryptic_remine.py confirm --gene AGRN  --dataset fdr_goterm_ad_pd --condition Condition --target PDCN --covariate batch
```
Fits a limma-style moderated-t (empirical-Bayes variance shrinkage, Smyth 2004, implemented
in-package — no R needed) for `condition==target` vs all other levels, adjusting for
`--covariate` (batch/plate) as fixed effects. Prints the **batch × condition balance table**
first (a confounded contrast is flagged before you trust the p), then per-gene **Simes** p +
per-peptide BH; writes a peptide CSV to `remine_runs/limma/`. Default matrix is the normalized
non-ComBat `peptides_log2_internal.parquet` (falls back to the dataset's `peptide_matrix`).
*Note:* treats each matrix column as independent — correct for one-injection-per-subject
datasets; collapse technical injection replicates first for designs that have them.

### 5. Mine data from Panorama / any Skyline document
Point the pipeline at Skyline output (e.g. a collaborator's public Panorama folder) without
re-running PRISM. Export/pull a Skyline **transition-level report**, then:
```bash
# optional: pull from Panorama over WebDAV (or just Map Network Drive to the _webdav URL)
python panorama_pull.py get "https://panoramaweb.org/_webdav/Lab/Study/@files/report.csv" ./report.csv --apikey $PANORAMA_APIKEY

# convert the Skyline report -> a mineable PRISM-shaped dir
python cryptic_remine.py skyline-import --report ./report.csv --out-dir ./study_prism
```
`skyline-import` rolls transitions up to a `peptides.parquet` matrix (PRISM-style area sum),
writes `merged_data.parquet` (for detection), `protein_groups.csv`, and `sample_metadata.csv`.
Register the out-dir in `datasets.yaml` with `peptide_matrix: peptides.parquet` and
`merged_data: merged_data.parquet` (+ `name_conditions` or a clinical join for the contrast),
then `mine` / `confirm` it like any PRISM dataset. Report columns are auto-detected
(modern Skyline `PeptideModifiedSequenceUnimodIds` / `ProteinGene` / `Area` / `ReplicateName`,
or older `Peptide Modified Sequence` / `Protein` — gene falls back to protein when absent).

**Getting a report out of a `.sky.zip` (Panorama docs) — validated recipe.** Panorama publishes
Skyline documents, not reports. To export one headlessly with `SkylineCmd.exe`:
```bash
# 1. SkylineCmd can't open a .sky.zip directly -> unzip first, point --in at the .sky
unzip doc.sky.zip -d doc_extracted
# 2. export the BUILT-IN "Transition Results" report (has Peptide/Protein/Replicate/Area)
SkylineCmd --in="doc_extracted/doc.sky" --report-name="Transition Results" \
           --report-file=report.csv --report-format=CSV --report-invariant
# 3. that report's peptide column is bare "Peptide" -> rename it so the adapter sees it
#    (one-line header swap: Peptide -> "Peptide Modified Sequence"), then skyline-import
```
Gotcha: the ClickOnce `SkylineCmd.exe` lives under `AppData\Local\Apps\2.0\...`; pick the newest
folder that ALSO contains `Skyline.exe`/`Skyline-daily.exe` (not merely the newest `SkylineCmd.exe`).
Or, if a Skyline instance is open, pull the report live via the Skyline MCP instead of `SkylineCmd`.

### 6. Within-protein peptide discordance (cryptic/proteoform stoichiometry)
```bash
# scan one dataset: does any peptide's disease response differ from its protein's siblings?
python cryptic_remine.py discordance --dataset huad_smtg_cleandx --condition "Cognitive Status" --target Dementia
# after scanning several, build the cross-dataset reproducible table
python cryptic_remine.py discordance --aggregate --ad-regions huad_smtg_cleandx,huad_ipl_cleandx,huad_hipp_cleandx
```
A **PeCorA-style condition x peptide interaction** (batch-adjusted) per peptide: a peptide that
drops (or rises) unlike the rest of its protein is the stoichiometry footprint of a cryptic exon
/ proteoform — detectable even when the cryptic peptide itself was never searched. Because dense
(Skyline-imputed) matrices make single-dataset hits unreliable, the deliverable is
`discordance/REPRODUCIBLE_discordant_peptides.csv` — peptides discordant in the **same direction
across >= 2 datasets**. (Validated: independently recovers the amyloid-beta APP peptide going up
discordantly across AD regions.) Whole-protein losses (neuronal death, NMD-coupled cryptics) are
concordant and correctly do **not** appear here — use `mine`/`confirm` for those.

## Two tiers: strict vs dive-deeper

Every `mine` reports two complementary layers:

- **STRICT** — the conservative screen: a *concordant* canonical decrease across a
  condition level (robust KW test, effect ≤ −1 log2, BH-adjusted). Few false positives;
  this is what ranks datasets.
- **DIVE DEEPER** (liberal) — deliberately generous leads the strict test averages away:
  1. **Subset outliers** ("some patients") — individual samples that are strong low-outliers
     (robust z ≤ −2.5 vs the peptide's own median/MAD), broken down by condition level. This
     catches a signal in *part* of a group — e.g. some ALS but not PD patients dropping a
     given ACTN1 peptide — which a group-median test cannot see. Written to a companion
     `*_divedeeper.csv` (per peptide × level: `n_outliers`, `n_level`, `outlier_frac`, and the
     exact sample IDs).
  2. **Suggestive trends** — whole-group decreases at a looser, uncorrected bar
     (effect ≤ −0.6, nominal p ≤ 0.05).

Dive-deeper hits are **leads, not findings** — the point is to not miss things. Confirm any
lead with the strict path + limma before believing it.

## Metadata-agnostic subpopulation discovery

The most important layer needs **no labels at all** — so it can surface subgroups the
metadata can't even describe. For each canonical peptide it computes a robust z per sample;
a sample that is a low-outlier across **many** of the protein's peptides at once is a
candidate **subpopulation with reduced canonical protein** (`n_low`/`n_testable` =
concordance). Runs on every dataset, including label-less ones. Written to `*_subpops.csv`.

**Specificity control (built in):** a sample low across many of a protein's peptides could
just be globally low-input. So each candidate's low-rate for *this* protein is compared to
its low-rate across a random proteome background (`global_low_frac`); a sample is only
`specific=True` when the protein drop exceeds its global drop by `SPECIFICITY_MARGIN` (0.25).
Globally-low samples are down-ranked as loading artifacts, not biology — without using any
metadata to do it.

Whatever metadata exists is then attached to the discovered samples as **description only**
(never an input), so you can notice post-hoc structure (e.g. the subgroup is enriched for
control/prep-variant samples) after the data has spoken.

## Interpreting output

- **Top condition / level / #down / median log2**: the stratum with the strongest
  concordant canonical decrease, how many peptides dropped, and by how much (negative =
  down). `discordant=True` means some peptides dropped while others held — the specific,
  interesting pattern.
- **detect frac**: fraction of the down peptides still truly detected in the affected
  level. Near 1.0 → a real abundance decrease; near 0 → the peptides are essentially
  absent there (an on/off signal, possibly stronger but check it isn't all-baseline).
- The detail CSV has one row per (peptide × condition) with `effect_log2`, `direction`,
  `kw_p`, BH-adjusted `kw_p_bh`, and a `notable` flag.

## Caveats (read before trusting a hit)

- A hit is a **screen**, not a result. Confirm with limma/DEqMS + batch covariate.
- **Plate confound (important).** The plasma PRISM run is a **single batch** (ComBat did
  not run), and a prior QC-calibrated analysis found the **PLT1/PLT2 plate** drives ~33%
  of panel-peptide variance and that within-protein discordance did **not** track
  TDP-43/LATE. `Plate` is therefore scanned as a condition: **if `Plate` tops the ranking
  for a gene, the hit is a plate artifact, not biology.** The screen does not yet
  deconfound plate — a real biological hit must beat plate and survive limma with a plate
  covariate.
- Clinical columns are used verbatim; some (e.g. `Thal`) may contain Excel-mangled values
  — sanity-check the affected level.
- Technical injections (QC / pool / batch-reference) are excluded via the registry's
  `patient_filter_*` keys; verify these for a new dataset.
- Many neuronal cryptic targets (HDGFL2, UNC13A, ELAVL3, …) are **not detected in plasma
  EVs** at all — the miner truthfully reports "no canonical peptides found".

## Adding a new dataset

Append an entry to `cryptic_db/datasets.yaml`. Only `name` and `prism_dir` are required;
everything else is optional and combinable. The miner is **schema-tolerant** — it auto-
detects modern (`PeptideModifiedSequenceUnimodIds`, `ProteinGene`) vs older
(`Peptide Modified Sequence`, `Protein Gene`) PRISM column names, works with or without
`protein_groups.csv` (falls back to `merged_data`), and reconstructs the sample list from
the matrix when there's no `sample_metadata.csv`.

Registry keys:

| Key | Meaning |
|-----|---------|
| `prism_dir` | PRISM output folder (required). |
| `clinical_metadata`, `clinical_join_key`, `condition_columns` | Join a clinical CSV on a key (matched to the base sample name) and scan those columns. |
| `name_conditions` | Derive condition columns from the **sample name** by ordered regex rules — for datasets whose condition lives in the filename (siRNA target, ALS/HC, prep method). `{cond: [{label, pattern}, ...]}`; first matching rule wins. |
| `scan_batch` | Scan `sample_metadata.batch` as a confound-check condition (default `true`; only added when >1 batch). |
| `exclude_sample_regex` | Drop samples whose name matches (QC/pool/reference injections). |
| `patient_filter_column`, `patient_filter_exclude` | Drop samples whose value in a **clinical** column matches (technical labels). |

Three worked patterns are in the shipped registry: a clinical-CSV cohort (Latimer plasma;
Yubin CSF EV), a name-encoded disease cohort (Yubin CSF `dx: ALS/PD/HC/HD`), and a
name-encoded treatment/method (iNeuron `siRNA: TDP43_KD/control`; tALS `prep:`).

**Validation.** On the iNeuron TDP-43-knockdown data, `mine --gene STMN2` recovers the
textbook result — 11 canonical STMN2 peptides down ~8-fold (−3.05 log2) in `TDP43_KD`,
discordant — the canonical collapse caused by cryptic exon 2a inclusion. That is the
positive control for this whole approach.

**Note on method-comparison datasets** (tALS Lab_Comparison): these have no disease axis —
their conditions are EV-isolation methods, so a "hit" reflects prep, not biology.
