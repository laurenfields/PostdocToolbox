# known_ND_markers

A curated panel of **known neurodegeneration biomarkers** (ALS / FTD / AD), with
literature provenance. Maintained as a single CSV so it can be extended over time
and pulled directly (raw URL) by analysis pipelines.

## File
`known_ND_markers.csv` — one row per literature observation. Columns:

| Column | Meaning |
|---|---|
| `Gene` | **HGNC gene symbol** — the only field used for matching to proteomics data. Fill this for every row. |
| `Protein` | Name as cited in the source (e.g. "TDP-43", "NfL") — reference only. |
| `Cryptic?` | `Y` if the marker is a cryptic/non-canonical product, else `N`/blank. |
| `Location` | Biofluid / compartment reported (e.g. CSF, Plasma EVs). |
| `Notes` | Disease context (e.g. ALS, ALS-FTD, AD). |
| `DOI` | Source citation (the DOI URL is auto-extracted by consumers). |

## To add a marker
Append a row and fill **`Gene`** with its HGNC symbol (the rest is optional metadata).
No code changes are needed in any consumer.

## Use from an analysis (raw URL)
```python
url = "https://raw.githubusercontent.com/laurenfields/PostdocToolbox/main/proteomics/known_ND_markers/known_ND_markers.csv"
import pandas as pd
panel = pd.read_csv(url)
```
The EV/cryptic-peptide reporting pipeline reads this via `load_nd_markers(url, gene_col="Gene", fallback=<local copy>)`.
