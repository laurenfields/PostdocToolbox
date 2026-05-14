# diann_lib_add_silac_heavy

Duplicate a DIA-NN spectral library (`.parquet`) with SILAC heavy-labeled
versions of every precursor so DIA-NN can search SILAC light/heavy pairs
against a library that was originally generated without heavy mods.

## Why this exists

Tools like carafe / DIA-NN's `--gen-spec-lib` produce a library with only
the fixed mods configured at generation time (typically just Carbamidomethyl
C). For a SILAC experiment, every K- and R-containing precursor needs a
heavy counterpart in the library. Rather than re-running carafe / DIA-NN
library generation with SILAC mods declared, this script appends the heavy
twin of each existing precursor directly to the parquet — much faster and
keeps fragment intensities/RTs from the carafe model intact.

## What it does

Default labels (configurable):

- `K` → Label:13C(6)15N(2), +8.014199 Da, `(UniMod:259)`
- `R` → Label:13C(6)15N(4), +10.008269 Da, `(UniMod:267)`

For each existing (light) precursor, the script adds a heavy duplicate:

- `Modified.Sequence`: inserts the UniMod tag after every K and every R
  (covers internal K/R from missed cleavages, not just the C-terminal one).
- `Precursor.Id`: rebuilt as `Modified.Sequence + Precursor.Charge`.
- `Precursor.Mz`: shifted by Σ(K/R deltas) / charge.
- `Product.Mz`: each b- or y-ion is shifted by the delta of whichever K/R
  residues *that specific fragment* contains, divided by `Fragment.Charge`.
  For pure C-terminal-K/R tryptic peptides this means b-ions are unchanged
  and y-ions all shift by the single label / charge. For missed-cleavage
  peptides, b-ions that cross an internal K/R are also shifted correctly.

Peptides with zero K and zero R are skipped — the heavy duplicate would be
identical to the light and would collide on `Precursor.Id`.

Output is a parquet file with the original light precursors plus the heavy
duplicates appended. Column order and dtypes are preserved so DIA-NN reads
it directly via `--lib`.

## Assumptions

- Input is a DIA-NN-format parquet library with columns including
  `Precursor.Id`, `Modified.Sequence`, `Stripped.Sequence`, `Precursor.Mz`,
  `Product.Mz`, `Fragment.Type` (`b`/`y`), `Fragment.Series.Number`,
  `Fragment.Charge`.
- The only existing modifications in `Modified.Sequence` are formatted as
  `<aa>(UniMod:<n>)` and never sit on a K or R — so naive
  `replace("K", "K(UniMod:259)")` is safe. This holds for carafe / DIA-NN
  libraries built with Carbamidomethyl C as the only fixed mod. If you
  plan to run this on a library that already carries a K/R modification,
  swap the `.str.replace` calls for something position-aware first.

## Dependencies

- `pyarrow`
- `pandas`
- `numpy`

```
pip install pyarrow pandas numpy
```

## Usage

```
python add_silac_heavy.py <library>.parquet
python add_silac_heavy.py <library>.parquet -o report-lib.heavy-light.parquet
```

Custom labels (e.g. mTRAQ, dimethyl, Lys6/Arg6):

```
python add_silac_heavy.py <library>.parquet \
  --k-delta 6.020129 --k-mod "(UniMod:188)" \
  --r-delta 6.020129 --r-mod "(UniMod:188)"
```

## Re-running DIA-NN with the output

```
diann.exe --lib report-lib.heavy-light.parquet \
          --f <your.mzML> --out report.parquet \
          --threads N --qvalue 0.01 --matrices --reanalyse ...
```

Drop the original `--fasta-search`, `--gen-spec-lib`, `--predictor`, and
`--unimod4` flags — the fixed mod is already baked into the library, and
you no longer want DIA-NN to regenerate. If your original `--max-pr-mz` is
close to the carafe library's upper limit, raise it slightly so heavy
precursors don't fall off the top edge.
