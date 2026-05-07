# diann_lib_to_openswath_tsv

Convert a DIA-NN spectral-library `.parquet` into an OpenSWATH-assay TSV that
Skyline's BlibBuild can read directly.

## Why this exists

DIA-NN's `.parquet.skyline.speclib` is a binary that requires a paired
report-parquet with a `Global.Q.Value` column. When DIA-NN is run with
`--gen-spec-lib --reanalyse` and no raw input files (e.g. just regenerating
a library), the companion parquet is a *library* parquet and lacks
`Global.Q.Value`, so BlibBuild fails with:

```
ERROR: file <name>.parquet does not have a column called Global.Q.Value
```

Skyline also rejects DIA-NN's raw library TSV format with:

```
Only OpenSWATH result, OpenSWATH assay, AlphaPepDeep, and Paser .tsv files
are supported
```

This script bypasses both: it remaps the DIA-NN library parquet's columns
into OpenSWATH-assay format, which BlibBuild ingests directly with no paired
parquet needed.

## Important: column-name compatibility

The OpenSWATH-assay column names BlibBuild looks for are exact and
case-sensitive. They were determined by string-grepping `BlibBuild.exe`'s
binary in Skyline 26.1.1.097.

Required output columns:

```
PrecursorMz, ProductMz, LibraryIntensity, NormalizedRetentionTime,
IonMobility, StrippedPeptide, ModifiedPeptideSequence, PrecursorCharge,
FragmentCharge, FragmentType, FragmentSeriesNumber, FragmentLossType,
ProteinId, ProteinName, GeneName, Decoy
```

Names like `Tr_recalibrated`, `FullUniModPeptideName`, `PeptideSequence`,
`transition_group_id` are **not** recognized by this BlibBuild version.

## Dependencies

- `pyarrow`
- `pandas`

```
pip install pyarrow pandas
```

## Usage

```
python convert_to_openswath.py <library>.parquet
python convert_to_openswath.py <library>.parquet -o assay.tsv
```

If `--out` / `-o` is omitted, the output path is the input with `.parquet`
replaced by `.tsv`.

## Importing into Skyline

1. In Skyline: **Settings -> Peptide Settings -> Library -> Build...**
2. Pick the generated `.tsv` as the input.
3. BlibBuild produces a `.blib` you can add to the document library list.

## Tested versions

- Skyline 26.1.1.097-922725ca01
- DIA-NN 2.2.0 spectral-library parquet

On a representative library, BlibBuild reads 22,721 PSMs at
`score_threshold=0.01` and produces a valid `.blib`.

## Known limitations

- Column names are tied to a specific BlibBuild build. Future Skyline
  versions may rename or add required tokens; if BlibBuild rejects the
  output, re-grep `BlibBuild.exe` for the current set.
- Assumes the DIA-NN parquet uses the standard report-style dotted column
  names (`Precursor.Mz`, `Modified.Sequence`, `Fragment.Type`, etc.).
  Custom DIA-NN builds with renamed columns will need the mapping in
  `convert_to_openswath.py` updated.
- The `Decoy` column is coerced to int; if your parquet uses a different
  encoding for decoy state, adjust accordingly.
