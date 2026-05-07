"""Convert DIA-NN library parquet to OpenSWATH-assay TSV that BlibBuild can read.

Column names match the tokens the BlibBuild binary searches for in OpenSWATH-assay
files (PrecursorMz, ProductMz, LibraryIntensity, NormalizedRetentionTime,
ModifiedPeptideSequence, StrippedPeptide, PrecursorCharge, FragmentCharge,
FragmentType, FragmentSeriesNumber, FragmentLossType, ProteinId/ProteinName,
Decoy, IonMobility).
"""
import argparse
import os

import pandas as pd
import pyarrow.parquet as pq


def convert(src: str, dst: str) -> None:
    df = pq.read_table(src).to_pandas()
    print(f"read {len(df):,} rows from {src}")

    out = pd.DataFrame({
        "PrecursorMz":             df["Precursor.Mz"],
        "ProductMz":                df["Product.Mz"],
        "LibraryIntensity":         df["Relative.Intensity"],
        "NormalizedRetentionTime":  df["RT"],
        "IonMobility":              df["IM"],
        "StrippedPeptide":          df["Stripped.Sequence"],
        "ModifiedPeptideSequence":  df["Modified.Sequence"],
        "PrecursorCharge":          df["Precursor.Charge"],
        "FragmentCharge":           df["Fragment.Charge"],
        "FragmentType":             df["Fragment.Type"],
        "FragmentSeriesNumber":     df["Fragment.Series.Number"],
        "FragmentLossType":         df["Fragment.Loss.Type"],
        "ProteinId":                df["Protein.Ids"],
        "ProteinName":              df["Protein.Names"],
        "GeneName":                 df["Genes"],
        "Decoy":                    df["Decoy"].astype(int),
    })

    out.to_csv(dst, sep="\t", index=False)
    print(f"wrote {dst} ({os.path.getsize(dst):,} bytes)")


def _default_dst(src: str) -> str:
    if src.lower().endswith(".parquet"):
        return src[: -len(".parquet")] + ".tsv"
    return src + ".tsv"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a DIA-NN spectral-library .parquet (report-style dotted "
            "column names) into an OpenSWATH-assay TSV that Skyline's BlibBuild "
            "can ingest directly."
        )
    )
    parser.add_argument("src", help="Path to the DIA-NN library .parquet file")
    parser.add_argument(
        "-o", "--out",
        default=None,
        help="Output TSV path (default: <src> with .parquet replaced by .tsv)",
    )
    args = parser.parse_args()
    dst = args.out if args.out is not None else _default_dst(args.src)
    convert(args.src, dst)


if __name__ == "__main__":
    main()
