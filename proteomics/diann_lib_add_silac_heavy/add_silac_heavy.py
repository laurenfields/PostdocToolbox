"""
Duplicate a DIA-NN spectral library (.parquet) with SILAC heavy-labeled
versions of every precursor.

Default labels (configurable via --k-mod / --r-mod / --k-delta / --r-delta):
  K -> Label:13C(6)15N(2)  +8.014199 Da  (UniMod:259)
  R -> Label:13C(6)15N(4)  +10.008269 Da (UniMod:267)

Every K and every R is labeled (metabolic SILAC labels all K/R, including
internal residues from missed cleavages). For each labeled precursor the
script:

  - inserts the UniMod tag after every K and R in Modified.Sequence
  - rebuilds Precursor.Id as Modified.Sequence + Precursor.Charge
  - shifts Precursor.Mz by sum(deltas) / charge
  - shifts each Product.Mz by the delta of whichever K/R residues that b- or
    y-ion actually contains, divided by Fragment.Charge

Peptides with zero K and zero R are skipped (the heavy duplicate would be
identical to the light precursor and would collide on Precursor.Id).

Output is a parquet file with the original light precursors plus the heavy
duplicates appended, preserving the source library's column order and dtypes.
DIA-NN reads it directly via --lib.

Tested against carafe-generated DIA-NN 2.x libraries. The script makes the
assumption (true for carafe / DIA-NN parquet libraries) that the only mods
present in Modified.Sequence appear as `<aa>(UniMod:<n>)` and that K and R
never carry a pre-existing modification — i.e. naive
`replace("K", "K(UniMod:259)")` is safe. If you run this on a library that
already contains a K/R modification, adjust the insertion logic first.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def build_delta_tables(unique_seqs: pd.Series, delta_k: float, delta_r: float):
    records = []
    total_deltas: dict[str, float] = {}
    for seq in unique_seqs:
        L = len(seq)
        b = np.zeros(L + 1)
        y = np.zeros(L + 1)
        for i, aa in enumerate(seq):
            d = delta_k if aa == "K" else delta_r if aa == "R" else 0.0
            b[i + 1] = b[i] + d
        for i, aa in enumerate(reversed(seq)):
            d = delta_k if aa == "K" else delta_r if aa == "R" else 0.0
            y[i + 1] = y[i] + d
        total_deltas[seq] = float(b[L])
        for n in range(1, L + 1):
            records.append((seq, "b", n, float(b[n])))
            records.append((seq, "y", n, float(y[n])))
    delta_df = pd.DataFrame(
        records,
        columns=["Stripped.Sequence", "Fragment.Type", "Fragment.Series.Number", "frag_delta"],
    )
    return delta_df, total_deltas


def make_heavy(
    df: pd.DataFrame,
    *,
    delta_k: float,
    delta_r: float,
    mod_k: str,
    mod_r: str,
) -> pd.DataFrame:
    unique_seqs = df["Stripped.Sequence"].drop_duplicates()
    delta_df, total_deltas = build_delta_tables(unique_seqs, delta_k, delta_r)

    heavy = df.copy()
    heavy["_total_delta"] = heavy["Stripped.Sequence"].map(total_deltas)

    n_skipped = int((heavy["_total_delta"] == 0).sum())
    if n_skipped:
        skipped_seqs = (
            heavy.loc[heavy["_total_delta"] == 0, "Stripped.Sequence"].drop_duplicates().tolist()
        )
        print(
            f"[info] skipping heavy duplicates for {len(skipped_seqs)} peptides with no K/R "
            f"(would be identical to light): example {skipped_seqs[:3]}"
        )
        heavy = heavy[heavy["_total_delta"] > 0].copy()

    heavy["Modified.Sequence"] = (
        heavy["Modified.Sequence"].str.replace("K", "K" + mod_k, regex=False).str.replace(
            "R", "R" + mod_r, regex=False
        )
    )
    heavy["Precursor.Id"] = heavy["Modified.Sequence"] + heavy["Precursor.Charge"].astype(str)
    heavy["Precursor.Mz"] = (
        heavy["Precursor.Mz"].astype(np.float64) + heavy["_total_delta"] / heavy["Precursor.Charge"]
    ).astype(np.float32)

    heavy = heavy.merge(
        delta_df,
        on=["Stripped.Sequence", "Fragment.Type", "Fragment.Series.Number"],
        how="left",
        validate="m:1",
    )
    if heavy["frag_delta"].isna().any():
        n_na = int(heavy["frag_delta"].isna().sum())
        raise RuntimeError(f"frag_delta NaN for {n_na} rows — fragment index out of range?")
    heavy["Product.Mz"] = (
        heavy["Product.Mz"].astype(np.float64) + heavy["frag_delta"] / heavy["Fragment.Charge"]
    ).astype(np.float32)

    return heavy.drop(columns=["_total_delta", "frag_delta"])


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Duplicate a DIA-NN parquet library with SILAC heavy-labeled K/R."
    )
    ap.add_argument("in_lib", type=Path, help="Input DIA-NN spectral library (.parquet)")
    ap.add_argument(
        "-o",
        "--out-lib",
        type=Path,
        default=None,
        help="Output path. Default: <in_lib stem>.heavy-light.parquet next to input.",
    )
    ap.add_argument("--k-delta", type=float, default=8.014199, help="Mass shift for K (Da)")
    ap.add_argument("--r-delta", type=float, default=10.008269, help="Mass shift for R (Da)")
    ap.add_argument(
        "--k-mod",
        default="(UniMod:259)",
        help="String inserted after every K in Modified.Sequence",
    )
    ap.add_argument(
        "--r-mod",
        default="(UniMod:267)",
        help="String inserted after every R in Modified.Sequence",
    )
    args = ap.parse_args(argv)

    src: Path = args.in_lib
    dst: Path = args.out_lib or src.with_suffix(".heavy-light.parquet")

    print(f"[info] reading {src}")
    df = pd.read_parquet(src)
    print(f"[info] loaded {len(df):,} rows / {df['Precursor.Id'].nunique():,} precursors")

    heavy = make_heavy(
        df,
        delta_k=args.k_delta,
        delta_r=args.r_delta,
        mod_k=args.k_mod,
        mod_r=args.r_mod,
    )
    print(f"[info] heavy: {len(heavy):,} rows / {heavy['Precursor.Id'].nunique():,} precursors")

    overlap = set(df["Precursor.Id"].unique()) & set(heavy["Precursor.Id"].unique())
    if overlap:
        raise RuntimeError(f"{len(overlap)} Precursor.Id collisions between light and heavy")

    combined = pd.concat([df, heavy], ignore_index=True)[df.columns.tolist()]
    print(
        f"[info] combined: {len(combined):,} rows / {combined['Precursor.Id'].nunique():,} precursors"
    )

    print(f"[info] writing {dst}")
    combined.to_parquet(dst, index=False)
    print("[info] done")


if __name__ == "__main__":
    main()
