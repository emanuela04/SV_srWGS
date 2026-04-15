import argparse
import pandas as pd

SVTYPES = ["DEL", "DUP", "INS", "INV"]

FINAL_COLS = [
    "MANTA","DELLY","SMOOVE","MAXQV","CG_CONTENT","COVERAGE_MOSDEPTH",
    "coverage_inside","mean_insert_inside","sd_insert_inside","mean_mapq_inside",
    "n_clipped_inside","n_split_inside","n_discordant_inside",
    "coverage_left","coverage_right",
    "mean_flank_insert","DELTA_insert","mean_flank_mapq","DELTA_mapq",
    "clipped_ratio","split_ratio","discordant_ratio",
    "SVLEN_CALLER"
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare final datasets for VarTrustML"
    )
    parser.add_argument("input_file", help="Input TSV file (tab-separated).")
    parser.add_argument("-o", "--out_prefix", default="dataset",
                        help="Output prefix (default: dataset).")
    args = parser.parse_args()

    df = pd.read_csv(args.input_file, sep="\t", low_memory=False)

    df = df.replace({".": pd.NA, "": pd.NA, "NA": pd.NA, "NaN": pd.NA, "nan": pd.NA})

    if "SVLEN_CALLER" in df.columns:
        s = (
            df["SVLEN_CALLER"]
            .astype(str)
            .str.replace(r"[()\s,]", "", regex=True)
            .replace({"": pd.NA, "None": pd.NA, "nan": pd.NA})
        )
        df["SVLEN_CALLER"] = pd.to_numeric(s, errors="coerce")

    keep = [
        "VARIANT_KEY",
        "MANTA_feature", "DELLY_feature", "SMOOVE_feature",
        "MAXQV", "CG_CONTENT", "COVERAGE_MOSDEPTH",
        "coverage_inside", "mean_insert_inside", "sd_insert_inside", "mean_mapq_inside",
        "n_clipped_inside", "n_split_inside", "n_discordant_inside",
        "coverage_left", "coverage_right",
        "mean_flank_insert", "Δinsert", "mean_flank_mapq", "Δmapq",
        "clipped_ratio", "split_ratio", "discordant_ratio",
        "SVLEN_CALLER",
        "SVTYPE",
    ]

    missing = [c for c in keep if c not in df.columns]
    if missing:
        raise SystemExit(
            "Colonne mancanti nel file:\n" + "\n".join(f"- {c}" for c in missing)
        )

    df = df[keep].rename(columns={
        "MANTA_feature": "MANTA",
        "DELLY_feature": "DELLY",
        "SMOOVE_feature": "SMOOVE",
        "Δinsert": "DELTA_insert",
        "Δmapq": "DELTA_mapq",
    })

    if df["VARIANT_KEY"].isna().any():
        raise SystemExit("VARIANT_KEY with missing values")

    if not df["VARIANT_KEY"].is_unique:
        raise SystemExit("VARIANT_KEY duplicated")

    for c in FINAL_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    missing_counts = df[FINAL_COLS].isna().sum()
    print("\nMissing:")
    print(missing_counts[missing_counts > 0])

    print("\nNumber of missing:")
    print(df[FINAL_COLS].isna().any(axis=1).sum())

    df["SVTYPE"] = df["SVTYPE"].astype(str).str.upper()

    # VARIANT_KEY diventa row_id
    df = df.set_index("VARIANT_KEY")
    df.index.name = "row_id"

    for sv in SVTYPES:
        out = f"{args.out_prefix}_{sv}.csv"
        sub = df[df["SVTYPE"] == sv][FINAL_COLS]
        sub.to_csv(out, index=True, index_label="row_id")

    print("Files:")
    for sv in SVTYPES:
        print(f"  {args.out_prefix}_{sv}.csv")


if __name__ == "__main__":
    main()
