#!/usr/bin/env python3

import os
import sys
import argparse
import gzip
import math
import pandas as pd
import pysam
pysam.set_verbosity(0)

import pyranges as pr

from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer


def log_error(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# Load inputs
def load_repeat_masker(file_path):
    # prendiamo solo chr/start/end (il resto non serve per overlap boolean)
    df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[5, 6, 7],
        names=["Chromosome", "Start", "End"]
    )
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return pr.PyRanges(df).merge()  # merge per ridurre tantissimo le regioni


def load_segmental_duplications(file_path):
    df = pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End"]
    )
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return pr.PyRanges(df).merge()


def load_gc_content(file_path):
    df = pd.read_csv(file_path, sep="\t", comment="#", header=None)
    df = df.iloc[:, [0, 1, 2, 6]].copy()
    df.columns = ["Chromosome", "Start", "End", "GC_CONTENT"]
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["GC_CONTENT"] = pd.to_numeric(df["GC_CONTENT"], errors="coerce")
    return df


def load_coverage(file_path):
    with gzip.open(file_path, "rt") as f:
        df = pd.read_csv(
            f,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 4],
            names=["Chromosome", "Start", "End", "COVERAGE_MOSDEPTH"]
        )
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["COVERAGE_MOSDEPTH"] = pd.to_numeric(df["COVERAGE_MOSDEPTH"], errors="coerce")
    return df


def load_feature_matrix(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.rename(columns={"CHROM": "Chromosome", "BPl": "Start", "BPr": "End"})
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df = df.drop_duplicates(subset=["Chromosome", "Start"])
    return df


#vcf open

def process_vcf(vcf_path):
    rows = []
    with pysam.VariantFile(vcf_path) as vcf:
        sample = list(vcf.header.samples)[0] if len(vcf.header.samples) else "NA"

        for rec in vcf:
            chrom = rec.chrom
            start = int(rec.pos)
            end = rec.stop or rec.info.get("END") or start
            end = int(end)

            svtype = rec.info.get("SVTYPE", "NA")
            svlen = rec.info.get("SVLEN", "NA")

            supp = rec.info.get("SUPP_VEC", "000")
            supp = supp if isinstance(supp, str) and len(supp) == 3 else "000"
            smoove, manta, delly = supp[0], supp[1], supp[2]

            #qvs = []
            #for s in rec.samples:
            #    qv = rec.samples[s].get("QV")
            #    if qv is None:
            #        continue
            #    try:
            #        if isinstance(qv, (list, tuple)):
            #            for x in qv:
            #                fx = float(x)
            #                if not math.isnan(fx):
            #                    qvs.append(fx)
            #        else:
            #            fx = float(qv)
            #            if not math.isnan(fx):
            #                qvs.append(fx)
            #    except Exception:
            #        pass

            #max_qv = int(max(qvs)) if qvs else "."
            gq = rec.samples[sample].get("GQ")
            max_qv = int(gq) if gq not in (None, ".") else "."
            rows.append({
                "SAMPLE_NAME": sample,
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "SVTYPE_CALLER": svtype,
                "SVLEN_CALLER": svlen,
                "MANTA_feature": int(manta) if str(manta).isdigit() else 0,
                "DELLY_feature": int(delly) if str(delly).isdigit() else 0,
                "SMOOVE_feature": int(smoove) if str(smoove).isdigit() else 0,
                "MAXQV": max_qv
            })

    return pd.DataFrame(rows)


# RAM-safe annotations
def add_overlap_flag(sv_pr, anno_pr, out_col):
    # join produce solo le righe che overlappano -> poi set boolean in modo leggero
    hits = sv_pr.join(anno_pr).as_df()[["Chromosome", "Start", "End"]].drop_duplicates()
    sv_df = sv_pr.as_df()
    sv_df[out_col] = sv_df.merge(
        hits.assign(_hit=True),
        on=["Chromosome", "Start", "End"],
        how="left"
    )["_hit"].fillna(False)
    return pr.PyRanges(sv_df)


def _mean_on_overlaps_chunk(sv_chunk_df, anno_df, value_col, max_mean_span):
    """
    """
    sv_chunk_df = sv_chunk_df.copy()
    sv_chunk_df["_span"] = (sv_chunk_df["End"] - sv_chunk_df["Start"]).abs()

    short = sv_chunk_df[sv_chunk_df["_span"] <= max_mean_span]
    long_ = sv_chunk_df[sv_chunk_df["_span"] > max_mean_span]

    out = sv_chunk_df[["Chromosome", "Start", "End"]].copy()
    out[value_col] = pd.NA

    if not short.empty:
        chroms = short["Chromosome"].unique().tolist()
        min_s = int(short["Start"].min())
        max_e = int(short["End"].max())

        a = anno_df[
            (anno_df["Chromosome"].isin(chroms)) &
            (anno_df["End"] >= min_s) &
            (anno_df["Start"] <= max_e)
        ].copy()

        if not a.empty:
            sv_pr = pr.PyRanges(short[["Chromosome", "Start", "End"]])
            a_pr = pr.PyRanges(a[["Chromosome", "Start", "End", value_col]])

            j = sv_pr.join(a_pr).as_df()
            if not j.empty:
                means = (
                    j.groupby(["Chromosome", "Start", "End"], as_index=False)[value_col]
                    .mean()
                )
                out = out.merge(means, on=["Chromosome", "Start", "End"], how="left", suffixes=("", "_m"))
                out[value_col] = out[value_col].fillna(out[f"{value_col}_m"])
                out.drop(columns=[f"{value_col}_m"], inplace=True)

    if not long_.empty:
        # per ogni SV lungo, prendiamo la finestra che contiene Start
        # (se i file sono finestre BED standard, funziona e non esplode)
        for idx, r in long_.iterrows():
            chrom = r["Chromosome"]
            pos = int(r["Start"])
            hit = anno_df[
                (anno_df["Chromosome"] == chrom) &
                (anno_df["Start"] <= pos) &
                (anno_df["End"] > pos)
            ]
            if not hit.empty:
                val = hit.iloc[0][value_col]
                out.loc[(out["Chromosome"] == chrom) & (out["Start"] == r["Start"]) & (out["End"] == r["End"]), value_col] = val

    return out.drop(columns=[], errors="ignore")


def add_numeric_annotation_ram_safe(sv_df, anno_df, value_col, out_col,
                                   chunk_size=20000, max_mean_span=5_000_000):
    results = []

    # processa per cromosoma (riduce tantissimo memoria)
    for chrom, sub in sv_df.groupby("Chromosome", sort=False):
        sub = sub.sort_values("Start").reset_index(drop=True)

        # prendi solo annotation del cromosoma
        anno_chr = anno_df[anno_df["Chromosome"] == chrom].copy()
        if anno_chr.empty:
            tmp = sub[["Chromosome", "Start", "End"]].copy()
            tmp[out_col] = pd.NA
            results.append(tmp)
            continue

        for i in range(0, len(sub), chunk_size):
            chunk = sub.iloc[i:i+chunk_size]
            tmp = _mean_on_overlaps_chunk(chunk, anno_chr, value_col, max_mean_span)
            tmp = tmp.rename(columns={value_col: out_col})
            results.append(tmp)

    ann = pd.concat(results, ignore_index=True)
    merged = sv_df.merge(ann, on=["Chromosome", "Start", "End"], how="left")
    return merged


# Main
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--repeat_masker", required=True)
    parser.add_argument("--segmental_duplications", required=True)
    parser.add_argument("--vcf_file", required=True)
    parser.add_argument("--gc_content", required=True)
    parser.add_argument("--coverage_mosdepth", required=True)
    parser.add_argument("--sv_feature_matrix", required=True)
    parser.add_argument("--output_path", required=True)

    # RAM-safe knobs (opzionali)
    parser.add_argument("--chunk_size", type=int, default=20000)
    parser.add_argument("--max_mean_span", type=int, default=5_000_000,
                        help="SV più lunghi di questo usano lookup su Start (RAM-safe)")

    args = parser.parse_args()

    # 1) VCF
    sv_df = process_vcf(args.vcf_file)
    if sv_df.empty:
        log_error("VCF produced 0 variants.")

    # 2) Overlap repeats 
    sv_pr = pr.PyRanges(sv_df[["Chromosome", "Start", "End"]])
    repeats_pr = load_repeat_masker(args.repeat_masker)
    segdup_pr = load_segmental_duplications(args.segmental_duplications)

    sv_pr2 = add_overlap_flag(sv_pr, repeats_pr, "OVERLAPS_REPEATS")
    sv_pr3 = add_overlap_flag(sv_pr2, segdup_pr, "OVERLAPS_SEG_DUP")
    sv_df = sv_pr3.as_df().merge(sv_df, on=["Chromosome", "Start", "End"], how="left")

    # 3) GC + coverage (RAM-safe chunking)
    gc_df = load_gc_content(args.gc_content)
    cov_df = load_coverage(args.coverage_mosdepth)

    sv_df = add_numeric_annotation_ram_safe(
        sv_df, gc_df, value_col="GC_CONTENT", out_col="CG_CONTENT",
        chunk_size=args.chunk_size, max_mean_span=args.max_mean_span
    )

    sv_df = add_numeric_annotation_ram_safe(
        sv_df, cov_df, value_col="COVERAGE_MOSDEPTH", out_col="COVERAGE_MOSDEPTH",
        chunk_size=args.chunk_size, max_mean_span=args.max_mean_span
    )

    encoder = ColumnTransformer(
        [("ohe", OneHotEncoder(sparse_output=False), ["SVTYPE_CALLER"])],
        remainder="passthrough",
        verbose_feature_names_out=False
    ).set_output(transform="pandas")

    encoded_df = encoder.fit_transform(sv_df)

    # 5) Merge feature matrix (semplice, Chrom+Start)
    feat_df = load_feature_matrix(args.sv_feature_matrix)

    final_df = encoded_df.merge(
        feat_df,
        on=["Chromosome", "Start"],
        how="left"
    )

    # 6) Rename + key
    final_df = final_df.rename(columns={
        "Chromosome": "CHROM_CALLER",
        "Start": "POS_CALLER",
        "End": "END_CALLER"
    })

    final_df["merge_key_chrom_start_feature"] = (
        final_df["CHROM_CALLER"].astype(str) + "_" + final_df["POS_CALLER"].astype(str)
    )
    # prima riga per chiave
    before = len(final_df)
    final_df = final_df.drop_duplicates(
        subset=["merge_key_chrom_start_feature"],
        keep="first"
    )
    after = len(final_df)
    
    # 7) Output: TSV + CSV
    base = os.path.splitext(args.output_path)[0]
    final_df.to_csv(f"{base}.tsv", sep="\t", index=False)
    final_df.to_csv(f"{base}.csv", sep=",", index=False)


if __name__ == "__main__":
    main()

