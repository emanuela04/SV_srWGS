import os
import sys
import argparse
import gzip
import pandas as pd
import pysam
import pyranges as pr
import math
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer

def log_error(message):
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def validate_file(file_path, description):
    if not os.path.isfile(file_path):
        log_error(f"{description} file not found at: {file_path}")
    if not os.access(file_path, os.R_OK):
        log_error(f"{description} file is not readable: {file_path}")


def load_repeat_masker(file_path):
    columns = ["Chromosome", "Start", "End", "repClass", "repFamily"]
    repeat_masker = pd.read_csv(file_path, sep="\t", comment="#", header=None,
                                usecols=[5, 6, 7, 11, 12], names=columns)
    repeat_masker["Start"] = repeat_masker["Start"].astype(int)
    repeat_masker["End"] = repeat_masker["End"].astype(int)
    return pr.PyRanges(repeat_masker)


def load_segmental_duplications(file_path):
    columns = ["Chromosome", "Start", "End"]
    duplication_data = pd.read_csv(file_path, sep="\t", header=None, names=columns)
    duplication_data["Start"] = duplication_data["Start"].astype(int)
    duplication_data["End"] = duplication_data["End"].astype(int)
    return pr.PyRanges(duplication_data)

def extract_representative_sample(sample_names):
    return sample_names[0]


def process_vcf_to_pyranges(vcf_path, repeats, duplications):
    with pysam.VariantFile(vcf_path) as vcf:
        representative_sample = extract_representative_sample(list(vcf.header.samples))
        vcf_data = []

        for record in vcf:
            chrom = record.chrom
            start = record.pos
            end = record.stop or record.info.get("END") or start
            svtype = record.info.get("SVTYPE", "NA")
            svlen = record.info.get("SVLEN", "NA")
            supp_vec = record.info.get("SUPP_VEC", "000")
            if len(supp_vec) < 3:
                manta = delly = smoove = '0'
            else:
                manta, delly, smoove = supp_vec[1], supp_vec[2], supp_vec[0]
             
            # adding QV from FORMAT VCF: for each tool
            all_qv = []
            for sample in record.samples:
                qv_value = record.samples[sample].get("QV")
                if qv_value is None:
                    continue
                if isinstance(qv_value, (int, float)):
                    if not math.isnan(qv_value):
                        all_qv.append(qv_value)
                elif isinstance(qv_value, (list, tuple)):
                    for v in qv_value:
                        try:
                            f = float(v)
                            if not math.isnan(f):
                                all_qv.append(f)
                        except (ValueError, TypeError):
                            continue
                elif isinstance(qv_value, str):
                    parts = qv_value.replace(",", " ").split()
                    for v in parts:
                        try:
                            f = float(v)
                            if not math.isnan(f):
                                all_qv.append(f)
                        except (ValueError, TypeError):
                            continue
            
            max_qv = "."
            if all_qv:
                max_qv = int(max(all_qv))
            

 
            
            vcf_data.append({
                "SAMPLE_NAME": representative_sample,
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "SVTYPE_CALLER": svtype,
                "SVLEN_CALLER": svlen,
                "MANTA": manta,
                "DELLY": delly,
                "SMOOVE": smoove,
                "MAXQV": max_qv
            })

    vcf_ranges = pr.PyRanges(pd.DataFrame(vcf_data))
    repeats_overlap = vcf_ranges.join(repeats)
    duplications_overlap = vcf_ranges.join(duplications)
    vcf_df = vcf_ranges.as_df()

    vcf_keys = vcf_df[['Chromosome', 'Start', 'End']].apply(tuple, axis=1)
    repeats_keys = set(repeats_overlap.as_df()[['Chromosome', 'Start', 'End']].apply(tuple, axis=1))
    duplications_keys = set(duplications_overlap.as_df()[['Chromosome', 'Start', 'End']].apply(tuple, axis=1))

    vcf_df['OVERLAPS_REPEATS'] = vcf_keys.apply(lambda x: "True" if x in repeats_keys else ".")
    vcf_df['OVERLAPS_SEG_DUP'] = vcf_keys.apply(lambda x: "True" if x in duplications_keys else ".")

    return pr.PyRanges(vcf_df)


def process_gc_content_to_pyranges(input_file):
    data = pd.read_csv(input_file, sep="\t", comment='#', header=None)
    processed_data = data.iloc[:, [0, 1, 2, 6]].copy()
    processed_data.columns = ["Chromosome", "Start", "End", "GC_Content"]
    processed_data["Start"] = processed_data["Start"].astype(int)
    processed_data["End"] = processed_data["End"].astype(int)
    return pr.PyRanges(processed_data)

def load_coverage_to_pyranges(input_file):
    with gzip.open(input_file, 'rt') as gz_file:
        df = pd.read_csv(gz_file, sep="\t", header=None, usecols=[0, 1, 2, 4],
                         names=["Chromosome", "Start", "End", "Coverage"])
        df["Start"] = df["Start"].astype(int)
        df["End"] = df["End"].astype(int)
    return {"sample": pr.PyRanges(df)}

#def load_coverage_to_pyranges(input_file):
#    sample_name = os.path.basename(input_file).replace(".regions.bed.gz", "")
#    with gzip.open(input_file, 'rt') as gz_file:
#        df = pd.read_csv(gz_file, sep="\t", header=None, usecols=[0, 1, 2, 4],
#                         names=["Chromosome", "Start", "End", "Coverage"])
#        df["Start"] = df["Start"].astype(int)
#        df["End"] = df["End"].astype(int)
#    return {sample_name: pr.PyRanges(df)}


def combine_bed_files_single_tsv(gc_pr, coverage_dict, sv_df, output_file):
    gc_dict = {(row.Chromosome, row.Start, row.End): row.GC_Content for _, row in gc_pr.df.iterrows()}
    unique_samples = sv_df["SAMPLE_NAME"].unique()

    header = [
        "SAMPLE_NAME", "CHROM_CALLER", "POS_CALLER", "END_CALLER", "SVTYPE_CALLER",
        "SVLEN_CALLER", "MANTA", "DELLY", "SMOOVE","MAXQV" , "OVERLAPS_REPEATS",
        "OVERLAPS_SEG_DUP", "CG_CONTENT", "COVERAGE_MOSDEPTH"
    ]

    with open(output_file, "w") as out_file:
        out_file.write("\t".join(header) + "\n")

        for sample_name in unique_samples:
            #coverage_pr = coverage_dict.get(sample_name)
            coverage_pr = list(coverage_dict.values())[0]  
            if coverage_pr is None:
                print(f"Warning: No coverage data found for sample {sample_name}")
                continue

            sample_sv_df = sv_df[sv_df["SAMPLE_NAME"] == sample_name].copy()
            sample_sv_df.rename(
                columns={"Chromosome": "CHROM_CALLER", "Start": "POS_CALLER", "End": "END_CALLER"},
                inplace=True, errors="ignore"
            )

            for _, row in sample_sv_df.iterrows():
                chrom = row["CHROM_CALLER"]
                start = row["POS_CALLER"] - 1
                end = row["END_CALLER"]
                gc_content = gc_dict.get((chrom, start, end), ".")

                coverage_subset = coverage_pr[
                    (coverage_pr.Chromosome == chrom) &
                    (coverage_pr.Start == start) &
                    (coverage_pr.End == end)
                ]
                coverage = coverage_subset.Coverage.values[0] if len(coverage_subset) == 1 else "."

                row_data = [
                    row["SAMPLE_NAME"], chrom, start, end, row["SVTYPE_CALLER"], row["SVLEN_CALLER"],
                    row.get("MANTA", 0), row.get("DELLY", 0), row.get("SMOOVE", 0),row.get("MAXQV", "."),
                    row.get("OVERLAPS_REPEATS", "."), row.get("OVERLAPS_SEG_DUP", "."),
                    gc_content, coverage
                ]

                out_file.write("\t".join(map(str, row_data)) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Efficiently filter structural variants using PyRanges.")
    parser.add_argument("--repeat_masker", required=True, help="Path to RepeatMasker TSV file")
    parser.add_argument("--segmental_duplications", required=True, help="Path to Segmental Duplications TSV file")
    parser.add_argument("--vcf_file", required=True, help="Path to a single VCF file")
    parser.add_argument("--cg_content", required=True, help="Path to CG content file (_merged_gc_content.txt)")
    parser.add_argument("--coverage_mosdepth", required=True, help="Path to mosdepth coverage file (.regions.bed.gz)")
    parser.add_argument("--output_path", required=True, help="Path to output TSV file")
    args = parser.parse_args()

    validate_file(args.repeat_masker, "RepeatMasker")
    validate_file(args.segmental_duplications, "Segmental Duplications")
    validate_file(args.vcf_file, "VCF")
    validate_file(args.cg_content, "GC content file")
    validate_file(args.coverage_mosdepth, "Coverage mosdepth file")

    try:
        repeats_pr = load_repeat_masker(args.repeat_masker)
        segdup_pr = load_segmental_duplications(args.segmental_duplications)
        sv_pr = process_vcf_to_pyranges(args.vcf_file, repeats_pr, segdup_pr)
        sv_df = sv_pr.as_df()
        gc_pr = process_gc_content_to_pyranges(args.cg_content)
        coverage_dict = load_coverage_to_pyranges(args.coverage_mosdepth)

        combine_bed_files_single_tsv(gc_pr, coverage_dict, sv_df, args.output_path)
        print(f"Finished processing. Merged TSV saved at: {args.output_path}")

        # --- One-Hot Encoding ---
        df = pd.read_csv(args.output_path, sep="\t")
        encoder: ColumnTransformer = ColumnTransformer(
            [("ohe", OneHotEncoder(drop="if_binary", sparse_output=False), ["SVTYPE_CALLER"])],
            remainder="passthrough",
            verbose_feature_names_out=False
        ).set_output(transform="pandas")
        df = encoder.fit_transform(df)
        # --- Subset of columns required ---
        final_columns = [
            "SAMPLE_NAME",
            "SVTYPE_CALLER_DEL",
            "SVTYPE_CALLER_DUP",
            "SVTYPE_CALLER_INS",
            "SVTYPE_CALLER_INV",
            "SVLEN_CALLER",
            "MANTA",
            "DELLY",
            "SMOOVE",
            "MAXQV",
            "CG_CONTENT",
            "COVERAGE_MOSDEPTH"
        ]
        
        for col in final_columns:
            if col not in df.columns:
                df[col] = "."
        
        df1 = df[final_columns]
        df.to_csv(args.output_path, sep="\t", index=False)
        #df1.to_csv("for_ml.csv", sep ="\t", index= False)
        base, ext = os.path.splitext(args.output_path)
        classifier_path = f"{base}_for_classifier{ext}"

        df1.to_csv(classifier_path, sep="\t", index=False)
        print(f"Final TSV with One-Hot Encoding saved at: {args.output_path}")

    except Exception as e:
        log_error(f"An unexpected error occurred during processing: {e}")

if __name__ == "__main__":
    main()
