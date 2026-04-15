import pysam
import pandas as pd
import argparse
import os

SVTYPES = ["del", "dup", "ins", "inv"]


def load_predictions(pred_dir, sample_name):
    preds = {}
    thresholds = {}

    for sv in SVTYPES:
        f = os.path.join(pred_dir, f"{sample_name}_prediction_{sv}.csv")

        if not os.path.exists(f):
            raise FileNotFoundError(f"Prediction file non trovato: {f}")

        df = pd.read_csv(f)

        required = {
            "prob_class_0",
            "prob_class_1",
            "predicted_class",
            "optimal_threshold",
        }
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{f}: colonne mancanti: {sorted(missing)}")

        thresholds[sv] = float(df["optimal_threshold"].iloc[0])
        preds[sv] = df

    return preds, thresholds


def main():
    parser = argparse.ArgumentParser(
        description="Annotate VCF with classifier predictions using record order per SVTYPE."
    )
    parser.add_argument("-i", "--input_vcf", required=True, help="Input VCF")
    parser.add_argument("-o", "--output_vcf", required=True, help="Output VCF")
    parser.add_argument(
        "-p", "--pred_dir", required=True,
        help="Directory containing HG002_prediction_del.csv, HG002_prediction_dup.csv, HG002_prediction_ins.csv, HG002_prediction_inv.csv"
    )
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.input_vcf)

    samples = list(vcf.header.samples)
    if len(samples) != 1:
        raise ValueError(f"multisample_vcf: {samples}")

    sample_name = samples[0]

    preds, thresholds = load_predictions(args.pred_dir, sample_name)
    counters = {sv: 0 for sv in SVTYPES}
    annotated = {sv: 0 for sv in SVTYPES}

    if "STRANDS" not in vcf.header.info:
        vcf.header.add_line(
            '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Strand orientation">'
        )

    for sv, thr in thresholds.items():
        vcf.header.add_line(f"##vartrust_{sv}_threshold={thr}")

    if "PC0" not in vcf.header.formats:
        vcf.header.formats.add("PC0", 1, "Float", "Probability class 0 (False positive)")
    if "PC1" not in vcf.header.formats:
        vcf.header.formats.add("PC1", 1, "Float", "Probability class 1 (True positive)")
    if "PRED" not in vcf.header.formats:
        vcf.header.formats.add("PRED", 1, "Integer", "Predicted class")

    out = pysam.VariantFile(args.output_vcf, "w", header=vcf.header)

    for rec in vcf:
        svtype = rec.info.get("SVTYPE")

        if svtype is None:
            out.write(rec)
            continue

        svtype = str(svtype).lower()

        if svtype not in SVTYPES:
            out.write(rec)
            continue

        idx = counters[svtype]
        pred_df = preds[svtype]

        if idx >= len(pred_df):
            out.write(rec)
            continue

        row = pred_df.iloc[idx]

        pc0 = float(row["prob_class_0"])
        pc1 = float(row["prob_class_1"])
        pred = int(row["predicted_class"])

        for sample in rec.samples:
            rec.samples[sample]["PC0"] = pc0
            rec.samples[sample]["PC1"] = pc1
            rec.samples[sample]["PRED"] = pred

        counters[svtype] += 1
        annotated[svtype] += 1

        out.write(rec)

    out.close()

    print("\nAnnotated variants:")
    for sv in SVTYPES:
        print(f"{sv}: {annotated[sv]}")


if __name__ == "__main__":
    main()
