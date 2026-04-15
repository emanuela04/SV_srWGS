#!/usr/bin/env python3
import pysam
import argparse
import os

# ------------------------------
# Parser degli argomenti
# ------------------------------
parser = argparse.ArgumentParser(
    description="Corregge l'END delle INS in un VCF/BCF usando POS + SVLEN"
)
parser.add_argument(
    "-i", "--input", required=True, help="File VCF o BCF di input"
)
parser.add_argument(
    "-o", "--output", required=True, help="File VCF o BCF di output"
)
args = parser.parse_args()

# ------------------------------
# Determina il formato di output
# ------------------------------
ext = os.path.splitext(args.output)[1].lower()
mode = "wb" if ext == ".bcf" else "w"  # write BCF o VCF

# ------------------------------
# Apri VCF/BCF
# ------------------------------
vcf_in = pysam.VariantFile(args.input)
vcf_out = pysam.VariantFile(args.output, mode, header=vcf_in.header)

# ------------------------------
# Aggiorna END delle INS
# ------------------------------
for rec in vcf_in:
    if rec.info.get("SVTYPE") == "INS":
        svlen = rec.info.get("SVLEN")
        if svlen is not None:  
            if isinstance(svlen, tuple):
                svlen = svlen[0]  # prendi il primo valore se è multiplo
            rec.stop = rec.pos + int(svlen)  # aggiorna END
    vcf_out.write(rec)

vcf_in.close()
vcf_out.close()

print(f"File corretto salvato in {args.output}")

