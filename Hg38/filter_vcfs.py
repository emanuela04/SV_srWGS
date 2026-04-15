#!/usr/bin/env python3
import pysam
import argparse
import os

# ------------------------------
# Parser degli argomenti
# ------------------------------
parser = argparse.ArgumentParser(
    description="Corregge END delle INS, calcola SVLEN se mancante e applica filtro SVLEN/PASS"
)
parser.add_argument("-i", "--input", required=True, help="File VCF o BCF di input")
parser.add_argument("-o", "--output", required=True, help="File VCF o BCF di output")
parser.add_argument("-s", "--skipped", required=True, help="TSV delle varianti scartate")
args = parser.parse_args()

# ------------------------------
# Determina modalità scrittura
# ------------------------------
ext = os.path.splitext(args.output)[1].lower()
mode = "wb" if ext == ".bcf" else "w"

# ------------------------------
# Apri VCF/BCF
# ------------------------------
vcf_in = pysam.VariantFile(args.input)
vcf_out = pysam.VariantFile(args.output, mode, header=vcf_in.header)

# Controlla se il VCF è Smoove
is_smoove = "smoove" in os.path.basename(args.input).lower()

# Apri file TSV per varianti scartate
skipped_file = open(args.skipped, "w")
skipped_file.write("CHROM\tPOS\tEND\tSVTYPE\tSVLEN\tFILTER\n")

# Filtro SVLEN
min_svlen = 50
max_neg_svlen = -50

for rec in vcf_in:
    svtype = rec.info.get("SVTYPE")
    svlen = rec.info.get("SVLEN")
    
    # ------------------------------
    # 1) Calcolo SVLEN se mancante (tutte le SVTYPE)
    # ------------------------------
    if svlen is None and rec.stop is not None:
        svlen = rec.stop - rec.pos
        rec.info["SVLEN"] = svlen
    elif isinstance(svlen, tuple):
        svlen = svlen[0]

    # ------------------------------
    # 2) Correzione END per INS
    # ------------------------------
    if svtype == "INS" and svlen is not None:
        rec.stop = rec.pos + int(svlen)

    # ------------------------------
    # 3) Controllo filtro PASS
    # ------------------------------
    filter_status = ",".join(rec.filter.keys()) if rec.filter.keys() else "PASS"
    write_record = is_smoove or filter_status == "PASS"

    # ------------------------------
    # 4) Filtro SVLEN
    # ------------------------------
    if write_record:
        svlen_val = svlen if svlen is not None else "."
        if svlen_val != "." and not (svlen_val >= min_svlen or svlen_val <= max_neg_svlen):
            # fuori dai limiti -> scartato
            write_record = False

    # ------------------------------
    # 5) Scrittura
    # ------------------------------
    if write_record:
        vcf_out.write(rec)
    else:
        end_val = rec.stop if rec.stop else "."
        svlen_val = svlen if svlen is not None else "."
        skipped_file.write(f"{rec.chrom}\t{rec.pos}\t{end_val}\t{svtype}\t{svlen_val}\t{filter_status}\n")

vcf_in.close()
vcf_out.close()
skipped_file.close()

print(f"File corretto salvato in {args.output}")
print(f"Varianti scartate salvate in {args.skipped}")

