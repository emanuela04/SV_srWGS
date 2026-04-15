
import argparse
import json
from collections import Counter
from typing import Optional, Tuple, List, Any

import pysam


FINAL_FORMAT_KEYS = ["GT", "GQ", "FT", "PE", "SR", "AD", "VAF"]

CORE_INFO = [
    "SVTYPE",
    "SVLEN",
    "SUPP",
    "SUPP_VEC",
    "IDLIST",
    "SVMETHOD",
]

PLUS_INFO = [
    "STRANDS",
    "STARTVARIANCE",
    "ENDVARIANCE",
    "AVG_LEN",
    "AVG_START",
    "AVG_END",
    "CIPOS",
]

OPTIONAL_INFO = [
    "CIGAR",
    "HOMLEN",
    "HOMSEQ",
]


def parse_args():
    p = argparse.ArgumentParser(description="Normalize Jasmine SV VCF")
    p.add_argument("-i", "--input", required=True, help="Input Jasmine VCF/BCF")
    p.add_argument("-o", "--output", required=True, help="Output normalized VCF/BCF")
    p.add_argument(
        "--symbolic-ins",
        action="store_true",
        help="Force INS records to REF=N and ALT=<INS>",
    )
    p.add_argument(
        "--report-json",
        default="normalization_report.json",
        help="JSON report with before/after summary",
    )
    p.add_argument(
        "--changes-tsv",
        default="normalization_changes.tsv",
        help="TSV listing record-level changes",
    )
    return p.parse_args()


def add_info_line_if_missing(header, id_, number, type_, desc):
    if id_ not in header.info:
        header.add_line(
            f'##INFO=<ID={id_},Number={number},Type={type_},Description="{desc}">'
        )


def add_format_line_if_missing(header, id_, number, type_, desc):
    if id_ not in header.formats:
        header.add_line(
            f'##FORMAT=<ID={id_},Number={number},Type={type_},Description="{desc}">'
        )


def add_filter_line_if_missing(header, id_, desc):
    if id_ not in header.filters:
        header.add_line(f'##FILTER=<ID={id_},Description="{desc}">')


def build_output_header(in_header: pysam.VariantHeader) -> pysam.VariantHeader:
    out = pysam.VariantHeader()

    # copy generic metadata where possible
    for rec in in_header.records:
        if rec.key in ("INFO", "FORMAT", "FILTER", "contig"):
            continue
        try:
            out.add_line(str(rec).strip())
        except Exception:
            pass

    # contigs
    for ctg in in_header.contigs:
        c = in_header.contigs[ctg]
        if c.length is not None:
            out.contigs.add(ctg, length=c.length)
        else:
            out.contigs.add(ctg)

    # copy filters
    for flt in in_header.filters:
        f = in_header.filters[flt]
        desc = f.description if f.description is not None else ""
        try:
            add_filter_line_if_missing(out, flt, desc)
        except Exception:
            pass

    add_filter_line_if_missing(out, "PASS", "All filters passed")

    # add canonical INFO lines
    add_info_line_if_missing(out, "END", "1", "Integer",
                             "End position of the variant described in this record")
    add_info_line_if_missing(out, "SVTYPE", "1", "String",
                             "Type of structural variant")
    add_info_line_if_missing(out, "SVLEN", ".", "Integer",
                             "Difference in length between REF and ALT alleles")
    add_info_line_if_missing(out, "SUPP", "1", "Integer",
                             "Number of supporting input calls")
    add_info_line_if_missing(out, "SUPP_VEC", "1", "String",
                             "Support vector across input files")
    add_info_line_if_missing(out, "IDLIST", ".", "String",
                             "IDs of original merged calls")
    add_info_line_if_missing(out, "SVMETHOD", "1", "String",
                             "Method used to produce this SV record")

    add_info_line_if_missing(out, "STRANDS", "1", "String",
                             "Strand orientation")
    add_info_line_if_missing(out, "STARTVARIANCE", "1", "Float",
                             "Variance of start positions in merged cluster")
    add_info_line_if_missing(out, "ENDVARIANCE", "1", "Float",
                             "Variance of end positions in merged cluster")
    add_info_line_if_missing(out, "AVG_LEN", "1", "Float",
                             "Average SV length across merged calls")
    add_info_line_if_missing(out, "AVG_START", "1", "Float",
                             "Average start position across merged calls")
    add_info_line_if_missing(out, "AVG_END", "1", "Float",
                             "Average end position across merged calls")
    add_info_line_if_missing(out, "CIPOS", "2", "Integer",
                             "Confidence interval around POS")
    add_info_line_if_missing(out, "CIGAR", "A", "String",
                             "CIGAR alignment for alternate allele")
    add_info_line_if_missing(out, "HOMLEN", ".", "Integer",
                             "Length of microhomology")
    add_info_line_if_missing(out, "HOMSEQ", ".", "String",
                             "Sequence of microhomology")

    # bookkeeping info
    add_info_line_if_missing(out, "OLD_END", "1", "Integer",
                             "Original END before normalization")
    add_info_line_if_missing(out, "OLD_REF", "1", "String",
                             "Original REF before normalization")
    add_info_line_if_missing(out, "OLD_ALT", ".", "String",
                             "Original ALT before normalization")
    add_info_line_if_missing(out, "NORM_TAG", ".", "String",
                             "Normalization operations applied")

    # final FORMAT only
    add_format_line_if_missing(out, "GT", "1", "String", "Genotype")
    add_format_line_if_missing(out, "GQ", "1", "Integer", "Genotype quality")
    add_format_line_if_missing(out, "FT", "1", "String", "Sample filter")
    add_format_line_if_missing(out, "PE", "1", "Integer", "Paired-end alternate support")
    add_format_line_if_missing(out, "SR", "1", "Integer", "Split-read alternate support")
    add_format_line_if_missing(out, "AD", "R", "Integer", "Allelic depths for REF and ALT")
    add_format_line_if_missing(out, "VAF", "1", "Float", "Variant allele fraction")

    for s in in_header.samples:
        out.add_sample(s)

    return out


def infer_svtype(rec) -> str:
    if "SVTYPE" in rec.info:
        try:
            return str(rec.info["SVTYPE"])
        except Exception:
            pass

    if rec.alts and len(rec.alts) > 0:
        alt = rec.alts[0]
        if alt.startswith("<") and alt.endswith(">"):
            return alt[1:-1]

    return "UNK"


def get_scalar(v: Any) -> Any:
    if isinstance(v, (tuple, list)):
        if len(v) == 0:
            return None
        return v[0]
    return v


def get_int(v: Any) -> Optional[int]:
    v = get_scalar(v)
    if v is None:
        return None
    try:
        return int(v)
    except Exception:
        try:
            return int(float(v))
        except Exception:
            return None


def get_float(v: Any) -> Optional[float]:
    v = get_scalar(v)
    if v is None:
        return None
    try:
        return float(v)
    except Exception:
        return None


def safe_sample_get(sample, key):
    try:
        return sample.get(key, None)
    except Exception:
        return None


def normalize_gt(gt):
    if gt is None:
        return (None, None)
    if isinstance(gt, tuple):
        if len(gt) == 1:
            return (gt[0], None)
        return tuple(gt[:2])
    return (None, None)


def infer_ins_len(rec) -> Optional[int]:
    # Prefer INFO/SVLEN
    if "SVLEN" in rec.info:
        svlen = rec.info["SVLEN"]
        if isinstance(svlen, (tuple, list)):
            if len(svlen) > 0:
                iv = get_int(svlen[0])
                if iv is not None:
                    return abs(iv)
        else:
            iv = get_int(svlen)
            if iv is not None:
                return abs(iv)

    # Explicit ALT sequence
    if rec.alts and len(rec.alts) > 0:
        alt = rec.alts[0]
        if not (alt.startswith("<") and alt.endswith(">")):
            diff = len(alt) - len(rec.ref)
            if diff > 0:
                return diff

    # Sequence in INFO
    for k in ("SVINSSEQ", "CONSENSUS"):
        if k in rec.info:
            val = rec.info[k]
            val = get_scalar(val)
            if val is not None:
                seq = str(val)
                if seq and seq != ".":
                    return len(seq)

    return None


def map_sample_fields(sample) -> dict:
    out = {
        "GT": (None, None),
        "GQ": None,
        "FT": None,
        "PE": None,
        "SR": None,
        "AD": None,
        "VAF": None,
    }

    # Common
    gt = normalize_gt(safe_sample_get(sample, "GT"))
    if gt != (None, None):
        out["GT"] = gt

    gq = get_int(safe_sample_get(sample, "GQ"))
    if gq is not None:
        out["GQ"] = gq

    ft = safe_sample_get(sample, "FT")
    if ft is not None:
        if isinstance(ft, (tuple, list)):
            ft = ";".join(str(x) for x in ft if x is not None)
        out["FT"] = str(ft)

    # Manta-like
    pr = safe_sample_get(sample, "PR")
    sr_manta = safe_sample_get(sample, "SR")

    if isinstance(pr, (tuple, list)) and len(pr) >= 2:
        pr_ref = get_int(pr[0]) or 0
        pr_alt = get_int(pr[1]) or 0
        out["PE"] = pr_alt
    else:
        pr_ref = pr_alt = None

    if isinstance(sr_manta, (tuple, list)) and len(sr_manta) >= 2:
        sr_ref = get_int(sr_manta[0]) or 0
        sr_alt = get_int(sr_manta[1]) or 0
        out["SR"] = sr_alt
    else:
        sr_ref = sr_alt = None

    if pr_ref is not None or sr_ref is not None:
        ref_total = (pr_ref or 0) + (sr_ref or 0)
        alt_total = (pr_alt or 0) + (sr_alt or 0)
        out["AD"] = (ref_total, alt_total)
        denom = ref_total + alt_total
        if denom > 0:
            out["VAF"] = round(alt_total / denom, 6)
        return out

    # Delly-like
    dr = get_int(safe_sample_get(sample, "DR"))
    dv = get_int(safe_sample_get(sample, "DV"))
    rr = get_int(safe_sample_get(sample, "RR"))
    rv = get_int(safe_sample_get(sample, "RV"))

    if any(x is not None for x in (dr, dv, rr, rv)):
        ref_total = (dr or 0) + (rr or 0)
        alt_total = (dv or 0) + (rv or 0)
        out["PE"] = dv if dv is not None else 0
        out["SR"] = rv if rv is not None else 0
        out["AD"] = (ref_total, alt_total)
        denom = ref_total + alt_total
        if denom > 0:
            out["VAF"] = round(alt_total / denom, 6)
        return out

    # Smoove-like
    ro = safe_sample_get(sample, "RO")
    ao = safe_sample_get(sample, "AO")
    ap = safe_sample_get(sample, "AP")
    ass = safe_sample_get(sample, "AS")
    ab = safe_sample_get(sample, "AB")

    ro_i = get_int(ro)
    ao_i = get_int(ao)
    ap_i = get_int(ap)
    as_i = get_int(ass)
    ab_f = get_float(ab)

    if any(x is not None for x in (ro_i, ao_i, ap_i, as_i, ab_f)):
        out["PE"] = ap_i if ap_i is not None else 0
        out["SR"] = as_i if as_i is not None else 0

        if ro_i is not None or ao_i is not None:
            ref_total = ro_i or 0
            alt_total = ao_i or 0
            out["AD"] = (ref_total, alt_total)
            denom = ref_total + alt_total
            if denom > 0:
                out["VAF"] = round(alt_total / denom, 6)
        elif ab_f is not None:
            out["VAF"] = round(ab_f, 6)

        return out

    return out


def summarize_vcf(path: str) -> dict:
    vf = pysam.VariantFile(path)

    summary = {
        "records_total": 0,
        "records_by_svtype": Counter(),
        "format_keys_global": Counter(),
        "info_keys_global": Counter(),
        "sample_count_distribution": Counter(),
        "records_with_ad": 0,
        "records_with_vaf": 0,
        "ins_records": 0,
        "ins_point_like": 0,
        "ins_interval_like": 0,
    }

    for rec in vf:
        summary["records_total"] += 1
        svtype = infer_svtype(rec)
        summary["records_by_svtype"][svtype] += 1

        for k in rec.info.keys():
            summary["info_keys_global"][k] += 1

        sample_names = list(rec.samples.keys())
        summary["sample_count_distribution"][len(sample_names)] += 1

        seen_fmt = set()
        for s in sample_names:
            smp = rec.samples[s]
            for k, v in smp.items():
                if v is not None:
                    seen_fmt.add(k)
            if "AD" in smp and smp["AD"] is not None:
                summary["records_with_ad"] += 1
            if "VAF" in smp and smp["VAF"] is not None:
                summary["records_with_vaf"] += 1

        for k in seen_fmt:
            summary["format_keys_global"][k] += 1

        if svtype == "INS":
            summary["ins_records"] += 1
            end = rec.stop if rec.stop is not None else rec.pos
            svlen = None
            if "SVLEN" in rec.info:
                svlen = rec.info["SVLEN"]
                if isinstance(svlen, (tuple, list)) and len(svlen) > 0:
                    svlen = get_int(svlen[0])
                else:
                    svlen = get_int(svlen)

            if end == rec.pos:
                summary["ins_point_like"] += 1
            if svlen is not None and end == rec.pos + abs(svlen):
                summary["ins_interval_like"] += 1

    vf.close()

    for key in (
        "records_by_svtype",
        "format_keys_global",
        "info_keys_global",
        "sample_count_distribution",
    ):
        summary[key] = dict(summary[key])

    return summary


def copy_selected_info(src_rec, dst_rec):
    for key in CORE_INFO + PLUS_INFO + OPTIONAL_INFO:
        if key == "END":
            continue
        if key in src_rec.info:
            try:
                dst_rec.info[key] = src_rec.info[key]
            except Exception:
                pass


def main():
    args = parse_args()

    before = summarize_vcf(args.input)

    invcf = pysam.VariantFile(args.input)
    out_header = build_output_header(invcf.header)
    outvcf = pysam.VariantFile(args.output, "w", header=out_header)

    change_stats = {
        "records_written": 0,
        "ins_seen": 0,
        "ins_end_changed": 0,
        "ins_refalt_rewritten": 0,
        "records_with_ad": 0,
        "records_with_vaf": 0,
    }

    with open(args.changes_tsv, "w") as tsv:
        tsv.write(
            "\t".join(
                [
                    "CHROM",
                    "POS",
                    "ID",
                    "SVTYPE",
                    "OLD_END",
                    "NEW_END",
                    "OLD_REF",
                    "NEW_REF",
                    "OLD_ALT",
                    "NEW_ALT",
                    "NORM_TAG",
                ]
            )
            + "\n"
        )

        for rec in invcf:
            svtype = infer_svtype(rec)

            old_end = rec.stop if rec.stop is not None else rec.pos
            old_ref = rec.ref if rec.ref else "N"
            old_alt = ",".join(rec.alts) if rec.alts else "."

            new_ref = old_ref
            new_alts = tuple(rec.alts) if rec.alts else (f"<{svtype}>",)
            new_end = old_end
            norm_tags: List[str] = []

            if svtype == "INS":
                change_stats["ins_seen"] += 1
                ins_len = infer_ins_len(rec)
                if ins_len is not None:
                    proposed_end = rec.pos + abs(ins_len)
                    if proposed_end != old_end:
                        new_end = proposed_end
                        norm_tags.append("INS_END_FROM_SVLEN")
                        change_stats["ins_end_changed"] += 1

                if args.symbolic_ins:
                    if new_ref != "N" or new_alts != ("<INS>",):
                        new_ref = "N"
                        new_alts = ("<INS>",)
                        norm_tags.append("INS_SYMBOLIC_ALT")
                        change_stats["ins_refalt_rewritten"] += 1

            # Create normalized record
            out_rec = outvcf.new_record(
                contig=rec.contig,
                start=rec.start,
                stop=new_end,
                id=rec.id,
                qual=rec.qual,
                alleles=(new_ref,) + tuple(new_alts),
            )

            # FILTER
            filters = list(rec.filter.keys())
            if filters:
                for f in filters:
                    out_rec.filter.add(f)
            else:
                out_rec.filter.add("PASS")

            # INFO restricted + canonical
            copy_selected_info(rec, out_rec)
            out_rec.info["SVTYPE"] = svtype
            if "SVLEN" in rec.info:
                try:
                    out_rec.info["SVLEN"] = rec.info["SVLEN"]
                except Exception:
                    pass
            out_rec.info["OLD_END"] = old_end
            out_rec.info["OLD_REF"] = old_ref
            if rec.alts:
                out_rec.info["OLD_ALT"] = tuple(rec.alts)
            if norm_tags:
                out_rec.info["NORM_TAG"] = tuple(norm_tags)

            # Samples
            for s in invcf.header.samples:
                src = rec.samples[s]
                dst = out_rec.samples[s]
                mapped = map_sample_fields(src)

                dst["GT"] = mapped["GT"]

                if mapped["GQ"] is not None:
                    dst["GQ"] = mapped["GQ"]

                dst["FT"] = mapped["FT"] if mapped["FT"] is not None else "PASS"

                if mapped["PE"] is not None:
                    dst["PE"] = mapped["PE"]

                if mapped["SR"] is not None:
                    dst["SR"] = mapped["SR"]

                if mapped["AD"] is not None:
                    dst["AD"] = mapped["AD"]
                    change_stats["records_with_ad"] += 1

                if mapped["VAF"] is not None:
                    dst["VAF"] = mapped["VAF"]
                    change_stats["records_with_vaf"] += 1

            outvcf.write(out_rec)
            change_stats["records_written"] += 1

            new_alt_str = ",".join(new_alts) if new_alts else "."
            tsv.write(
                "\t".join(
                    [
                        str(rec.contig),
                        str(rec.pos),
                        rec.id if rec.id else ".",
                        svtype,
                        str(old_end),
                        str(new_end),
                        old_ref,
                        new_ref,
                        old_alt,
                        new_alt_str,
                        ",".join(norm_tags) if norm_tags else ".",
                    ]
                )
                + "\n"
            )

    invcf.close()
    outvcf.close()

    after = summarize_vcf(args.output)

    report = {
        "before": before,
        "after": after,
        "changes": change_stats,
    }

    with open(args.report_json, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
