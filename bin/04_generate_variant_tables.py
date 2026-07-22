#!/usr/bin/env python3
"""
Generate variant tables from the canonical merged multi-sample VCF.

The genotypes, variant set, TYPE and snpEff annotation are all read from a
single biallelic, snpEff-annotated VCF derived from
${POP}.merged.vcf.gz (built by 04_run_generate_variant_vcf.sh). That VCF is the
shared source of truth for both these tables and the TreeTime VCF, so the two
artifacts agree on every shared (CHR,POS,REF,ALT) key by construction, and the
self-alignment repeat exclusion applied once when building it carries through
here automatically (excluded sites are simply absent from the VCF).

All per-sample filtering (QUAL/DP/AO-frac/RO-frac thresholds, het masking) and
the mosdepth low-coverage mask were already baked into the GT field upstream, so
this script does no genotype calling of its own:
    GT 1/1 (or 1)   -> "1"   confident ALT call
    GT 0/0 (or 0)   -> "0"   covered reference
    GT ./. / .      -> "-"   missing / masked

Usage:
    python 04_generate_variant_tables.py --pop <POP> --merged-vcf <PATH> --gff <PATH>
        [--out-dir PATH]

Outputs:
  var.annot.tab         all variants with at least one ALT call
  disruptive.annot.tab  disruptive variants only
  snp.annot.tab         SNPs/MNPs with at least one ALT call
  snp.gt.tab            SNP genotype matrix
  indel.annot.tab       INDELs with at least one ALT call
  indel.gt.tab          INDEL genotype matrix
"""

import sys
import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

from pipeline_helpers import classify_isolates, is_missing_gt

ANNOT_COLS = [
    "FTYPE", "STRAND", "NT_POS", "AA_POS", "EFFECT", "IMPACT",
    "NT_CHANGE", "AA_CHANGE", "ID", "LOCUS_TAG", "GENE", "PRODUCT",
]

TYPE_PRIORITY = {"snp": 0, "mnp": 1, "ins": 2, "del": 2, "complex": 3}

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def flip_nt_change_to_ref_strand(nt_change):
    """Rewrite snpEff HGVS.c-style NT_CHANGE (gene strand) onto the reference
    strand: positions are preserved, base/sequence tokens are reverse-complemented
    so the bases match the REF/ALT columns. Returns the input unchanged if the
    notation doesn't match a known pattern."""
    if not nt_change:
        return nt_change
    # Substitution: "124G>A"
    m = re.match(r'^(\d+)([ACGTacgt]+)>([ACGTacgt]+)$', nt_change)
    if m:
        pos, ref, alt = m.groups()
        return f"{pos}{_revcomp(ref)}>{_revcomp(alt)}"
    # Deletion: "161delC", "578_579delCC", "161del"
    m = re.match(r'^(\d+(?:_\d+)?)del([ACGTacgt]*)$', nt_change)
    if m:
        pos, bases = m.groups()
        return f"{pos}del{_revcomp(bases)}"
    # Insertion: "286703_286704insGGC"
    m = re.match(r'^(\d+_\d+)ins([ACGTacgt]+)$', nt_change)
    if m:
        pos, bases = m.groups()
        return f"{pos}ins{_revcomp(bases)}"
    # Duplication: "123dupA", "123_125dupABC", "123dup"
    m = re.match(r'^(\d+(?:_\d+)?)dup([ACGTacgt]*)$', nt_change)
    if m:
        pos, bases = m.groups()
        return f"{pos}dup{_revcomp(bases)}"
    # Delins: "123_124delinsAB", "123delinsAB"
    m = re.match(r'^(\d+(?:_\d+)?)delins([ACGTacgt]+)$', nt_change)
    if m:
        pos, bases = m.groups()
        return f"{pos}delins{_revcomp(bases)}"
    return nt_change


# ── Isolate classification ─────────────────────────────────────────────────────
# classify_isolates lives in pipeline_helpers so the canonical-VCF pass can reuse
# the identical PR/PO split for its per-site missingness QC.


# ── GFF lookup ────────────────────────────────────────────────────────────────

def load_gff_lookup(gff_path):
    """
    Parse ref.gff and return dict keyed by ID= (snpEff Gene_ID).
    Each entry: {"id": ID, "locus_tag": locus_tag, "bakta_id": bakta_id,
                 "gene": gene, "product": product, "strand": strand}.
    URL-decodes %2C -> , in product strings.
    """
    lookup = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            strand = parts[6]
            attrs = parts[8]
            feat_id = ""
            locus_tag = ""
            bakta_id = ""
            gene = ""
            product = ""
            for attr in attrs.split(";"):
                if attr.startswith("ID="):
                    feat_id = attr[3:]
                elif attr.startswith("locus_tag="):
                    locus_tag = attr[10:]
                elif attr.startswith("bakta_id="):
                    bakta_id = attr[9:]
                elif attr.startswith("gene="):
                    gene = attr[5:]
                elif attr.startswith("product="):
                    product = attr[8:].replace("%2C", ",")
            if feat_id:
                lookup[feat_id] = {"id": feat_id, "locus_tag": locus_tag, "bakta_id": bakta_id,
                                   "gene": gene, "product": product, "strand": strand}
    return lookup


# ── VCF parsing ───────────────────────────────────────────────────────────────

def parse_ann_field(info, gff_lookup, target_alt=None):
    """Parse snpEff ANN= field from INFO. Coding annotations preferred over intergenic.

    When target_alt is given, only consider ANN entries whose allele field (field 0)
    matches the target ALT allele. This correctly handles multi-allelic VCF lines.
    """
    empty = {col: "" for col in ANNOT_COLS}

    ann_str = ""
    for kv in info.split(";"):
        if kv.startswith("ANN="):
            ann_str = kv[4:]
            break

    if not ann_str:
        return empty

    best = None
    for ann_entry in ann_str.split(","):
        parts = ann_entry.split("|")
        if len(parts) < 10:
            continue
        if target_alt is not None and parts[0] != target_alt:
            continue
        effect = parts[1]
        if best is None:
            best = parts
        elif "intergenic" in best[1] and "intergenic" not in effect:
            best = parts

    if best is None:
        return empty

    effect = best[1]
    impact = best[2] if len(best) > 2 else ""
    # best[3] = Gene_Name (matches ID= in GFF, used as lookup key)
    gene_id_ann = best[3] if len(best) > 3 else ""
    feature_type = best[5] if len(best) > 5 else ""
    ftype = "CDS" if feature_type not in ("", "intergenic_region") else ""
    nt_pos = best[12] if len(best) > 12 else ""
    aa_pos = best[13] if len(best) > 13 else ""
    hgvs_c = best[9] if len(best) > 9 else ""
    nt_change = re.sub(r'^[cn]\.', '', hgvs_c)
    hgvs_p = best[10] if len(best) > 10 else ""
    aa_change = re.sub(r'^p\.', '', hgvs_p)

    gff_info = gff_lookup.get(gene_id_ann, {})
    strand = gff_info.get("strand", "")
    if strand == "-":
        nt_change = flip_nt_change_to_ref_strand(nt_change)

    return {
        "FTYPE": ftype,
        "STRAND": strand,
        "NT_POS": nt_pos,
        "AA_POS": aa_pos,
        "EFFECT": effect,
        "IMPACT": impact,
        "NT_CHANGE": nt_change,
        "AA_CHANGE": aa_change,
        "ID": gff_info.get("bakta_id", ""),
        "LOCUS_TAG": gff_info.get("locus_tag", gene_id_ann),
        "GENE": gff_info.get("gene", ""),
        "PRODUCT": gff_info.get("product", ""),
    }


# ── Canonical merged VCF parsing ───────────────────────────────────────────────

def _opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def gt_to_cell(gt):
    """Map a (biallelic) GT to a table cell: "1" hom-ALT, "0" hom-REF, "-" missing.

    Missing/het classification is delegated to pipeline_helpers.is_missing_gt so
    the canonical-VCF missingness filter and these tables agree on "-"."""
    if is_missing_gt(gt):
        return "-"
    alleles = re.split(r"[/|]", gt)
    return "0" if set(alleles) == {"0"} else "1"


def _info_value(info, key):
    prefix = key + "="
    for kv in info.split(";"):
        if kv.startswith(prefix):
            return kv[len(prefix):]
    return ""


def pick_type(existing, new):
    ep = TYPE_PRIORITY.get(existing, 99)
    np_ = TYPE_PRIORITY.get(new, 99)
    return new if np_ < ep else existing


def parse_merged_vcf(merged_vcf_path, gff_lookup):
    """Read the biallelic, snpEff-annotated canonical VCF.

    Returns (isolates, variants, var_types, var_annots), where variants is
    {(chrom,pos,ref,alt): {iso: cell}}, keyed exactly as the per-row REF/ALT in
    the VCF (so keys match the TreeTime VCF's snp/mnp rows)."""
    isolates = []
    variants = defaultdict(lambda: defaultdict(lambda: "0"))
    var_types = {}
    var_annots = {}

    with _opener(merged_vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                isolates = line.rstrip("\n").split("\t")[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            fmt_keys = fields[8].split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                continue

            key = (chrom, pos, ref, alt)
            var_type = _info_value(info, "TYPE")
            if key in var_types:
                var_types[key] = pick_type(var_types[key], var_type)
            else:
                var_types[key] = var_type

            annot = parse_ann_field(info, gff_lookup, target_alt=alt)
            if key not in var_annots or (
                not var_annots[key]["LOCUS_TAG"] and annot["LOCUS_TAG"]
            ):
                var_annots[key] = annot

            for i, iso in enumerate(isolates):
                cell = fields[9 + i]
                gt = cell.split(":")[gt_idx] if gt_idx < len(cell.split(":")) else "./."
                variants[key][iso] = gt_to_cell(gt)

    if not isolates:
        print(f"ERROR: no #CHROM line in {merged_vcf_path}", file=sys.stderr)
        sys.exit(1)
    print(f"  Total unique variants: {len(variants)}")
    return isolates, variants, var_types, var_annots


# ── Disruptive filter ───────────────────────────────────────────────────────────

def is_disruptive(effect):
    if not effect:
        return False
    if effect in ("synonymous_variant", "missense_variant"):
        return False
    if effect.startswith("intragenic_"):
        return False
    if "_retained_variant" in effect:
        return False
    if effect.startswith("intergenic_"):
        return False
    if effect.startswith("conservative_"):
        return False
    return True


def filter_disruptive(variants, var_annots, isolates):
    disruptive = {}
    for key, gts in variants.items():
        effect = var_annots.get(key, {}).get("EFFECT", "")
        if not is_disruptive(effect):
            continue
        if any(gts.get(iso, "0") == "1" for iso in isolates):
            disruptive[key] = gts
    print(f"  Disruptive variants: {len(disruptive)}")
    return disruptive


# ── Frequency computation ─────────────────────────────────────────────────────

def compute_freq(gts, group_isolates):
    if not group_isolates:
        return 0, 0.0
    ct = sum(1 for iso in group_isolates if gts.get(iso, "0") == "1")
    total = sum(1 for iso in group_isolates if gts.get(iso, "0") != "-")
    frq = ct / total if total > 0 else 0.0
    return ct, frq


# ── Output ────────────────────────────────────────────────────────────────────

def sort_key(variant_key):
    chrom, pos, ref, alt = variant_key
    return (chrom, int(pos), ref, alt)


def write_tables(sorted_keys, variants_dict, var_types, var_annots,
                 isolates, pr_isolates, po_isolates,
                 gt_path, annot_path, label):
    header_gt = ["CHR", "POS", "TYPE", "REF", "ALT"] + isolates
    header_annot = (["CHR", "POS", "TYPE", "REF", "ALT"] + ANNOT_COLS +
                    ["PR_CT", "PO_CT", "FRQ", "PR_FRQ", "PO_FRQ"])

    with open(gt_path, "w", newline="") as fg, open(annot_path, "w", newline="") as fa:
        wg = csv.writer(fg, delimiter="\t")
        wa = csv.writer(fa, delimiter="\t")
        wg.writerow(header_gt)
        wa.writerow(header_annot)

        for key in sorted_keys:
            chrom, pos, ref, alt = key
            gts = variants_dict[key]
            t = var_types.get(key, "")
            annot = var_annots.get(key, {col: "" for col in ANNOT_COLS})
            _, frq = compute_freq(gts, isolates)
            pr_ct, pr_frq = compute_freq(gts, pr_isolates)
            po_ct, po_frq = compute_freq(gts, po_isolates)

            wg.writerow(
                [chrom, pos, t, ref, alt] + [gts.get(iso, "0") for iso in isolates]
            )
            wa.writerow(
                [chrom, pos, t, ref, alt] +
                [annot.get(col, "") for col in ANNOT_COLS] +
                [pr_ct, po_ct, f"{frq:.4f}", f"{pr_frq:.4f}", f"{po_frq:.4f}"]
            )

    print(f"  Wrote {label}: {gt_path.name} + {annot_path.name} ({len(sorted_keys)} variants)")


def write_annot_table(sorted_keys, variants, var_types, var_annots,
                      isolates, pr_isolates, po_isolates, annot_path, label):
    header = (["CHR", "POS", "TYPE", "REF", "ALT"] + ANNOT_COLS +
              ["PR_CT", "PO_CT", "FRQ", "PR_FRQ", "PO_FRQ"])

    with open(annot_path, "w", newline="") as fa:
        wa = csv.writer(fa, delimiter="\t")
        wa.writerow(header)
        for key in sorted_keys:
            chrom, pos, ref, alt = key
            gts = variants[key]
            t = var_types.get(key, "")
            annot = var_annots.get(key, {col: "" for col in ANNOT_COLS})
            _, frq = compute_freq(gts, isolates)
            pr_ct, pr_frq = compute_freq(gts, pr_isolates)
            po_ct, po_frq = compute_freq(gts, po_isolates)
            wa.writerow(
                [chrom, pos, t, ref, alt] +
                [annot.get(col, "") for col in ANNOT_COLS] +
                [pr_ct, po_ct, f"{frq:.4f}", f"{pr_frq:.4f}", f"{po_frq:.4f}"]
            )

    print(f"  Wrote {label}: {annot_path.name} ({len(sorted_keys)} variants)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate variant tables from the canonical merged multi-sample VCF."
    )
    parser.add_argument("--pop", required=True,
                        help="Population name (e.g. IA01)")
    parser.add_argument("--merged-vcf", required=True,
                        help="Biallelic, snpEff-annotated canonical VCF "
                             "(bcftools norm -m- of <POP>.merged.vcf.gz, then snpEff ann)")
    parser.add_argument("--gff", required=True,
                        help="Path to reference GFF3 for annotation lookup")
    parser.add_argument("--out-dir", type=str, default=None,
                        help="Output directory (default: <merged-vcf dir>/../pop_tables)")
    args = parser.parse_args()

    pop_name = args.pop
    merged_vcf = Path(args.merged_vcf)
    gff_path = Path(args.gff)

    print(f"Population: {pop_name}")
    print(f"Canonical VCF: {merged_vcf}")

    print(f"Loading GFF: {gff_path}")
    gff_lookup = load_gff_lookup(gff_path)
    print(f"  Loaded {len(gff_lookup)} locus tags from GFF")

    print("Reading genotypes/TYPE/annotation from canonical VCF...")
    isolates, variants, var_types, var_annots = parse_merged_vcf(merged_vcf, gff_lookup)
    print(f"Found {len(isolates)} isolates in {merged_vcf.name}")

    pr_isolates, po_isolates = classify_isolates(pop_name, isolates)

    out_dir = (Path(args.out_dir) if args.out_dir
               else merged_vcf.resolve().parent.parent / "pop_tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    # All variants with at least one ALT call
    called_variants = {k: gts for k, gts in variants.items()
                       if any(gts.get(iso, "0") == "1" for iso in isolates)}
    print(f"  Variants with at least one ALT call: {len(called_variants)} / {len(variants)}")
    all_sorted_keys = sorted(called_variants.keys(), key=sort_key)
    write_annot_table(all_sorted_keys, called_variants, var_types, var_annots,
                      isolates, pr_isolates, po_isolates,
                      out_dir / "var.annot.tab", "all variants")

    # SNPs: snp/mnp, plus complex variants that are missense/synonymous (substitutions
    # masquerading as complex due to neighboring variation). len(REF)==len(ALT) is
    # verified as a sanity check. These complex variants are excluded from the INDEL
    # set so there is no overlap between snp.* and indel.* outputs.
    def is_snp_like_complex(key):
        if var_types.get(key, "") != "complex":
            return False
        effect = var_annots.get(key, {}).get("EFFECT", "")
        if effect not in ("missense_variant", "synonymous_variant"):
            return False
        _, _, ref, alt = key
        if len(ref) != len(alt):
            print(f"  WARNING: complex {effect} at {key[0]}:{key[1]} has "
                  f"len(REF)={len(ref)} != len(ALT)={len(alt)}; excluding from SNP set",
                  file=sys.stderr)
            return False
        return True

    snp_keys = [k for k in all_sorted_keys
                if var_types.get(k, "") in ("snp", "mnp") or is_snp_like_complex(k)]
    print(f"  SNPs: {len(snp_keys)}")
    write_tables(snp_keys, called_variants, var_types, var_annots,
                 isolates, pr_isolates, po_isolates,
                 out_dir / "snp.gt.tab", out_dir / "snp.annot.tab", "SNPs/MNPs")

    # INDELs: ins/del/complex, excluding the complex variants reclassified as SNPs
    indel_keys = [k for k in all_sorted_keys
                  if var_types.get(k, "") in ("ins", "del", "complex")
                  and not is_snp_like_complex(k)]
    print(f"  INDELs: {len(indel_keys)}")
    write_tables(indel_keys, called_variants, var_types, var_annots,
                 isolates, pr_isolates, po_isolates,
                 out_dir / "indel.gt.tab", out_dir / "indel.annot.tab", "INDELs")

    # Disruptive variants
    print("Filtering disruptive variants...")
    disruptive = filter_disruptive(variants, var_annots, isolates)
    disruptive_keys = sorted(disruptive.keys(), key=sort_key)
    write_annot_table(disruptive_keys, disruptive, var_types, var_annots,
                      isolates, pr_isolates, po_isolates,
                      out_dir / "disruptive.annot.tab", "disruptive variants")

    print("Done!")


if __name__ == "__main__":
    main()
