#!/usr/bin/env python3
"""
Convert bakta GenBank (.gbff) to a flat snippy-style GFF3 for snpEff.

Produces CDS-only GFF3 (no gene/mRNA hierarchy) matching the format that
snippy generates and snpEff expects. snpEff uses:
  - ID=  -> Gene_ID in ANN field
  - Name= -> Gene_Name in ANN field

For each CDS feature in the GenBank file:
  - ID = bakta locus_tag (always unique, e.g., IA01_autocycler_00001)
  - locus_tag = UserProtein ID if available (e.g., PA0001), else bakta locus_tag
  - Name = gene name if available, else locus_tag
  - gene = gene name if available
  - product = product description

If gene_name_overrides.tsv exists alongside this script (columns:
bakta_locus_tag, override_locus_tag, gene_name), matching loci get their
locus_tag/gene_name replaced -- e.g. to reconcile a gene split across bakta
annotations under one stable name/ID for cross-population comparison. This is
optional and dataset-specific: with no such file (or no matching rows for
your loci), nothing changes.

Usage:
    python gbff_to_snpeff_gff.py <input.gbff> > genes.gff
"""

import sys
import re
from pathlib import Path


def parse_location(loc_str):
    """Parse GenBank location string -> (start, end, strand).

    Handles: 123..456, complement(123..456), join(1..100,200..300).
    Returns 1-based coordinates as strings.
    """
    strand = "+"
    s = loc_str.strip()
    if s.startswith("complement(") and s.endswith(")"):
        strand = "-"
        s = s[11:-1]
    # Handle join — take outer bounds
    if s.startswith("join(") and s.endswith(")"):
        s = s[5:-1]
    # Remove < and > (partial features)
    s = s.replace("<", "").replace(">", "")
    # Get all coordinate pairs
    coords = re.findall(r"(\d+)\.\.(\d+)", s)
    if not coords:
        return None, None, strand
    start = coords[0][0]
    end = coords[-1][1]
    return start, end, strand


def parse_gbff_features(gbff_path):
    """Parse GenBank file, yield (seqid, features_list) per record."""
    with open(gbff_path) as f:
        seqid = None
        in_features = False
        current_feature = None
        current_qualifiers = {}
        current_qual_key = None
        current_qual_val = None
        current_location = None
        current_ftype = None

        def flush_feature():
            if current_ftype and current_location:
                return (current_ftype, current_location, dict(current_qualifiers))
            return None

        def flush_qualifier():
            nonlocal current_qual_key, current_qual_val
            if current_qual_key:
                val = current_qual_val.strip('"') if current_qual_val else ""
                if current_qual_key in current_qualifiers:
                    existing = current_qualifiers[current_qual_key]
                    if isinstance(existing, list):
                        existing.append(val)
                    else:
                        current_qualifiers[current_qual_key] = [existing, val]
                else:
                    current_qualifiers[current_qual_key] = val
                current_qual_key = None
                current_qual_val = None

        features = []

        for line in f:
            line = line.rstrip("\n")

            if line.startswith("LOCUS"):
                # New record — flush previous
                if seqid and features:
                    yield seqid, features
                parts = line.split()
                seqid = parts[1] if len(parts) > 1 else None
                features = []
                in_features = False
                current_feature = None
                continue

            if line.startswith("FEATURES"):
                in_features = True
                continue

            if line.startswith("ORIGIN") or line.startswith("//"):
                # End of features section
                if in_features:
                    flush_qualifier()
                    feat = flush_feature()
                    if feat:
                        features.append(feat)
                in_features = False
                if line.startswith("//") and seqid and features:
                    yield seqid, features
                    features = []
                    seqid = None
                continue

            if not in_features:
                continue

            # Feature line: "     CDS             complement(6536..7309)"
            # Qualifier line: "                     /gene="dnaA""
            # Continuation:   "                     AAGGCC..."

            if len(line) >= 6 and line[5] != " " and line[:5] == "     ":
                # New feature key
                flush_qualifier()
                feat = flush_feature()
                if feat:
                    features.append(feat)

                parts = line.strip().split(None, 1)
                current_ftype = parts[0]
                current_location = parts[1] if len(parts) > 1 else ""
                current_qualifiers = {}
                current_qual_key = None
                current_qual_val = None

            elif line.startswith("                     /"):
                # New qualifier
                flush_qualifier()
                qline = line.strip()[1:]  # remove leading /
                if "=" in qline:
                    current_qual_key, current_qual_val = qline.split("=", 1)
                else:
                    current_qual_key = qline
                    current_qual_val = ""

            elif line.startswith("                     "):
                # Continuation of qualifier value (e.g., translation)
                if current_qual_key:
                    current_qual_val += line.strip()

            elif line.startswith("     ") and current_location:
                # Continuation of location
                current_location += line.strip()

        # Final flush
        if in_features:
            flush_qualifier()
            feat = flush_feature()
            if feat:
                features.append(feat)
        if seqid and features:
            yield seqid, features


def extract_user_protein(qualifiers):
    """Extract UserProtein ID from db_xref qualifiers."""
    db_xref = qualifiers.get("db_xref", [])
    if isinstance(db_xref, str):
        db_xref = [db_xref]
    for xref in db_xref:
        if xref.startswith('"'):
            xref = xref.strip('"')
        if xref.startswith("UserProtein:"):
            return xref[12:]
    return ""


def load_overrides(tsv_path):
    """Load bakta_locus_tag -> (override_locus_tag, gene_name) from a TSV with
    header columns bakta_locus_tag, override_locus_tag, gene_name. Returns {}
    if tsv_path doesn't exist -- overrides are optional and dataset-specific."""
    if not tsv_path.is_file():
        return {}
    overrides = {}
    with open(tsv_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header != ["bakta_locus_tag", "override_locus_tag", "gene_name"]:
            print(f"ERROR: unexpected header in {tsv_path}: {header}",
                  file=sys.stderr)
            sys.exit(1)
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            bakta_id, override_lt, gene_name = line.split("\t")
            overrides[bakta_id] = (override_lt, gene_name)
    return overrides


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.gbff>", file=sys.stderr)
        sys.exit(1)

    gbff_path = sys.argv[1]
    feature_types = {"CDS"}

    # Optional, dataset-specific locus-tag/gene-name overrides -- e.g. this
    # cohort's gene_name_overrides.tsv bundles the panaroo cluster
    # mucA~~~rseA~~~algU under LOCUS_TAG=PA0763 (mucA) for cross-population
    # joining, since bakta doesn't set UserProtein for those loci. Irrelevant
    # (empty) for a fresh dataset unless you add your own overrides file.
    overrides_path = Path(__file__).resolve().parent / "gene_name_overrides.tsv"
    OVERRIDES = load_overrides(overrides_path)

    print("##gff-version 3")

    for seqid, features in parse_gbff_features(gbff_path):
        # Collect sequence length from the last feature for ##sequence-region
        max_end = 0
        for ftype, location, quals in features:
            start, end, strand = parse_location(location)
            if end:
                max_end = max(max_end, int(end))
        if max_end > 0:
            print(f"##sequence-region {seqid} 1 {max_end}")

        seen_ids = {}
        for ftype, location, quals in features:
            if ftype not in feature_types:
                continue

            start, end, strand = parse_location(location)
            if not start or not end:
                continue

            bakta_locus_tag = quals.get("locus_tag", "").strip('"')
            if not bakta_locus_tag:
                continue

            gene_name = quals.get("gene", "").strip('"')
            product = quals.get("product", "").strip('"')
            user_protein = extract_user_protein(quals)
            codon_start = quals.get("codon_start", "").strip('"')

            # Apply manual overrides if present
            if bakta_locus_tag in OVERRIDES:
                override_lt, override_gene = OVERRIDES[bakta_locus_tag]
                if not user_protein:
                    user_protein = override_lt
                gene_name = override_gene

            # locus_tag = UserProtein if available, else bakta locus_tag
            locus_tag = user_protein if user_protein else bakta_locus_tag
            # Name = gene name if available, else locus_tag
            name = gene_name if gene_name else locus_tag
            # ID = locus_tag (snpEff uses this as Gene_Name in ANN)
            # Deduplicate: IS elements share locus_tags; append _N for copies
            feat_id = locus_tag
            if feat_id in seen_ids:
                seen_ids[feat_id] += 1
                feat_id = f"{feat_id}_{seen_ids[feat_id]}"
            else:
                seen_ids[feat_id] = 0

            # Build GFF attributes
            attrs = [f"ID={feat_id}"]
            attrs.append(f"Name={name}")
            if gene_name:
                attrs.append(f"gene={gene_name}")
            attrs.append(f"locus_tag={locus_tag}")
            attrs.append(f"bakta_id={bakta_locus_tag}")
            if product:
                # URL-encode commas for GFF3 compliance
                product_safe = product.replace(",", "%2C")
                attrs.append(f"product={product_safe}")
            if codon_start:
                attrs.append(f"codon_start={codon_start}")

            # Frame: 0 for CDS, . for others
            phase = "0" if ftype == "CDS" else "."

            print("\t".join([
                seqid, "bakta", ftype, start, end, ".", strand, phase,
                ";".join(attrs)
            ]))


if __name__ == "__main__":
    main()
