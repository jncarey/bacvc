"""Shared helpers across pipeline step 04 scripts.

- build_callable_mask / is_callable:  bit-packed boolean mosdepth-coverage mask
  (1 bit per genome position per sample, ~genome_bp/8 bytes total). Used by
  04_generate_variant_vcf.py, which holds masks for *every* isolate in the
  population at once and needs to scale to ~1000 samples. Collapses
  "depth >= MASK_MIN_DP" into a single bit at build time.

- load_repeat_bed / overlaps_repeat:  self-alignment repeat-region exclusion,
  shared by 04_generate_variant_vcf.py and 04_generate_variant_tables.py.

- is_missing_gt:  the single definition of a "-" (missing) genotype cell, shared
  by the canonical-VCF missingness filter (04_generate_variant_vcf.py) and the
  table generator's gt_to_cell (04_generate_variant_tables.py).

- group_isolates_by_visit:  splits isolates into groups keyed by visit label
  (V1, V1A, V2, ...) parsed from each isolate's name, shared by the table
  generator and the canonical pass's per-site missingness QC.
"""

import bisect
import gzip
import re
from collections import defaultdict


# ── Bit-packed callable-region mask ──────────────────────────────────────────

def build_callable_mask(bed_gz_path, chrom_lengths, min_dp):
    """Read mosdepth per-base bedGraph and return a bit-packed callable mask.

    Returns {chrom: bytearray}. Bit i in chrom's bytearray is set iff the
    0-based position i had mosdepth depth >= min_dp.

    Memory: ceil(chrom_len/8) bytes per chromosome per sample, e.g. ~850 KB
    for a 6.8 Mb bacterial genome -- ~270x smaller than load_depth_index().
    """
    masks = {c: bytearray((L + 7) // 8) for c, L in chrom_lengths.items()}
    with gzip.open(bed_gz_path, "rt") as f:
        for line in f:
            chrom, start, end, depth = line.rstrip("\n").split("\t")
            if int(depth) < min_dp:
                continue
            arr = masks.get(chrom)
            if arr is None:
                continue
            _set_range(arr, int(start), int(end))
    return masks


def _set_range(arr, start, end):
    """OR-set bits [start, end) (0-based, half-open) in a bytearray treated
    as a packed bit array (bit i lives in byte i>>3 at position i&7)."""
    if start >= end:
        return
    fb = start >> 3
    lb = (end - 1) >> 3
    if fb == lb:
        # Whole range falls within a single byte.
        arr[fb] |= ((1 << (end - start)) - 1) << (start & 7)
        return
    # Multi-byte: partial first byte, full middle bytes, partial last byte.
    arr[fb] |= (0xFF << (start & 7)) & 0xFF
    if lb > fb + 1:
        arr[fb + 1:lb] = b'\xff' * (lb - fb - 1)
    arr[lb] |= (1 << (((end - 1) & 7) + 1)) - 1


def is_callable(masks, chrom, pos_1based):
    """True iff the 1-based POS is set in the chrom's callable mask. False
    when the chrom is unknown (treated as not callable -> ./. downstream)."""
    arr = masks.get(chrom)
    if arr is None:
        return False
    p0 = pos_1based - 1
    return bool(arr[p0 >> 3] >> (p0 & 7) & 1)


# ── Self-alignment repeat-region exclusion ───────────────────────────────────
# A single overlap implementation shared by the canonical-VCF builder
# (04_generate_variant_vcf.py) and the table builder (04_generate_variant_tables.py),
# so repeat exclusion is identical across every step-04 artifact.

def load_repeat_bed(bed_path):
    """Load a (0-based, half-open) repeat BED into {chrom: (starts, ends)} with
    starts sorted for bisect lookup. Intervals are assumed non-overlapping and
    coordinate-sorted (bedtools merge output), but we sort defensively."""
    chrom_data = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom_data[parts[0]].append((int(parts[1]), int(parts[2])))
    repeat_idx = {}
    for chrom, ivs in chrom_data.items():
        ivs.sort()
        repeat_idx[chrom] = ([s for s, _ in ivs], [e for _, e in ivs])
    return repeat_idx


def overlaps_repeat(repeat_idx, chrom, pos_1based, ref_len=1):
    """True iff the variant span [pos-1, pos-1+ref_len) (0-based, half-open)
    intersects any repeat interval on chrom. Reduces to POS-in-interval for
    SNPs (ref_len=1). False when the chrom has no intervals."""
    entry = repeat_idx.get(chrom)
    if entry is None:
        return False
    starts, ends = entry
    v_start = pos_1based - 1
    v_end = v_start + max(ref_len, 1)
    # Right-most interval whose start is < v_end; if it (or an earlier one)
    # ends past v_start, the spans intersect.
    i = bisect.bisect_left(starts, v_end) - 1
    return i >= 0 and ends[i] > v_start


# ── Genotype "missing" predicate ─────────────────────────────────────────────
# Single definition of a "-" cell, so the canonical-VCF missingness filter and
# the table generator's gt_to_cell agree on which genotypes count as missing.

def is_missing_gt(gt):
    """True iff a (biallelic) GT should be treated as missing ("-").

    Missing = any allele is "." (./., ., .|.), OR a het-with-ref (allele set
    contains "0" and a nonzero allele, e.g. 0/1; masked upstream but treated as
    missing if it survives). Hom-ref ("0") and hom-alt (1/1, 2/2, ...) are not
    missing; 1/2 (two distinct nonzero alleles, no ref) is not missing."""
    if gt in (".", "./.", ".|."):
        return True
    alleles = re.split(r"[/|]", gt)
    if any(a == "." for a in alleles):
        return True
    allele_set = set(alleles)
    return "0" in allele_set and len(allele_set) > 1


# ── Isolate grouping by visit ────────────────────────────────────────────────

def group_isolates_by_visit(patient_id, isolates, log=print):
    """Group isolates by visit label (V1, V1A, V2, ...) parsed directly from
    each isolate's real name, e.g. "009-007_V1A_10" -> "V1A".

    Each isolate's name must literally start with "{patient_id}_V<visit>_" to
    be included -- this both extracts the visit label and confirms the
    isolate actually belongs to patient_id, so isolates from a different
    patient can never be silently merged into this patient's visit groups.

    Returns {visit_label: [isolates]}. Isolates that don't match (wrong
    patient, or no parseable visit label) are reported via `log` and excluded
    from every group."""
    pattern = re.compile(r'^' + re.escape(patient_id) + r'_V(\d+[A-Za-z]?)_')

    groups = defaultdict(list)
    unmatched = []
    for iso in isolates:
        m = pattern.match(iso)
        if m:
            groups[f"V{m.group(1)}"].append(iso)
        else:
            unmatched.append(iso)

    for visit in sorted(groups):
        log(f"  {visit}: {len(groups[visit])} isolate(s)")
    if unmatched:
        log(f"  WARNING: {len(unmatched)} isolates didn't match patient "
            f"'{patient_id}' + a parseable visit label, and will be excluded "
            f"from frequency calculations:")
        for iso in unmatched[:5]:
            log(f"    {iso}")
        if len(unmatched) > 5:
            log(f"    ... and {len(unmatched) - 5} more")

    return dict(groups)
