#!/usr/bin/env python3
"""
bed_and_fasta_to_genbank.py

Goal: 
  To convert a FASTA and BED annotation file into a GenBank file.

This code was made for centromeric annotations. Please modify the features as required.
It can optionally extract a specific genomic region (e.g. chrI:3677528-3877528).

Usage:
    python bed_and_fasta_to_genbank.py genome.fasta annotations.bed output.gb [chr:start-end]
"""

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import re

FEATURE_TYPE_MAP = {
    "gene": "gene",
    "tRNA": "tRNA",
    "dh": "misc_feature",
    "dg": "misc_feature",
    "cnt1": "repeat_region",
    "cnt2": "repeat_region",
    "cnt3": "repeat_region",
    "imr_chr1": "centromere",
    "imr_chr2": "centromere",
    "imr_chr3": "centromere",
}


def parse_region(region_str):
    """Parse a region string like 'chrI:3677528-3877528'."""
    m = re.match(r"^(\S+):(\d+)-(\d+)$", region_str)
    if not m:
        sys.exit(f"Invalid region format: {region_str} (expected chr:start-end)")
    chrom, start, end = m.groups()
    return chrom, int(start), int(end)


def bed_to_genbank(fasta_file, bed_file, output_file, region=None):
    # Load all contigs
    records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    for rec in records.values():
        rec.annotations["molecule_type"] = "DNA"
        rec.features = []

    # If a region was provided, parse it
    if region:
        region_chr, region_start, region_end = parse_region(region)
        if region_chr not in records:
            sys.exit(f"Contig {region_chr} not found in FASTA file.")
        region_mode = True
    else:
        region_mode = False

    # Parse BED file
    with open(bed_file) as bed:
        header = bed.readline().strip().split()
        if not header[0].lower().startswith("chrom"):
            bed.seek(0)

        for line in bed:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                sys.stderr.write(f"Skipping malformed line: {line}")
                continue

            chrom, start, end, label = parts[0:4]
            strand = parts[5] if len(parts) > 5 else "+"
            tRNA_type = parts[6] if len(parts) > 6 else "None"
            tRNA_seq = parts[7] if len(parts) > 7 else "None"
            gene_name = parts[8] if len(parts) > 8 else "None"

            if chrom not in records:
                sys.stderr.write(f"Warning: contig '{chrom}' not in FASTA; skipping.\n")
                continue

            start, end = int(start), int(end)
            strand_val = 1 if strand == "+" else -1
            ftype = FEATURE_TYPE_MAP.get(label, "misc_feature")

            # Skip features outside the region if region mode
            if region_mode and chrom == region_chr:
                if end < region_start or start > region_end:
                    continue

            qualifiers = {"label": label}
            if gene_name and gene_name != "None":
                qualifiers["gene"] = gene_name
            if label == "tRNA":
                qualifiers["product"] = f"tRNA-{tRNA_type}({tRNA_seq})"
            elif label == "gene" and gene_name:
                qualifiers["note"] = f"protein-coding gene {gene_name}"

            feature = SeqFeature(
                FeatureLocation(start, end, strand=strand_val),
                type=ftype,
                qualifiers=qualifiers
            )
            records[chrom].features.append(feature)

    # Handle region extraction
    if region_mode:
        rec = records[region_chr]
        # Extract the subsequence and overlapping features
        sub_rec = rec[region_start:region_end]
        sub_rec.id = f"{region_chr}_{region_start}_{region_end}"
        sub_rec.name = f"{region_chr}_{region_start}_{region_end}"
        sub_rec.description = f"Subsequence from {region_chr}:{region_start}-{region_end}"
        SeqIO.write(sub_rec, output_file, "genbank")
        print(f"Extracted region {region_chr}:{region_start}-{region_end} -> {output_file}")
    else:
        # Write all contigs
        SeqIO.write(records.values(), output_file, "genbank")
        print(f"GenBank written: {output_file}")
        print(f"Contigs: {', '.join(records.keys())}")


if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        sys.exit("Usage: python bed_to_genbank.py genome.fasta annotations.bed output.gb [chr:start-end]")
    fasta, bed, out = sys.argv[1:4]
    region = sys.argv[4] if len(sys.argv) == 5 else None
    bed_to_genbank(fasta, bed, out, region)
