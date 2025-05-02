#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Combine viral region FASTAs into a single file per virus.")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing individual viral region FASTAs")
    parser.add_argument("--output_dir", required=True, help="Directory to save combined FASTAs")
    return parser.parse_args()

def group_fastas_by_virus(fasta_files):
    virus_groups = defaultdict(list)
    for fasta in fasta_files:
        basename = os.path.basename(fasta)
        if "_" in basename:
            virus_prefix = basename.split("_")[0]
            virus_groups[virus_prefix].append(fasta)
    return virus_groups

def combine_fastas(virus_groups, output_dir):
    for virus, files in virus_groups.items():
        combined_records = []
        for f in sorted(files):  # keep order consistent
            region_tag = os.path.splitext(os.path.basename(f))[0]  # e.g., HPV16_E6
            for record in SeqIO.parse(f, "fasta"):
                record.id = f"{region_tag}|{record.id}"
                record.description = ""  # Clean up description
                combined_records.append(record)
        output_path = os.path.join(output_dir, f"{virus}_combined.fa")
        SeqIO.write(combined_records, output_path, "fasta")
        print(f"✅ Combined {len(files)} FASTAs into: {output_path} ({len(combined_records)} sequences)")

def main():
    args = parse_args()
    fasta_files = [os.path.join(args.fasta_dir, f) for f in os.listdir(args.fasta_dir)
                   if f.endswith(".fa") or f.endswith(".fasta")]
    virus_groups = group_fastas_by_virus(fasta_files)
    os.makedirs(args.output_dir, exist_ok=True)
    combine_fastas(virus_groups, args.output_dir)

if __name__ == "__main__":
    main()
