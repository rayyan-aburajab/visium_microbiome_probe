#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from Bio.Seq import Seq
import click

@click.command()
@click.option('--input_dir', required=True, help="Directory with *_probes.tsv files")
@click.option('--output_file', required=True, help="TSV with combined hybrid sequences")
@click.option('--output_fasta', required=True, help="FASTA file of hybrid-only probes")
def main(input_dir, output_file, output_fasta):
    input_dir = Path(input_dir)
    all_combined = []
    fasta_lines = []

    for file in sorted(input_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t', header=0)
        df.columns = df.columns.str.lstrip('#')  # Strip '#' from column names

        species = file.stem.replace("_probes", "").replace(" ", "_")

        lhs_df = df[df['id'].str.contains('_lhs-')].copy()
        rhs_df = df[df['id'].str.contains('_rhs-')].copy()

        lhs_df['hyb_region'] = lhs_df['hyb_region'].str.upper()
        rhs_df['hyb_region'] = rhs_df['hyb_region'].str.upper()

        # Merge LHS and RHS by position
        merged = pd.merge(lhs_df, rhs_df, on='pos', suffixes=('_lhs', '_rhs'))
        merged['combined_probe'] = merged['hyb_region_lhs'] + merged['hyb_region_rhs']
        merged['species'] = species

        # FASTA-formatted column (now includes both LHS and RHS IDs)
        merged['fasta_entry'] = [
            f">{lhs_id}|{rhs_id}\n{seq}"
            for lhs_id, rhs_id, seq in zip(merged['id_lhs'], merged['id_rhs'], merged['combined_probe'])
        ]

        all_combined.append(merged)
        for _, row in merged.iterrows():
            fasta_lines.append(f">{row['id_lhs']}|{row['id_rhs']}")
            fasta_lines.append(row['combined_probe'])

    # Output TSV
    combined_df = pd.concat(all_combined, ignore_index=True)
    combined_df = combined_df[['species', 'id_lhs', 'id_rhs', 'pos', 'hyb_region_lhs', 'hyb_region_rhs', 'combined_probe', 'fasta_entry']]
    combined_df.to_csv(output_file, sep='\t', index=False)

    # Output FASTA
    with open(output_fasta, "w") as f:
        f.write("\n".join(fasta_lines) + "\n")

    print(f"✅ Combined TSV written: {output_file}")
    print(f"✅ FASTA file written:   {output_fasta}")

if __name__ == "__main__":
    main()
