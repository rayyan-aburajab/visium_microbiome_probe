#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import click
import re

@click.command()
@click.option('--input_dir', required=True, help="Directory with *_probes.tsv files")
@click.option('--output_excel', required=True, help="Excel file with probes split by sheet")
@click.option('--output_fasta', required=True, help="FASTA file of hybrid-only probes")
def main(input_dir, output_excel, output_fasta):
    input_dir = Path(input_dir)
    fasta_lines = []
    excel_writer = pd.ExcelWriter(output_excel, engine='xlsxwriter')

    for file in sorted(input_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t', header=0)
        df.columns = df.columns.str.lstrip('#')

        species = file.stem.replace("_probes", "").replace(" ", "_")

        # Extract metadata from IDs
        df[['gene', 'protein_id', 'location', 'tag']] = df['id'].str.extract(
            r'([^_]+)_([^_]+)_([^_]+)_(lhs|rhs)-\d+'
        )

        lhs_df = df[df['tag'] == 'lhs'].copy()
        rhs_df = df[df['tag'] == 'rhs'].copy()

        lhs_df['hyb_region'] = lhs_df['hyb_region'].str.upper()
        rhs_df['hyb_region'] = rhs_df['hyb_region'].str.upper()

        # Merge LHS and RHS probes
        merged = pd.merge(lhs_df, rhs_df,
                          on=['pos', 'gene', 'protein_id'],
                          suffixes=('_lhs', '_rhs'))
        merged['combined_probe'] = merged['hyb_region_lhs'] + merged['hyb_region_rhs']
        merged['species'] = species

        # Select and reorder columns for Excel
        output_df = merged[[
            'species', 'gene', 'protein_id',
            'id_lhs', 'id_rhs', 'pos',
            'hyb_region_lhs', 'hyb_region_rhs', 'combined_probe'
        ]]

        # Write this species' sheet
        output_df.to_excel(excel_writer, sheet_name=species[:31], index=False)

        # FASTA output
        for _, row in output_df.iterrows():
            fasta_lines.append(f">{row['id_lhs']}|{row['id_rhs']}")
            fasta_lines.append(row['combined_probe'])

    # Save Excel
    excel_writer.close()

    # Write FASTA
    with open(output_fasta, "w") as f:
        f.write("\n".join(fasta_lines) + "\n")

    print(f"✅ Excel file written: {output_excel}")
    print(f"✅ FASTA file written: {output_fasta}")

if __name__ == "__main__":
    main()
