#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import click
import re

@click.command()
@click.option('--input_dir', required=True, help="Directory with *_probes.tsv files")
@click.option('--output_excel', required=True, help="Excel file with probes split by sheet")
@click.option('--output_fasta_dir', required=True, help="Directory to save per-virus FASTA files")
def main(input_dir, output_excel, output_fasta_dir):
    input_dir = Path(input_dir)
    output_fasta_dir = Path(output_fasta_dir)
    output_fasta_dir.mkdir(parents=True, exist_ok=True)

    fasta_dict = {}
    excel_writer = pd.ExcelWriter(output_excel, engine='xlsxwriter')

    for file in sorted(input_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t')
        df.columns = df.columns.str.lstrip('#')

        virus = file.stem.replace("_probes", "").replace(" ", "_")

        # Parse ID like: EBV_LMP1_probe1_lhs
        pattern = r"(?P<virus>[^_]+)_(?P<gene>[^_]+)_(?P<probe_num>probe\d+)_(?P<tag>lhs|rhs)"
        extracted = df['id'].str.extract(pattern)
        df = pd.concat([df, extracted], axis=1)

        lhs_df = df[df['tag'] == 'lhs'].copy()
        rhs_df = df[df['tag'] == 'rhs'].copy()

        lhs_df['hyb_region'] = lhs_df['hyb_region'].str.upper()
        rhs_df['hyb_region'] = rhs_df['hyb_region'].str.upper()

        # Merge LHS and RHS on matching probe ID fields
        merged = pd.merge(lhs_df, rhs_df,
                          on=['virus', 'gene', 'probe_num'],
                          suffixes=('_lhs', '_rhs'))
        
        merged['protein'] = merged['protein_lhs']
        merged['combined_probe'] = merged['hyb_region_lhs'] + merged['hyb_region_rhs']
        merged['probe_id'] = merged['virus'] + "_" + merged['gene'] + "_" + merged['probe_num']

        # Excel sheet (per virus)
        sheet_df = merged[[
            'probe_id', 'gene', 'protein', 'pos_lhs',
            'hyb_region_lhs', 'hyb_region_rhs', 'combined_probe'
        ]].rename(columns={'pos_lhs': 'pos'})
        sheet_df.to_excel(excel_writer, sheet_name=virus[:31], index=False)

        # FASTA lines
        fasta_lines = []
        for _, row in merged.iterrows():
            fasta_lines.append(f">{row['probe_id']}")
            fasta_lines.append(row['combined_probe'])

        # Write per-virus FASTA
        with open(output_fasta_dir / f"{virus}.fa", "w") as f:
            f.write("\n".join(fasta_lines) + "\n")

    excel_writer.close()
    print(f"Excel file written: {output_excel}")
    print(f"Per-virus FASTA files written to: {output_fasta_dir}")

if __name__ == "__main__":
    main()
