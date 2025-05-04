#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import click
import re

@click.command()
@click.option('--input_dir', required=True, help="Directory with *_probes.tsv files")
@click.option('--output_excel', required=True, help="Excel file with one sheet per genus")
@click.option('--output_fasta_dir', required=True, help="Directory to store FASTA files per genus")
def main(input_dir, output_excel, output_fasta_dir):
    input_dir = Path(input_dir)
    output_fasta_dir = Path(output_fasta_dir)
    output_fasta_dir.mkdir(parents=True, exist_ok=True)
    excel_writer = pd.ExcelWriter(output_excel, engine='xlsxwriter')

    genus_groups = {}

    for file in sorted(input_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t')
        df.columns = df.columns.str.lstrip('#')
        species = file.stem.replace("_probes", "")
        genus = species.split("_")[0]

        # Extract probe number and tag (lhs/rhs)
        pattern = r"^(?P<species>[A-Za-z0-9_]+)_probe(?P<number>\d+)_(?P<tag>lhs|rhs)$"
        extracted = df['id'].str.extract(pattern)
        df = pd.concat([df, extracted], axis=1)

        lhs_df = df[df['tag'] == 'lhs'].copy()
        rhs_df = df[df['tag'] == 'rhs'].copy()

        lhs_df['hyb_region'] = lhs_df['hyb_region'].str.upper()
        rhs_df['hyb_region'] = rhs_df['hyb_region'].str.upper()

        merged = pd.merge(lhs_df, rhs_df, on=['species', 'number'], suffixes=('_lhs', '_rhs'))
        merged['probe_id'] = merged['species'] + "_probe" + merged['number']
        merged['combined_probe'] = merged['hyb_region_lhs'] + merged['hyb_region_rhs']

        # Add to genus group
        genus_groups.setdefault(genus, []).append(merged)

    # Write Excel and FASTAs
    for genus, dfs in genus_groups.items():
        df_all = pd.concat(dfs, ignore_index=True)
        df_all['pos'] = df_all['pos_lhs']

        # Excel columns
        excel_df = df_all[['probe_id', 'pos', 'hyb_region_lhs', 'hyb_region_rhs', 'combined_probe']]
        excel_df.to_excel(excel_writer, sheet_name=genus[:31], index=False)

        # FASTA
        fasta_path = output_fasta_dir / f"{genus}.fa"
        with open(fasta_path, 'w') as f:
            for _, row in df_all.iterrows():
                f.write(f">{row['probe_id']}\n{row['combined_probe']}\n")

    excel_writer.close()
    print(f"✅ Excel file written: {output_excel}")
    print(f"✅ FASTA files written to: {output_fasta_dir}")

if __name__ == "__main__":
    main()
