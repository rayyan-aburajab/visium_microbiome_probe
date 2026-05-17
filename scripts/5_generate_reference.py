#!/usr/bin/env python3
from pathlib import Path
import re
import click
import pandas as pd

@click.command()
@click.option('--input', 'input_file', required=True, help='probe_qc.tsv from step 4')
@click.option('--output_dir', required=True, help='Directory to write reference files')
@click.option('--output_prefix', required=True, help='Prefix for output files')
def main(input_file, output_dir, output_prefix):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_out = output_dir / f"{output_prefix}.fa"
    gtf_out = output_dir / f"{output_prefix}.gtf"
    meta_out = output_dir / f"{output_prefix}_metadata.tsv"

    df = pd.read_csv(input_file, sep='\t')

    with open(fasta_out, 'w') as fasta, \
         open(gtf_out, 'w') as gtf, \
         open(meta_out, 'w') as meta:

        meta.write("probe_id\tsequence\tcomposite_id\tincluded\ttranscript_type\n")

        for _, row in df.iterrows():
            probe_id = row['probe_id']
            gene_name = row['gene_name']
            gene_symbol = row['gene_symbol']
            composite_group = row['composite_group']
            sequence = row['hyb_region_left'] + row['hyb_region_right']
            length = len(sequence)

            # FASTA
            fasta.write(f">{probe_id}\n{sequence}\n")

            # GTF
            attr = f'gene_id "{gene_name}"; gene_name "{gene_symbol}"; transcript_id "{probe_id}";'
            gtf.write(f"{probe_id}\tCustom\ttranscript\t1\t{length}\t.\t+\t.\t{attr}\n")
            gtf.write(f"{probe_id}\tCustom\texon\t1\t{length}\t.\t+\t.\t{attr}\n")

            # Metadata
            match = re.search(r'(\d+)$', probe_id)
            if not match:
                print(f"Warning: could not parse probe number from {probe_id}, skipping.")
                continue
            probe_number = int(match.group(1))
            padded_num = f"{probe_number:03}"
            composite_id = f"{gene_name}|{gene_symbol}|{composite_group}{gene_symbol.lower()}{padded_num}"
            meta.write(f"{probe_id}\t{sequence}\t{composite_id}\tTRUE\tunspliced\n")


if __name__ == '__main__':
    main()