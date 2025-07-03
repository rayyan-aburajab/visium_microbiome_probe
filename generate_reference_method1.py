#!/usr/bin/env python3

import csv
import click
import re

@click.command()
@click.option('--input_file', required=True, help='TSV input file with probe_id, gene_name, gene_symbol, hyb_region_lhs, hyb_region_rhs')
@click.option('--output_prefix', required=True, help='Prefix for output files')
def main(input_file, output_prefix):
    fasta_file = f"{output_prefix}.fa"
    gtf_file = f"{output_prefix}_gtf.tsv"
    meta_file = f"{output_prefix}_metadata.csv"

    with open(input_file, newline='') as infile, \
         open(fasta_file, 'w') as fasta_out, \
         open(gtf_file, 'w', newline='') as gtf_out, \
         open(meta_file, 'w', newline='') as meta_out:

        reader = csv.DictReader(infile, delimiter=',')
        print("Detected columns:", reader.fieldnames)
        gtf_writer = csv.writer(gtf_out, delimiter='\t')
        meta_writer = csv.writer(meta_out)
        meta_writer.writerow(['probe_id', 'sequence', 'composite_id', 'included', 'transcript_type'])

        for row in reader:
            probe_id = row['probe_id'].strip()
            gene_name = row['gene_name'].strip()
            gene_symbol = row['gene_symbol'].strip()
            lhs = row['hyb_region_lhs'].strip().upper()
            rhs = row['hyb_region_rhs'].strip().upper()
            sequence = lhs + rhs
            length = len(sequence)

            # FASTA
            fasta_out.write(f">{probe_id}\n{sequence}\n")

            # GTF
            attr = f'gene_id "{gene_name}"; gene_name "{gene_symbol}"; transcript_id "{probe_id}";'
            gtf_writer.writerow([probe_id, 'Custom', 'transcript', 1, length, '.', '+', '.', attr])
            gtf_writer.writerow([probe_id, 'Custom', 'exon', 1, length, '.', '+', '.', attr])

            # Metadata
            match = re.search(r'(\d+)$', probe_id)
            if match:
                probe_number = int(match.group(1))
                padded_num = f"{probe_number:03}"
                composite_id = f"{probe_id}|{gene_symbol}|1{gene_symbol.lower()}{padded_num}"
                meta_writer.writerow([probe_id, sequence, composite_id, 'TRUE', 'unspliced'])
            else:
                continue

    print(f"Files created:\n- {fasta_file}\n- {gtf_file}\n- {meta_file}")

if __name__ == '__main__':
    main()
