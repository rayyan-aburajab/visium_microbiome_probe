#!/usr/bin/env python3
from collections import Counter, defaultdict
from pathlib import Path
from itertools import groupby
import re
import click
import pandas as pd
from Bio import SeqIO

LHS = 'CCTTGGCACCCGAGAATTCCA'
RHS = 'A' * 30
NUCLEOTIDES = set("AGTC")
MAX_GC = 0.72
MIN_GC = 0.44
MAX_REPEATS = 3

@click.command()
@click.option('--input_dir', required=True, help='Directory containing input FASTA files')
@click.option('--output_dir', required=True, help='Directory to write probe TSVs to')
@click.option('--output_excel', required=True, help='Excel file with one sheet per genus')
@click.option('--output_fasta_dir', required=True, help='Directory to store FASTA files per genus')
def main(input_dir, output_dir, output_excel, output_fasta_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_fasta_dir = Path(output_fasta_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_fasta_dir.mkdir(parents=True, exist_ok=True)

    for fasta_file in sorted(input_dir.glob("*.fa")):
        species_name = fasta_file.stem
        output_file = output_dir / f"{species_name}_probes.tsv"

        candidates = defaultdict(list)
        for record in SeqIO.parse(fasta_file, "fasta"):
            rc_seq = str(record.seq.reverse_complement().upper())
            for i in range(len(rc_seq) - 50 + 1):
                left = rc_seq[i:i+25]
                right = rc_seq[i+25:i+50]
                if passes_qc(left, right):
                    candidates[record.id].append((left, right, i))

        selected = defaultdict(list)
        for seq_id, probes in candidates.items():
            last_pos = None
            for probe in probes:
                if last_pos and abs(probe[2] - last_pos) < 50:
                    continue
                selected[seq_id].append(probe)
                last_pos = probe[2]

        final_probes = [p for probe_list in selected.values() for p in probe_list]
        if not final_probes:
            print(f"Warning: {species_name}: no valid probes generated.")
            continue

        with open(output_file, 'w') as out:
            out.write("#id\tpos\thyb_region\tprobe\n")
            for i, (left, right, position) in enumerate(final_probes, 1):
                probe_id = f"{species_name}_probe{i}"
                out.write(f"{probe_id}_left\t{position}\t{left}\t{LHS + left}\n")
                out.write(f"{probe_id}_right\t{position}\t{right}\t{right + RHS}\n")

    # Aggregate all TSVs into Excel and FASTA
    pattern = r"^(?P<species>[A-Za-z0-9_]+)_probe(?P<number>\d+)_(?P<tag>left|right)$"
    genus_groups = {}

    for file in sorted(output_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t')
        df.columns = df.columns.str.lstrip('#')
        species = file.stem.replace("_probes", "")
        genus = species.split("_")[0]

        extracted = df['id'].str.extract(pattern)
        df = pd.concat([df, extracted], axis=1)

        left_df = df[df['tag'] == 'left'].copy()
        right_df = df[df['tag'] == 'right'].copy()

        merged = pd.merge(left_df, right_df, on=['species', 'number'], suffixes=('_left', '_right'))
        merged['probe_id'] = merged['species'] + "_probe" + merged['number']
        merged['combined_probe'] = merged['hyb_region_left'] + merged['hyb_region_right']

        genus_groups.setdefault(genus, []).append(merged)

    excel_writer = pd.ExcelWriter(output_excel, engine='xlsxwriter')
    for genus, dfs in genus_groups.items():
        df_all = pd.concat(dfs, ignore_index=True)
        df_all['pos'] = df_all['pos_left']

        excel_df = df_all[['probe_id', 'pos', 'hyb_region_left', 'hyb_region_right', 'combined_probe']]
        excel_df.to_excel(excel_writer, sheet_name=genus[:31], index=False)

        fasta_path = output_fasta_dir / f"{genus}.fa"
        with open(fasta_path, 'w') as f:
            for _, row in df_all.iterrows():
                f.write(f">{row['probe_id']}\n{row['combined_probe']}\n")

    excel_writer.close()


def passes_qc(left, right):
    for seq in [left, right]:
        ct = Counter(seq)
        if not set(ct).issubset(NUCLEOTIDES):
            return False
        gc = (ct['G'] + ct['C']) / sum(ct.values())
        if not (MIN_GC <= gc <= MAX_GC):
            return False
        if not all(sum(1 for _ in group) <= MAX_REPEATS for _, group in groupby(seq)):
            return False
    return left[-1] == 'T'

if __name__ == '__main__':
    main()