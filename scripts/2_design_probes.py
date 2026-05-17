#!/usr/bin/env python3
from collections import Counter, defaultdict
from pathlib import Path
from itertools import groupby
import click
import pandas as pd
from Bio import SeqIO

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
        final_probes = list({left + right: (left, right, position) for left, right, position in final_probes}.values())
        if not final_probes:
            print(f"Warning: {species_name}: no valid probes generated.")
            continue

        with open(output_file, 'w') as out:
            out.write("probe_id\tpos\thyb_region_left\thyb_region_right\thyb_region_full\n")
            for i, (left, right, position) in enumerate(final_probes, 1):
                probe_id = f"{species_name}_probe{i}"
                out.write(f"{probe_id}\t{position}\t{left}\t{right}\t{left + right}\n")

    # Aggregate all TSVs into Excel and FASTA
    genus_groups = {}

    for file in sorted(output_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t')
        species = file.stem.replace("_probes", "")
        genus = species.split("_")[0]
        genus_groups.setdefault(genus, []).append(df)

    excel_writer = pd.ExcelWriter(output_excel, engine='xlsxwriter')
    for genus, dfs in genus_groups.items():
        df_all = pd.concat(dfs, ignore_index=True)

        excel_df = df_all[['probe_id', 'pos', 'hyb_region_left', 'hyb_region_right', 'hyb_region_full']]
        excel_df.to_excel(excel_writer, sheet_name=genus[:31], index=False)

        fasta_path = output_fasta_dir / f"{genus}.fa"
        with open(fasta_path, 'w') as f:
            for _, row in df_all.iterrows():
                f.write(f">{row['probe_id']}\n{row['hyb_region_full']}\n")

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