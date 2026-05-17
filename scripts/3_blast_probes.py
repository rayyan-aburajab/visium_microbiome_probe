#!/usr/bin/env python3
from pathlib import Path
import re
import subprocess
import click
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

@click.command()
@click.option('--input_dir', required=True, help='Directory containing probe TSV files')
@click.option('--blastdb', required=True, help='Path to 16S BLAST database (no extension)')
@click.option('--output_dir', required=True, help='Directory to write BLAST outputs')
def main(input_dir, blastdb, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_file = output_dir / "probes.fa"
    blast_out = output_dir / "blast_results.tsv"
    blast_excel = output_dir / "blast_results.xlsx"

    # Combine all probe TSVs into one FASTA for BLAST
    records = []
    probe_lengths = {}
    for file in sorted(input_dir.glob("*_probes.tsv")):
        df = pd.read_csv(file, sep='\t')
        for _, row in df.iterrows():
            rec = SeqRecord(Seq(row['hyb_region_full']), id=row['probe_id'], description="")
            records.append(rec)
            probe_lengths[row['probe_id']] = len(row['hyb_region_full'])
    SeqIO.write(records, fasta_file, "fasta")

    # Run BLAST
    cmd = [
        "blastn",
        "-query", str(fasta_file),
        "-db", blastdb,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
        "-out", str(blast_out),
        "-evalue", "1e-5",
        "-max_target_seqs", "100",
    ]
    subprocess.run(cmd, check=True)

    # Parse BLAST results and add header
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
    ]
    df_blast = pd.read_csv(blast_out, sep='\t', header=None, names=cols)
    df_blast["qcov"] = df_blast.apply(lambda row: row["length"] / probe_lengths[row["qseqid"]], axis=1)
    df_filtered = df_blast[(df_blast["pident"] >= 90.0) & (df_blast["qcov"] >= 0.9)].copy()
    df_filtered["group"] = df_filtered["qseqid"].apply(lambda x: x.split("_")[0])
    df_filtered["species"] = df_filtered["stitle"].apply(
        lambda x: re.findall(r"\b[A-Z][a-z]+ [a-z]+", x)[0] if re.findall(r"\b[A-Z][a-z]+ [a-z]+", x) else "Unknown"
    )

    # Write header to TSV
    df_blast.to_csv(blast_out, sep='\t', index=False)

    # Write filtered results to Excel by genus
    excel_writer = pd.ExcelWriter(blast_excel, engine='xlsxwriter')
    for group_name in sorted(df_filtered["group"].unique()):
        df_group = df_filtered[df_filtered["group"] == group_name]
        df_group.to_excel(excel_writer, sheet_name=group_name[:31], index=False)
    excel_writer.close()


if __name__ == '__main__':
    main()