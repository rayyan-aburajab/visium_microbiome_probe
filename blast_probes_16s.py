#!/usr/bin/env python3

import pandas as pd
import subprocess
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import click

@click.command()
@click.option('--csv', required=True, help='Input CSV with probe_id, hyb_region_lhs, hyb_region_rhs')
@click.option('--fasta', required=True, help='Output FASTA file for BLAST input')
@click.option('--blastdb', required=True, help='Path to 16S BLAST database (no extension)')
@click.option('--blastout', required=True, help='Output file for BLAST results')
@click.option('--output', required=True, help='Output Excel file (species vs probe group)')
def main(csv, fasta, blastdb, blastout, output):

    # === Step 1: Read CSV and build full probe sequence ===
    df = pd.read_csv(csv)
    df["sequence"] = df["hyb_region_lhs"] + df["hyb_region_rhs"]

    # Write FASTA
    records = [
        SeqRecord(Seq(row["sequence"]), id=row["probe_id"], description="")
        for _, row in df.iterrows()
    ]
    SeqIO.write(records, fasta, "fasta")

    # === Step 2: Run BLAST ===
    cmd = [
        "blastn",
        "-query", fasta,
        "-db", blastdb,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
        "-out", blastout,
        "-evalue", "1e-5",
        "-max_target_seqs", "100"
    ]
    print("🔍 Running BLAST...")
    subprocess.run(cmd, check=True)

    # === Step 3: Parse BLAST results and filter ===
    print("Parsing BLAST results...")
    probe_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}
    df_blast = pd.read_csv(blastout, sep="\t", header=None)
    df_blast.columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "stitle"
    ]
    df_blast["qcov"] = df_blast.apply(lambda row: row["length"] / probe_lengths[row["qseqid"]], axis=1)

    df_filtered = df_blast[(df_blast["pident"] >= 90.0) & (df_blast["qcov"] >= 0.9)]

    # Extract group (genus/family) from probe_id (e.g., "Alistipes_probe1" -> "Alistipes")
    df_filtered["group"] = df_filtered["qseqid"].apply(lambda x: x.split("_")[0])

    # Extract species name from stitle
    df_filtered["species"] = df_filtered["stitle"].apply(
        lambda x: re.findall(r"\b[A-Z][a-z]+ [a-z]+", x)[0] if re.findall(r"\b[A-Z][a-z]+ [a-z]+", x) else "Unknown"
    )

    # === Step 4: Build species-by-group matrix ===
    print("Building species matrix...")
    species_to_groups = defaultdict(set)
    for _, row in df_filtered.iterrows():
        species_to_groups[row["species"]].add(row["group"])

    all_species = sorted(species_to_groups.keys())
    all_groups = sorted(set(g for gs in species_to_groups.values() for g in gs))

    matrix = pd.DataFrame(0, index=all_species, columns=all_groups)
    for sp, groups in species_to_groups.items():
        for g in groups:
            matrix.loc[sp, g] = 1

    matrix.to_excel(output)
    print(f"Done! Matrix saved to: {output}")

if __name__ == "__main__":
    main()
