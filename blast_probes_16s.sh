#!/bin/bash

# === Base directory ===
BASE_DIR="/Users/raburajab/Documents/microbiome/sequences"
DOC_DIR="$BASE_DIR/probe_output_docs"
DB_DIR="$BASE_DIR/16S_refseq"

# === File paths ===
PROBE_CSV="$DOC_DIR/ProbeList_bacteria.csv"
FASTA_FILE="$DOC_DIR/ProbeList_bacteria.fasta"
BLAST_DB="$DB_DIR/16S_ribosomal_RNA"
BLAST_OUT="$DOC_DIR/blast_results.tsv"
OUTPUT_MATRIX="$DOC_DIR/species_probe_matrix.xlsx"

# === Run the Python pipeline ===
python3 blast_probes_stringent_16s.py \
  --csv "$PROBE_CSV" \
  --fasta "$FASTA_FILE" \
  --blastdb "$BLAST_DB" \
  --blastout "$BLAST_OUT" \
  --output "$OUTPUT_MATRIX"

