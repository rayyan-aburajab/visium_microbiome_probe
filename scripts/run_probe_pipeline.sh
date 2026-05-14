#!/bin/bash
DATA_DIR="$WD_SEQ/data"

FASTA_DIR="$DATA_DIR/FASTA_Bacteria_NCBI"
TSV_FILE="$DATA_DIR/accessions.tsv"
CONSENSUS_DIR="$DATA_DIR/FASTA_Bacteria_consensus"
PROBE_DIR="$DATA_DIR/probes/tsv"
PROBE_FASTA_DIR="$DATA_DIR/probes/fasta"
PROBE_EXCEL="$DATA_DIR/probes/bacteria_probes.xlsx"

SCRIPT_CONSENSUS="$(dirname "$0")/1_consensus_species_seq.py"
SCRIPT_DESIGN="$(dirname "$0")/2_design_probes.py"


echo "[Step 1] Generating consensus FASTAs..."
python "$SCRIPT_CONSENSUS" \
    --fasta_dir "$FASTA_DIR" \
    --tsv_file "$TSV_FILE" \
    --output_dir "$CONSENSUS_DIR"

echo "[Step 2] Designing and aggregating probes..."
python "$SCRIPT_DESIGN" \
    --input_dir "$CONSENSUS_DIR/consensus" \
    --output_dir "$PROBE_DIR" \
    --output_excel "$PROBE_EXCEL" \
    --output_fasta_dir "$PROBE_FASTA_DIR"

echo "Done"