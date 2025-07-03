#!/bin/bash

WD_SEQ="/Users/raburajab/Documents/microbiome/sequences"
WD_SCRIPTS="/Users/raburajab/Documents/microbiome/visium_microbiome/bacteria"

ALL_PROBE_DIR="$WD_SEQ/probe_outputs/bacteria_probes"
OUT_DIR="$WD_SEQ/probe_outputs"
FASTA_DIR="$OUT_DIR/bacteria_fasta"
SCRIPT_COMBINE="$WD_SCRIPTS/3_combine_probes.py"

COMBINED_EXCEL="$OUT_DIR/bacteria_probes.xlsx"

echo "[Step 3] Writing FASTA + Excel..."

mkdir -p "$FASTA_DIR"

python "$SCRIPT_COMBINE" \
  --input_dir "$ALL_PROBE_DIR" \
  --output_excel "$COMBINED_EXCEL" \
  --output_fasta_dir "$FASTA_DIR"

echo "✅ Step 3 complete"
