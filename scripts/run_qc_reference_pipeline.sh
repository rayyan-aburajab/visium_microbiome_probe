#!/bin/bash
DATA_DIR="$WD_SEQ/data"

SELECTION="$WD_SEQ/inputs/selected_probes.tsv"
PROBE_DIR="$DATA_DIR/probes/tsv"
QC_DIR="$DATA_DIR/qc"
REF_DIR="$DATA_DIR/reference"

SCRIPT_QC="$(dirname "$0")/4_select_and_qc.py"
SCRIPT_REF="$(dirname "$0")/5_generate_reference.py"

echo "[Step 1] Running probe selection and QC..."
python "$SCRIPT_QC" \
    --selection "$SELECTION" \
    --probe_dir "$PROBE_DIR" \
    --output_dir "$QC_DIR"

echo "[Step 2] Generating reference files..."
python "$SCRIPT_REF" \
    --input "$QC_DIR/probe_qc.tsv" \
    --output_dir "$REF_DIR" \
    --output_prefix "bacteria_probes"

echo "Done"