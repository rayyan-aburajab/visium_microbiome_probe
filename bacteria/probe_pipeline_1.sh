#!/bin/bash

WD_SEQ="/Users/raburajab/Documents/microbiome/sequences"
WD_SCRIPTS="/Users/raburajab/Documents/microbiome/visium_microbiome/bacteria"

FASTA_DIR="$WD_SEQ/FASTA_Bacteria_NCBI"
TSV_FILE="$WD_SEQ/bacteria_accessions_all.tsv"
CONSENSUS_DIR="$WD_SEQ/FASTA_Bacteria_consensus"
ALL_PROBE_DIR="$WD_SEQ/probe_outputs/bacteria_probes"

SCRIPT_CONSENSUS="$WD_SCRIPTS/1_consensus_species_seq.py"
SCRIPT_DESIGN="$WD_SCRIPTS/2_design_probes.py"

mkdir -p "$CONSENSUS_DIR"/combined "$CONSENSUS_DIR"/aligned "$CONSENSUS_DIR"/consensus
mkdir -p "$ALL_PROBE_DIR" "$QC_PROBE_DIR"

echo "[Step 1] Generating consensus FASTAs..."
python "$SCRIPT_CONSENSUS" \
  --fasta_dir "$FASTA_DIR" \
  --tsv_file "$TSV_FILE" \
  --output_dir "$CONSENSUS_DIR"

echo "[Step 2] Designing probes from consensus FASTAs..."
for fasta in "$CONSENSUS_DIR/consensus"/*_consensus.fa; do
    [ -e "$fasta" ] || continue  # skip if no files
    species_tag=$(basename "$fasta" _consensus.fa)
    echo "Designing for $species_tag..."
    python "$SCRIPT_DESIGN" \
      "$fasta" \
      --species_name "$species_tag" \
      --output_dir "$ALL_PROBE_DIR"
done

echo "✅ Steps 1 and 2 complete"
