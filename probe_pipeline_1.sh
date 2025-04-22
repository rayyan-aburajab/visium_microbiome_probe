#!/bin/bash

# === Configuration ===
WD_SEQ="/Users/raburajab/Documents/microbiome/sequences"
WD_SCRIPTS="/Users/raburajab/Documents/microbiome/visium-probe-microbiome"

FASTA_DIR="$WD_SEQ/FASTA_NCBI"
TSV_FILE="$WD_SEQ/species_accessions_all.tsv"
CONSENSUS_DIR="$WD_SEQ/species_consensus_FASTAs"
ALL_PROBE_DIR="$WD_SEQ/probe_outputs/all_probes"

SCRIPT_CONSENSUS="$WD_SCRIPTS/1_consensus_species_seq.py"
SCRIPT_DESIGN="$WD_SCRIPTS/2_design_probes.py"

# === Create output directories ===
mkdir -p "$CONSENSUS_DIR"/combined "$CONSENSUS_DIR"/aligned "$CONSENSUS_DIR"/consensus
mkdir -p "$ALL_PROBE_DIR" "$QC_PROBE_DIR"

# === Step 1: Build consensus sequences per species ===
echo "🧬 [Step 1] Generating consensus FASTAs..."
python "$SCRIPT_CONSENSUS" \
  --fasta_dir "$FASTA_DIR" \
  --tsv_file "$TSV_FILE" \
  --output_dir "$CONSENSUS_DIR"

# === Step 2: Design probes from consensus FASTAs ===
echo "🧪 [Step 2] Designing probes from consensus FASTAs..."
for fasta in "$CONSENSUS_DIR/consensus"/*_consensus.fa; do
    [ -e "$fasta" ] || continue  # skip if no files
    species_tag=$(basename "$fasta" _consensus.fa)
    echo "  📐 Designing for $species_tag..."
    python "$SCRIPT_DESIGN" \
      "$fasta" \
      --species_name "$species_tag" \
      --output_dir "$ALL_PROBE_DIR"
done

# === Completion ===
echo "✅ Full pipeline complete: Consensus ➜ Probe Design ➜ QC"
