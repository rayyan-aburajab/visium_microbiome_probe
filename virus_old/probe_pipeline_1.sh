#!/bin/bash

# === Configuration ===
WD_SEQ="/Users/raburajab/Documents/microbiome/sequences"
WD_SCRIPTS="/Users/raburajab/Documents/microbiome/visium_microbiome/virus"

FASTA_DIR="$WD_SEQ/FASTA_Virus_NCBI"
COMBINED_DIR="$WD_SEQ/FASTA_Virus_combined"
ALL_PROBE_DIR="$WD_SEQ/probe_outputs/virus_probes"

SCRIPT_COMBINE="$WD_SCRIPTS/1_combine_virus_seq.py"
SCRIPT_DESIGN="$WD_SCRIPTS/2_design_probes.py"

# === Create output directories ===
mkdir -p "$COMBINED_DIR" "$ALL_PROBE_DIR"

# === Step 1: Combine viral regions into a single FASTA per virus ===
echo "🧬 [Step 1] Combining viral FASTAs..."
python "$SCRIPT_COMBINE" \
  --fasta_dir "$FASTA_DIR" \
  --output_dir "$COMBINED_DIR"

# === Step 2: Design probes from combined FASTAs ===
echo "🧪 [Step 2] Designing probes from combined FASTAs..."
for fasta in "$COMBINED_DIR"/*_combined.fa; do
    [ -e "$fasta" ] || continue  # skip if no files
    virus_tag=$(basename "$fasta" _combined.fa)
    echo "  📐 Designing for $virus_tag..."
    python "$SCRIPT_DESIGN" \
      "$fasta" \
      --species_name "$virus_tag" \
      --output_dir "$ALL_PROBE_DIR"
done

# === Completion ===
echo "✅ Full virus pipeline complete: Combine ➜ Probe Design"
