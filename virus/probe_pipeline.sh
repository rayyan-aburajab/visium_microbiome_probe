#!/bin/bash

# === Configuration ===
WD_SEQ="/Users/raburajab/Documents/microbiome/sequences"
WD_SCRIPTS="/Users/raburajab/Documents/microbiome/visium_microbiome/virus"

FASTA_DIR="$WD_SEQ/FASTA_Virus_CDS"
ALL_PROBE_DIR="$WD_SEQ/probe_outputs/virus_probes"
COMBINED_EXCEL="$WD_SEQ/probe_outputs/virus_all_combined.xlsx"
COMBINED_FASTA="$WD_SEQ/probe_outputs/virus_all_probes.fa"

SCRIPT_DESIGN="$WD_SCRIPTS/1_design_probes.py"
SCRIPT_COMBINE="$WD_SCRIPTS/2_combine_probes.py"

# === Create output directory ===
mkdir -p "$ALL_PROBE_DIR"

# === Step 1: Design probes from viral FASTAs ===
echo "🧪 [Step 1] Designing probes from viral FASTAs..."
for fasta in "$FASTA_DIR"/*.fa; do
    [ -e "$fasta" ] || continue  # skip if no files
    virus_tag=$(basename "$fasta" .fa)
    echo "  📐 Designing for $virus_tag..."
    python "$SCRIPT_DESIGN" \
      "$fasta" \
      --species_name "$virus_tag" \
      --output_dir "$ALL_PROBE_DIR"
done

# === Step 2: Combine and format probes ===
echo "📦 [Step 2] Combining probe files into TSV and FASTA..."
python "$SCRIPT_COMBINE" \
  --input_dir "$ALL_PROBE_DIR" \
  --output_excel "$COMBINED_EXCEL" \
  --output_fasta "$COMBINED_FASTA"

# === Completion ===
echo "✅ Full virus pipeline complete: FASTAs ➜ Probes ➜ Combined Output"