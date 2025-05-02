#!/bin/bash

# === Step 3: Combine probes and generate FASTA ===
echo "🔬 [Step 3] Combining probe regions and writing FASTA..."

# Paths
ALL_PROBE_DIR="/Users/raburajab/Documents/microbiome/sequences/probe_outputs/virus_probes"

SCRIPT_COMBINE="/Users/raburajab/Documents/microbiome/visium_microbiome/virus/3_combine_probes.py"

COMBINED_TSV="/Users/raburajab/Documents/microbiome/sequences/probe_outputs/virus_probes.tsv"
COMBINED_FASTA="/Users/raburajab/Documents/microbiome/sequences/probe_outputs/virus_probes.fasta"

# Run combination script
python "$SCRIPT_COMBINE" \
  --input_dir "$ALL_PROBE_DIR" \
  --output_file "$COMBINED_TSV" \
  --output_fasta "$COMBINED_FASTA"

echo "✅ Step 3 complete: probe file and FASTA created"
