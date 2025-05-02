#!/bin/bash

# === Step 3: Combine hybrid-only probes and generate FASTA ===
echo "🔬 [Step 3] Combining probe hybrid regions and writing FASTA..."

# Paths
ALL_PROBE_DIR="/Users/raburajab/Documents/microbiome/sequences/probe_outputs/bacteria_probes"

SCRIPT_COMBINE="/Users/raburajab/Documents/microbiome/visium-probe-microbiome/3_combine_probes.py"

COMBINED_TSV="/Users/raburajab/Documents/microbiome/sequences/probe_outputs/bacteria_probes.tsv"
COMBINED_FASTA="/Users/raburajab/Documents/microbiome/sequences/probe_outputs/bacteria_probes.fasta"

# Run combination script
python "$SCRIPT_COMBINE" \
  --input_dir "$ALL_PROBE_DIR" \
  --output_file "$COMBINED_TSV" \
  --output_fasta "$COMBINED_FASTA"

echo "✅ Step 3 complete: Hybrid probe file and FASTA created"
