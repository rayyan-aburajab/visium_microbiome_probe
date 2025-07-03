#!/bin/bash

INPUT_DIR="/Users/raburajab/Documents/microbiome/sequences/probe_output_docs"
QC_SCRIPT="evaluate_QC.py"

for csv in "$INPUT_DIR"/*.csv; do
    filename=$(basename "$csv")
    
    # Skip files that already end in _QC.csv
    if [[ "$filename" == *_QC.csv ]]; then
        echo "⚠️ Skipping already processed file: $filename"
        continue
    fi

    base="${filename%.csv}"
    output_csv="${INPUT_DIR}/${base}_QC.csv"

    echo "🔬 Running QC on: $csv"
    python "$QC_SCRIPT" \
        --input_csv "$csv" \
        --output_csv "$output_csv"
done

echo "✅ All QC files saved with _QC.csv in: $INPUT_DIR"
