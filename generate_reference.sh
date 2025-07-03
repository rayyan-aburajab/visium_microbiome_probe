#!/bin/bash

BASE_DIR="/Users/raburajab/Documents/microbiome/sequences/Reference_custom"

INPUT_FILE="$BASE_DIR/REF_INPUT.csv"
OUTPUT_PREFIX="$BASE_DIR/MicrobiomeV2-2_Ref"

python3 generate_reference_v2.py --input_file "$INPUT_FILE" --output_prefix "$OUTPUT_PREFIX"
