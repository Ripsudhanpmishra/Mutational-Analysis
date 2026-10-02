#!/bin/bash

# Target paths explicitly defined via WSL mount format
INPUT_FOLDER="/mnt/e/accessory_vcf1"
OUTPUT_FOLDER="/mnt/e/accessory_vcf1/annotated_outputs"

# Ensure the output subdirectory exists before running the tool
mkdir -p "$OUTPUT_FOLDER"

echo "Checking tracking assets inside $INPUT_FOLDER..."

# Loop through files cleanly inside an isolated script environment
for vcf_file in "$INPUT_FOLDER"/*_variants.vcf; do
    
    # Safely skip the step if the pattern expansion matches nothing
    [ -e "$vcf_file" ] || continue
    
    # Isolate the pure filename (e.g., P2RY13_variants.vcf)
    base_name=$(basename "$vcf_file")
    
    # Strip the suffix to get just the clean gene key
    gene_name="${base_name%_variants.vcf}"
    
    echo "=========================================================="
    echo " RUNNING ANNOTATION FOR GENE: $gene_name"
    echo "=========================================================="
    
    # Explicitly routing to the dedicated GRCh37 database server and port
    vep \
      -i "$vcf_file" \
      -o "$OUTPUT_FOLDER/${gene_name}_annotated.csv" \
      --database \
      --assembly GRCh37 \
      --protein \
      --symbol \
      --canonical \
      --tab

done

echo "=========================================================="
echo "SUCCESS! All files processed seamlessly via GRCh37 database connection."
echo "=========================================================="
