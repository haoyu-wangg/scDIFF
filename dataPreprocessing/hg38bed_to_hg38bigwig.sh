#!/bin/bash

# Check parameters
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    echo "Example: $0 /path/to/input /path/to/output"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
GENOME_SIZE="$(dirname "$0")/hg38.chrom.sizes"  # Use hg38 chromosome size file in script directory

# Check required tools
command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "bedGraphToBigWig tool is required"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools tool is required"; exit 1; }

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/temp"

# Process each bed.gz file
for file in "$INPUT_DIR"/*.bed.gz; do
    if [ ! -f "$file" ]; then
        echo "No .bed.gz files found"
        exit 1
    fi

    # Get filename (without path and extension)
    filename=$(basename "$file" .bed.gz)
    echo "Processing $filename..."

    # 1. Decompress file
    gunzip -c "$file" > "$OUTPUT_DIR/temp/$filename.bed"

    # 2. Extract required columns and convert to bedGraph format
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' "$OUTPUT_DIR/temp/$filename.bed" > "$OUTPUT_DIR/temp/${filename}.bedGraph"

    # 3. Check and sort bedGraph file
    sort -k1,1 -k2,2n "$OUTPUT_DIR/temp/${filename}.bedGraph" > "$OUTPUT_DIR/temp/${filename}.sorted.bedGraph"

    # 4. Handle overlapping regions (use bedtools merge and calculate maximum value)
    bedtools merge -i "$OUTPUT_DIR/temp/${filename}.sorted.bedGraph" -c 4 -o max > "$OUTPUT_DIR/temp/${filename}.merged.bedGraph"

    # 5. Create genome file
    cut -f1,2 "$GENOME_SIZE" > "$OUTPUT_DIR/temp/${filename}.genome"

    # 6. Generate genome-wide coverage file
    bedtools genomecov -i "$OUTPUT_DIR/temp/${filename}.merged.bedGraph" -g "$OUTPUT_DIR/temp/${filename}.genome" -bga > "$OUTPUT_DIR/temp/${filename}.full.bedGraph"

    # 7. Convert to bigWig format
    bedGraphToBigWig "$OUTPUT_DIR/temp/${filename}.full.bedGraph" "$GENOME_SIZE" "$OUTPUT_DIR/${filename}.bigWig"

    # 8. Clean up temporary files
    rm "$OUTPUT_DIR/temp/$filename.bed" \
       "$OUTPUT_DIR/temp/${filename}.bedGraph" \
       "$OUTPUT_DIR/temp/${filename}.sorted.bedGraph" \
       "$OUTPUT_DIR/temp/${filename}.merged.bedGraph" \
       "$OUTPUT_DIR/temp/${filename}.genome" \
       "$OUTPUT_DIR/temp/${filename}.full.bedGraph"

done

# Remove temporary directory
rm -r "$OUTPUT_DIR/temp"
echo "Conversion completed! Output files are located in: $OUTPUT_DIR"