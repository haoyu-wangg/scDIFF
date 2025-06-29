#!/bin/bash

# Check parameters
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    echo "Example: $0 /path/to/input /path/to/output"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
GENOME_SIZE="$(dirname "$0")/mm9.chrom.sizes" # Use chromosome size file from script directory
CHAIN_FILE="$(dirname "$0")/mm10ToMm9.over.chain" # Use chain file from script directory

# Check required tools
command -v liftOver >/dev/null 2>&1 || { echo "liftOver tool installation required"; exit 1; }
command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "bedGraphToBigWig tool installation required"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools installation required"; exit 1; }

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
    # MACS2 BED format typically contains: chr start end name score strand...
    # We extract chr, start, end, score columns
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' "$OUTPUT_DIR/temp/$filename.bed" > "$OUTPUT_DIR/temp/${filename}.bedGraph"

    # 3. Use liftOver for coordinate conversion
    liftOver "$OUTPUT_DIR/temp/${filename}.bedGraph" \
        "$CHAIN_FILE" \
        "$OUTPUT_DIR/temp/${filename}.mm9.bedGraph" \
        "$OUTPUT_DIR/temp/${filename}.unmapped"

    # 4. Check and sort converted file
    sort -k1,1 -k2,2n "$OUTPUT_DIR/temp/${filename}.mm9.bedGraph" > "$OUTPUT_DIR/temp/${filename}.mm9.sorted.bedGraph"

    # 4.1 Handle overlapping regions (use bedtools merge and calculate maximum values)
    bedtools merge -i "$OUTPUT_DIR/temp/${filename}.mm9.sorted.bedGraph" \
        -c 4 \
        -o max \
        > "$OUTPUT_DIR/temp/${filename}.mm9.merged.bedGraph"

    # 4.2 Create genome file
    cut -f1,2 "$GENOME_SIZE" > "$OUTPUT_DIR/temp/${filename}.genome"

    # 4.3 Generate genome-wide coverage file
    bedtools genomecov -i "$OUTPUT_DIR/temp/${filename}.mm9.merged.bedGraph" \
        -g "$OUTPUT_DIR/temp/${filename}.genome" \
        -bga > "$OUTPUT_DIR/temp/${filename}.mm9.full.bedGraph"
        #-bga > "$OUTPUT_DIR/temp/${filename}.mm9.full.bedGraph"

    # 5. Convert to bigwig format
    bedGraphToBigWig "$OUTPUT_DIR/temp/${filename}.mm9.full.bedGraph" \
        "$GENOME_SIZE" \
        "$OUTPUT_DIR/${filename}.bigWig"

    # 6. Clean up temporary files
    rm "$OUTPUT_DIR/temp/$filename.bed" \
        "$OUTPUT_DIR/temp/${filename}.bedGraph" \
        "$OUTPUT_DIR/temp/${filename}.mm9.bedGraph" \
        "$OUTPUT_DIR/temp/${filename}.unmapped" \
        "$OUTPUT_DIR/temp/${filename}.mm9.sorted.bedGraph" \
        "$OUTPUT_DIR/temp/${filename}.mm9.merged.bedGraph" \
        "$OUTPUT_DIR/temp/${filename}.genome" \
        "$OUTPUT_DIR/temp/${filename}.mm9.full.bedGraph"

done

# Remove temporary directory
rm -r "$OUTPUT_DIR/temp"
echo "Conversion completed! Output files located at: $OUTPUT_DIR"