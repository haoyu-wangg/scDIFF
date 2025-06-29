#!/bin/bash

# Script to convert mm10 (mouse) bed files to hg19 (human) using liftOver
# Usage: ./mm10_to_hg19_liftover.sh -i input.bed -c chain_file.chain.gz -o output_dir -m 0.8

# Exit on error
set -e

# Default values
INPUT_BED=""
CHAIN_FILE=""
OUTPUT_DIR="."
MIN_MATCH=0.95  # Default minimum match value (95%)

# Function to display usage
usage() {
    echo "Usage: $0 -i <input.bed> [-c <chain_file>] [-o <output_directory>] [-m <min_match>]"
    echo "  -i  Input BED file in mm10 coordinates"
    echo "  -c  Chain file for mm10 to hg19 conversion (optional, will download if not provided)"
    echo "  -o  Output directory (optional, defaults to current directory)"
    echo "  -m  Minimum match ratio (optional, defaults to 0.95, range 0.0-1.0)"
    echo "  -h  Display this help message"
    exit 1
}

# Parse command line arguments
while getopts "i:c:o:m:h" opt; do
    case ${opt} in
        i )
            INPUT_BED=$OPTARG
            ;;
        c )
            CHAIN_FILE=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        m )
            MIN_MATCH=$OPTARG
            ;;
        h )
            usage
            ;;
        \? )
            usage
            ;;
    esac
done

# Check if required parameters are provided
if [ -z "$INPUT_BED" ]; then
    echo "Error: Input BED file is required"
    usage
fi

# Validate minMatch parameter
if (( $(echo "$MIN_MATCH < 0.0" | bc -l) )) || (( $(echo "$MIN_MATCH > 1.0" | bc -l) )); then
    echo "Error: Minimum match value must be between 0.0 and 1.0"
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_BED" ]; then
    echo "Error: Input file $INPUT_BED does not exist"
    exit 1
fi

# Check if output directory exists, create if not
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Get the base name of the input file
INPUT_BASE=$(basename "$INPUT_BED")
INPUT_NAME="${INPUT_BASE%.bed}"

# Set output file paths
OUTPUT_BED="$OUTPUT_DIR/${INPUT_NAME}_hg19.bed"
UNMAPPED="$OUTPUT_DIR/${INPUT_NAME}_unmapped.bed"

# Check if liftOver is installed
if ! command -v liftOver &> /dev/null; then
    echo "Error: liftOver is not installed or not in PATH"
    echo "Please install liftOver from UCSC: http://hgdownload.soe.ucsc.edu/admin/exe/"
    exit 1
fi

# Check if chain file is provided, download if not
if [ -z "$CHAIN_FILE" ]; then
    CHAIN_FILE="$OUTPUT_DIR/mm10ToHg19.over.chain.gz"
    CHAIN_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg19.over.chain.gz"
    
    if [ ! -f "$CHAIN_FILE" ]; then
        echo "Downloading chain file from UCSC to $CHAIN_FILE..."
        wget -P "$OUTPUT_DIR" "$CHAIN_URL" || curl -o "$CHAIN_FILE" "$CHAIN_URL"
    else
        echo "Using existing chain file: $CHAIN_FILE"
    fi
else
    # Check if provided chain file exists
    if [ ! -f "$CHAIN_FILE" ]; then
        echo "Error: Chain file $CHAIN_FILE does not exist"
        exit 1
    fi
    echo "Using provided chain file: $CHAIN_FILE"
fi

# Run liftOver with minMatch parameter
echo "Converting $INPUT_BED from mm10 to hg19 with minimum match threshold of $MIN_MATCH..."
liftOver -minMatch=$MIN_MATCH "$INPUT_BED" "$CHAIN_FILE" "$OUTPUT_BED" "$UNMAPPED"

# Check results
TOTAL=$(wc -l < "$INPUT_BED")
MAPPED=$(wc -l < "$OUTPUT_BED")
FAILED=$(wc -l < "$UNMAPPED")

echo "Conversion complete!"
echo "Total records: $TOTAL"
echo "Successfully mapped to hg19: $MAPPED"
echo "Failed to map: $FAILED"
echo "Output file: $OUTPUT_BED"
echo "Unmapped regions: $UNMAPPED"
echo "Minimum match threshold used: $MIN_MATCH"