#!/bin/bash

# Binary bigWig file merging script (specifically handles 0/1 signals)
# Supports intersection and union operations, ensures output remains 0/1 binary signals
# Usage: ./binary_merge_bigwig.sh <input_dir> <output_dir> <merge_method>

# Check parameters
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <merge_method>"
    echo ""
    echo "merge_method options:"
    echo "  intersect - Intersection: only regions where all files are 1 will be 1"
    echo "  union     - Union: regions where any file is 1 will be 1"
    echo ""
    echo "Examples:"
    echo "  $0 /path/to/bigwig /path/to/output intersect"
    echo "  $0 /path/to/bigwig /path/to/output union"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
MERGE_METHOD=$3
GENOME_SIZE="$(dirname "$0")/mm9.chrom.sizes"

# Validate merge method
case $MERGE_METHOD in
    intersect|union)
        ;;
    *)
        echo "Error: Unsupported merge method '$MERGE_METHOD'"
        echo "Supported methods: intersect, union"
        exit 1
        ;;
esac

# Check required tools
for tool in bigWigToBedGraph bedGraphToBigWig bedtools awk; do
    command -v $tool >/dev/null 2>&1 || { 
        echo "$tool tool installation required"; 
        exit 1; 
    }
done

# Check genome size file
if [ ! -f "$GENOME_SIZE" ]; then
    echo "Genome size file not found: $GENOME_SIZE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/temp"

echo "Starting binary bigWig file merging..."
echo "Merge method: $MERGE_METHOD"
echo ""

# Define ENCODE ID to histone modification mapping
declare -A encode_to_modification

# chage as needed for your specific dataset

# H3k4me3
encode_to_modification["ENCFF218SDV"]="H3k4me3"  # Cortex
encode_to_modification["ENCFF609BZI"]="H3k4me3"  # Intestine
encode_to_modification["ENCFF26XUO"]="H3k4me3"   # Cerebellum

# H3k4me1
encode_to_modification["ENCFF407TGJ"]="H3k4me1"  # Cortex
encode_to_modification["ENCFF523POG"]="H3k4me1"  # Intestine
encode_to_modification["ENCFF934THC"]="H3k4me1"  # Cerebellum

# H3k27ac
encode_to_modification["ENCFF958PXM"]="H3k27ac"  # Cortex
encode_to_modification["ENCFF182RON"]="H3k27ac"  # Intestine
encode_to_modification["ENCFF763RNO"]="H3k27ac"  # Cerebellum

# Define histone modification types
declare -A modifications
modifications["H3k4me3"]="H3k4me3"
modifications["H3k4me1"]="H3k4me1" 
modifications["H3k27ac"]="H3k27ac"

# Function: Extract ENCODE ID from filename
extract_encode_id() {
    local filename=$(basename "$1")
    echo "$filename" | grep -o 'ENCFF[A-Z0-9]*'
}

# Function: Calculate binary union
calculate_binary_union() {
    local mod=$1
    local -a files=("${!2}")
    
    echo "Calculating binary union (regions where any file is 1 will be 1)..."
    
    # Use bedtools unionbedg to get values at all positions
    bedtools unionbedg -i "${files[@]}" > "$OUTPUT_DIR/temp/${mod}_union_raw.bedGraph"
    
    # For union: if any file has signal 1 in the region, result is 1
    awk 'BEGIN{OFS="\t"} {
        has_signal = 0
        for(i=4; i<=NF; i++) {
            if($i > 0) {
                has_signal = 1
                break
            }
        }
        if(has_signal) print $1, $2, $3, 1
    }' "$OUTPUT_DIR/temp/${mod}_union_raw.bedGraph" > "$OUTPUT_DIR/temp/${mod}_final.bedGraph"
}

# Function: Calculate binary intersection
calculate_binary_intersect() {
    local mod=$1
    local -a files=("${!2}")
    
    echo "Calculating binary intersection (only regions where all files are 1 will be 1)..."
    
    # Use bedtools unionbedg to get values at all positions
    bedtools unionbedg -i "${files[@]}" > "$OUTPUT_DIR/temp/${mod}_union_raw.bedGraph"
    
    # For intersection: only if all files have signal 1 in the region, result is 1
    awk -v num_files=${#files[@]} 'BEGIN{OFS="\t"} {
        signal_count = 0
        for(i=4; i<=NF; i++) {
            if($i > 0) {
                signal_count++
            }
        }
        if(signal_count == num_files) print $1, $2, $3, 1
    }' "$OUTPUT_DIR/temp/${mod}_union_raw.bedGraph" > "$OUTPUT_DIR/temp/${mod}_final.bedGraph"
}

# Process files for each modification type
for mod in "${!modifications[@]}"; do
    echo ""
    echo "=== Processing $mod modification ==="
    
    # Find all bigWig files for current modification type
    files=()
    for file in "$INPUT_DIR"/*.bigWig "$INPUT_DIR"/*.bw; do
        if [ -f "$file" ]; then
            encode_id=$(extract_encode_id "$file")
            if [ -n "$encode_id" ] && [ "${encode_to_modification[$encode_id]}" = "$mod" ]; then
                files+=("$file")
            fi
        fi
    done
    
    if [ ${#files[@]} -eq 0 ]; then
        echo "Warning: No bigWig files found for $mod"
        continue
    fi
    
    if [ ${#files[@]} -eq 1 ]; then
        echo "Only 1 $mod file found, copying directly..."
        cp "${files[0]}" "$OUTPUT_DIR/${mod}_mm9.bigWig"
        continue
    fi
    
    echo "Found ${#files[@]} $mod files:"
    for file in "${files[@]}"; do
        encode_id=$(extract_encode_id "$file")
        echo "  - $(basename "$file") (ENCODE ID: $encode_id)"
    done
    
    # Validate if input files are binary signals
    echo "Validating input file signal types..."
    for file in "${files[@]}"; do
        # Convert to bedGraph and check values in first 1000 lines
        bigWigToBedGraph "$file" "$OUTPUT_DIR/temp/check_$(basename "$file").bedGraph"
        non_binary=$(awk '{if($4 != 0 && $4 != 1) count++} END{print count+0}' "$OUTPUT_DIR/temp/check_$(basename "$file").bedGraph" | head -1000)
        if [ "$non_binary" -gt 0 ]; then
            echo "Warning: $(basename "$file") contains non-0/1 values, will be converted to binary"
            # Convert non-zero values to 1
            awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ($4 > 0 ? 1 : 0)}' "$OUTPUT_DIR/temp/check_$(basename "$file").bedGraph" > "$OUTPUT_DIR/temp/check_$(basename "$file")_binary.bedGraph"
            mv "$OUTPUT_DIR/temp/check_$(basename "$file")_binary.bedGraph" "$OUTPUT_DIR/temp/check_$(basename "$file").bedGraph"
        fi
        rm -f "$OUTPUT_DIR/temp/check_$(basename "$file").bedGraph"
    done
    
    # Convert all related bigWig files to bedGraph
    bedgraph_files=()
    for i in "${!files[@]}"; do
        file="${files[$i]}"
        filename=$(basename "$file" | sed 's/\.[^.]*$//')
        echo "Converting $(basename "$file") to bedGraph..."
        bigWigToBedGraph "$file" "$OUTPUT_DIR/temp/${mod}_${i}.bedGraph"
        
        # Ensure binary format conversion and sort
        awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ($4 > 0 ? 1 : 0)}' "$OUTPUT_DIR/temp/${mod}_${i}.bedGraph" | \
        sort -k1,1 -k2,2n > "$OUTPUT_DIR/temp/${mod}_${i}_sorted.bedGraph"
        
        bedgraph_files+=("$OUTPUT_DIR/temp/${mod}_${i}_sorted.bedGraph")
    done
    
    # Process according to merge method
    case $MERGE_METHOD in
        "union")
            calculate_binary_union "$mod" bedgraph_files[@]
            ;;
        "intersect")
            calculate_binary_intersect "$mod" bedgraph_files[@]
            ;;
    esac
    
    # Validate final result is binary
    echo "Validating final result..."
    non_binary_final=$(awk '{if($4 != 0 && $4 != 1) print "Non-binary value:", $0}' "$OUTPUT_DIR/temp/${mod}_final.bedGraph")
    if [ -n "$non_binary_final" ]; then
        echo "Error: Final result contains non-0/1 values:"
        echo "$non_binary_final" | head -5
        exit 1
    fi
    
    echo "Generating genome-wide coverage file..."
    # Generate genome-wide coverage (including empty regions as 0)
    bedtools genomecov -i "$OUTPUT_DIR/temp/${mod}_final.bedGraph" \
        -g "$GENOME_SIZE" \
        -bga > "$OUTPUT_DIR/temp/${mod}_full.bedGraph"
    
    echo "Converting to bigWig format..."
    # Convert back to bigWig format
    bedGraphToBigWig "$OUTPUT_DIR/temp/${mod}_full.bedGraph" \
        "$GENOME_SIZE" \
        "$OUTPUT_DIR/${mod}_mm9.bigWig"
    
    # Output statistics
    total_regions=$(wc -l < "$OUTPUT_DIR/temp/${mod}_final.bedGraph")
    total_positions=$(wc -l < "$OUTPUT_DIR/temp/${mod}_union_raw.bedGraph")
    
    # Calculate coverage statistics
    total_bp=$(awk '{sum += ($3-$2)} END{print sum+0}' "$OUTPUT_DIR/temp/${mod}_final.bedGraph")
    
    echo "Statistics:"
    echo "  - Analyzed genome positions: $total_positions"
    echo "  - Regions with signal: $total_regions" 
    echo "  - Total covered base pairs: $total_bp bp"
    
    if [ "$MERGE_METHOD" = "union" ]; then
        echo "  - Union strategy: regions with signal in any sample are retained"
    else
        echo "  - Intersection strategy: only regions with signal in all samples are retained"
    fi
    
    echo "Completed $mod processing, output file: ${mod}_mm9.bigWig"
done

# Clean up temporary files
echo ""
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/temp"

echo ""
echo "==============================================="
echo "Binary bigWig merging completed!"
echo "Output files located at: $OUTPUT_DIR"
echo "Merge method used: $MERGE_METHOD"
echo ""
echo "Generated files (all are 0/1 binary signals):"
for mod in "${!modifications[@]}"; do
    if [ -f "$OUTPUT_DIR/${mod}_mm9.bigWig" ]; then
        size=$(du -h "$OUTPUT_DIR/${mod}_mm9.bigWig" | cut -f1)
        echo "  - ${mod}_mm9.bigWig ($size)"
    fi
done
echo ""
echo "Note: All output files are binary signals (0 or 1)"
echo "==============================================="