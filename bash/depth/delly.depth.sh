#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="5 10 15 20 25 30 45 60"  # Different n values

delly_path="${BASE_DIR}/result/WGS-SVs/bamfile/delly"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"

# ========== Initialization ==========
mkdir -p "$delly_path" "$output_base"

# Activate Delly environment
source activate delly

# ========== Main processing function ==========
process_sample() {
    local bamfile="$1"
    local n="$2"
    
    # Create sample-specific directory
    local input_bam="${input_base}/bam_p${n}/${bamfile}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    mkdir -p "${delly_path}/${bamfile}.p${n}" "${output_dir}"

    # Check if input file exists
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file not found: $input_bam"
        return 1
    fi

    echo "Start processing sample: ${bamfile}.p${n}"
    
    # Step 1: Delly SV detection
    echo "Running Delly detection..."
    delly call \
        -g "$reference_genome" \
        -o "${delly_path}/${bamfile}.p${n}/${bamfile}.bcf" \
        "$input_bam"

    # Step 2: Filter PASS variants
    echo "Filtering PASS variants..."
    bcftools view \
        -f PASS \
        -r "$chromosomes" \
        "${delly_path}/${bamfile}.p${n}/${bamfile}.bcf" \
        > "${output_dir}/${bamfile}.delly.vcf"

    # Step 3: Compress and index
    echo "Compressing and indexing VCF..."
    bgzip -f "${output_dir}/${bamfile}.delly.vcf"
    tabix "${output_dir}/${bamfile}.delly.vcf.gz"

    echo "Sample ${bamfile}.p${n} completed"
    return 0
}

# ========== Main Process ==========

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi

# Read sample list
readarray -t samples < "$bamfile_list"

# Process each sample for each n value
for sample in "${samples[@]}"; do
    # Remove possible whitespace
    sample=$(echo "$sample" | xargs)
    
    # Skip empty lines
    if [[ -z "$sample" ]]; then
        continue
    fi
    
    # Process sample for each n value
    for n in $n_values; do
        echo "Processing sample: $sample using n=$n"
        if process_sample "$sample" "$n"; then
            echo "Successfully processed: $sample.p${n}"
        else
            echo "Processing failed: $sample.p${n}" >&2
        fi
    done
done

echo "All samples processed"