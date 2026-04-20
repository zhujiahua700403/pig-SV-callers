#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="5 10 15 20 25 30 45 60"  # Different n values

survindel_dir="${BASE_DIR}/result/WGS-SVs/bamfile/SurVIndel2"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"

# ========== Initialization ==========
mkdir -p "$survindel_dir" "$output_base"

# Activate SurVIndel2 environment
source activate SurVIndel2

# ========== Sample processing function ==========
process_sample() {
    local sample="$1"
    local n="$2"
    
    # Create output directories
    local input_bam="${input_base}/bam_p${n}/${sample}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    local sample_dir="${survindel_dir}/${sample}.p${n}"
    mkdir -p "$sample_dir" "$output_dir"

    # Check if input file exists
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file not found: $input_bam"
        return 1
    fi

    echo "Start processing sample: ${sample}.p${n}"

    # Step 1: Run SurVIndel2
    echo "Running SurVIndel2 detection..."
    survindel2.py \
        "$input_bam" \
        "$sample_dir" \
        "$reference_genome" \
        --threads 1 \
        --min_sv_size 50 \
        --samplename "${sample}.SurVIndel2"

    # Check output file
    if [ ! -f "${sample_dir}/out.pass.vcf.gz" ]; then
        echo "Error: SurVIndel2 output file not found"
        return 1
    fi

    # Step 2: Filter chromosomes
    echo "Filtering chromosomes..."
    bcftools view \
        -r "$chromosomes" \
        "${sample_dir}/out.pass.vcf.gz" \
        -Oz -o "${output_dir}/${sample}.SurVIndel2.vcf.gz"

    if [ $? -ne 0 ]; then
        echo "Error: bcftools view failed"
        return 1
    fi

    # Step 3: Create index
    tabix "${output_dir}/${sample}.SurVIndel2.vcf.gz"
    if [ $? -ne 0 ]; then
        echo "Error: tabix failed"
        return 1
    fi

    echo "Sample ${sample}.p${n} completed"
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