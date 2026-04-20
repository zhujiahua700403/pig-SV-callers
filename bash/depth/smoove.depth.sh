#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="5 10 15 20 25 30 45 60"  # Different n values

smoove_dir="${BASE_DIR}/result/WGS-SVs/bamfile/smoove"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
software_path="${BASE_DIR}/software"

# ========== Initialization ==========
mkdir -p "$smoove_dir" "$output_base" "$smoove_dir/results-genotyped"

# Activate smoove environment
source activate smoove

# ========== Sample processing function ==========
process_sample() {
    local sample="$1"
    local n="$2"
    
    # Create output directories
    local input_bam="${input_base}/bam_p${n}/${sample}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    local sample_dir="${smoove_dir}/${sample}.p${n}"
    mkdir -p "$sample_dir" "$output_dir"

    # Check if input file exists
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file not found: $input_bam"
        return 1
    fi

    echo "Start processing sample: ${sample}.p${n}"

    # Step 1: smoove call
    echo "Running smoove call..."
    $software_path/smoove call \
        --outdir "$sample_dir" \
        --name "$sample" \
        --fasta "$reference_genome" \
        -p 1 \
        --genotype "$input_bam"

    # Step 2: smoove genotype
    echo "Running smoove genotype..."
    $software_path/smoove genotype -d -x \
        -p 1 \
        --name "${sample}.smoove" \
        --outdir "$smoove_dir/results-genotyped" \
        --fasta "$reference_genome" \
        --vcf "${sample_dir}/${sample}-smoove.genotyped.vcf.gz" \
        "$input_bam"

    # Step 3: Filter chromosomes
    echo "Filtering chromosomes..."
    input_vcf="${smoove_dir}/results-genotyped/${sample}.smoove-smoove.genotyped.vcf.gz"
    
    bcftools view \
        -r "$chromosomes" \
        "$input_vcf" \
        -Oz -o "${output_dir}/${sample}.smoove.vcf.gz"
    
    tabix "${output_dir}/${sample}.smoove.vcf.gz"

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