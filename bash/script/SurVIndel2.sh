#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate SurVIndel2  # Activate conda environment

input_dir="${BASE_DIR}/WGS/bam"
survindel_dir="${BASE_DIR}/result/WGS-SVs/SurVIndel2"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_dir="${BASE_DIR}/result/WGS-SVs/Truvari"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
log_dir="${BASE_DIR}/error"

# ========== Initialization ==========
mkdir -p "$log_dir" "$survindel_dir" "$output_dir"

# Read sample list into array
readarray -t bamfiles < "$bamfile_list"

# ========== Sample processing function ==========
process_sample() {
    local bamfile="$1"

    # Create sample-specific directories
    mkdir -p "${survindel_dir}/${bamfile}" "${output_dir}/${bamfile}"
    bam_path="${input_dir}/${bamfile}.WGS.bam"

    # Step 1: Run SurVIndel2
    survindel2.py \
        "$bam_path" \
        "${survindel_dir}/${bamfile}" \
        "$reference_genome" \
        --threads 4 \
        --min_sv_size 50 \
        --samplename "${bamfile}.SurVIndel2"

    # Step 2: Filter by chromosomes
    local input_vcf="${survindel_dir}/${bamfile}/out.pass.vcf.gz"
    
    bcftools view \
        -r "$chromosomes" \
        "$input_vcf" \
        -Oz -o "${output_dir}/${bamfile}/${bamfile}.SurVIndel2.vcf.gz"

    tabix "${output_dir}/${bamfile}/${bamfile}.SurVIndel2.vcf.gz"
}

# Process each sample serially
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile"
done

wait