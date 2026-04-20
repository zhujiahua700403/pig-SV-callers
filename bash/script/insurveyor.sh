#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate insurveyor  # Activate conda environment

input_dir="${BASE_DIR}/WGS/bam"
ref_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_dir="${BASE_DIR}/result/WGS-SVs/Truvari"
insurveyor_dir="${BASE_DIR}/result/WGS-SVs/insurveyor"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
log_dir="${BASE_DIR}/error"

# ========== Initialization ==========
mkdir -p "$log_dir" "$insurveyor_dir" "$output_dir"

# ========== Sample processing function ==========
process_sample() {
    local bamfile="$1"

    # Create directory structure
    mkdir -p "$insurveyor_dir/$bamfile" "$output_dir/$bamfile"
    
    bam_path="$input_dir/$bamfile.WGS.bam"
    
    # Step 1: Insurveyor detection
    insurveyor.py \
        "$bam_path" \
        "$insurveyor_dir/$bamfile" \
        "$ref_genome" \
        --min-insertion-size 50 \
        --threads 4

    # Step 2: Filter by chromosomes
    bcftools view -r "$chromosomes" \
        "$insurveyor_dir/$bamfile/out.pass.vcf.gz" \
        -Oz -o "$output_dir/$bamfile/$bamfile.insurveyor.vcf.gz"
    tabix "$output_dir/$bamfile/$bamfile.insurveyor.vcf.gz"

    # Remove intermediate directory
    rm -f "$insurveyor_dir/$bamfile"
}

# ========== Main loop processing each sample ==========

# Read sample list into array
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Process each sample serially (no background)
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile"
done

wait