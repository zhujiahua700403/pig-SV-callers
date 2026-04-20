#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate wham  # Activate conda environment

input_dir="${BASE_DIR}/WGS/bam"
output_dir="${BASE_DIR}/result/WGS-SVs/Truvari"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
log_dir="${BASE_DIR}/error"

# ========== Initialization ==========
mkdir -p "$log_dir" "$output_dir"

# ========== Sample processing function ==========
process_sample() {
    local bamfile="$1"
        
    # Check input file
    local bam_path="$input_dir/${bamfile}.WGS.bam"

    if [[ ! -f "$bam_path" ]]; then
        echo "[ERROR] BAM file does not exist: $bam_path"
        return 1
    fi

    # Create output directory
    mkdir -p "${output_dir}/${bamfile}"

    # Step 1: Run WHAM
    whamg -x 2 \
        -f "$bam_path" \
        -a "$reference_genome" \
        -c "$chromosomes" \
        > "${output_dir}/${bamfile}/${bamfile}.wham.vcf"

    # Step 2: Compress and index VCF
    bgzip -f "${output_dir}/${bamfile}/${bamfile}.wham.vcf"
    tabix "${output_dir}/${bamfile}/${bamfile}.wham.vcf.gz"
}

# ========== Main loop processing each sample ==========

# Read sample list into array
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Process each sample in parallel
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile" &
done

wait