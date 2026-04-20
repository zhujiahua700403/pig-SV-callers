#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate delly  # Activate conda environment

input_path="${BASE_DIR}/WGS/bam"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari"
delly_path="${BASE_DIR}/result/WGS-SVs/delly"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
log_dir="${BASE_DIR}/error"  # Log directory
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# ========== Initialization ==========
mkdir -p "$log_dir" "$delly_path" "$output_path"

# ========== Main processing function ==========
process_sample() {
    local bamfile="$1"
    
    # Create sample-specific directories
    mkdir -p "${delly_path}/${bamfile}" "${output_path}/${bamfile}"

    # Step 1: Delly SV detection
    echo "[$(date '+%T')] Running Delly detection..."
    delly call \
        -g "$reference_genome" \
        -o "${delly_path}/${bamfile}/${bamfile}.bcf" \
        "${input_path}/${bamfile}.WGS.bam"

    # Step 2: Filter PASS variants
    bcftools view \
        -f PASS \
        -r "$chromosomes" \
        "${delly_path}/${bamfile}/${bamfile}.bcf" \
        > "${output_path}/${bamfile}/${bamfile}.delly.vcf"

    # Step 3: Compress and index
    bgzip -f "${output_path}/${bamfile}/${bamfile}.delly.vcf"
    tabix "${output_path}/${bamfile}/${bamfile}.delly.vcf.gz"
}

# Read sample list into array
readarray -t bamfiles < "$bamfile_list"

# Process each sample
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile" &
done
wait

exit 0