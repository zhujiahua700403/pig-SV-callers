#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate dysgu  # Activate conda environment

input_dir="${BASE_DIR}/WGS/bam"
output_dir="${BASE_DIR}/result/WGS-SVs/Truvari"
reference="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
dysgu_output="${BASE_DIR}/result/WGS-SVs/dysgu"
log_dir="${BASE_DIR}/error"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# ========== Initialization ==========
mkdir -p "$log_dir" "$dysgu_output" "$output_dir"

# ========== Sample processing function ==========
process_sample() {
    local sample="$1"
    bam_path="$input_dir/${sample}.WGS.bam"
    # Create directories
    mkdir -p "${output_dir}/${sample}" "${dysgu_output}/${sample}"

    # Step 1: Dysgu detection
    dysgu run \
        "$reference" \
        -x "${dysgu_output}/${sample}" \
        --min-size 50 \
        "$bam_path" \
        --overwrite > \
        "${dysgu_output}/${sample}/${sample}.dysgu.vcf"

    # Step 2: Convert TRA to BND
    python "${BASE_DIR}/script/WGS/tra2bnd.py" \
        -i "${dysgu_output}/${sample}/${sample}.dysgu.vcf" \
        -o "${dysgu_output}/${sample}/${sample}.dysgu2.vcf"

    # Step 3: Sort and compress
    bcftools sort \
        -o "${dysgu_output}/${sample}/${sample}.dysgu2.vcf.gz" \
        -O z \
        "${dysgu_output}/${sample}/${sample}.dysgu2.vcf"

    tabix -f "${dysgu_output}/${sample}/${sample}.dysgu2.vcf.gz"

    # Step 4: Filter by chromosomes
    bcftools view \
        -r "$chromosomes" \
        -f PASS \
        "${dysgu_output}/${sample}/${sample}.dysgu2.vcf.gz" \
        -Oz \
        -o "${output_dir}/${sample}/${sample}.dysgu.vcf.gz"

    tabix "${output_dir}/${sample}/${sample}.dysgu.vcf.gz"
    rm -f "${output_dir}/${sample}"  # Remove temporary directory
}

# ========== Main loop processing each sample ==========

# Read sample list into array
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Process each sample using a for loop
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile" &
done
wait