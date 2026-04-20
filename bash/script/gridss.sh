#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate java  # Activate conda environment

input_path="${BASE_DIR}/WGS/bam"
gridss_path="${BASE_DIR}/result/WGS-SVs/gridss"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari"
software="${BASE_DIR}/software"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
log_dir="${BASE_DIR}/error"
gridss_jar="$software/gridss-2.13.2-gridss-jar-with-dependencies.jar"

# ========== Initialization ==========
mkdir -p "$log_dir" "$gridss_path" "$output_path"

# ========== Sample processing function ==========
process_sample() {
    local bamfile="$1"

    # Create output directories
    mkdir -p "${gridss_path}/${bamfile}" "${output_path}/${bamfile}"
    
    # 1. Preprocessing
    cd "${gridss_path}/${bamfile}"

    $software/gridss \
        -r "$reference_genome" \
        -j "$gridss_jar" \
        -s preprocess \
        -t 8 \
        "${input_path}/${bamfile}.WGS.bam"

    # 2. Distributed assembly
    for jobindex in {0..2}; do
        $software/gridss \
            -r "$reference_genome" \
            -j "$gridss_jar" \
            -t 8 \
            -s assemble \
            -a assembly.bam \
            --jobnodes 3 \
            --jobindex "$jobindex" \
            "${input_path}/${bamfile}.WGS.bam"
    done

    # 3. Final assembly and calling
    $software/gridss \
        -r "$reference_genome" \
        -j "$gridss_jar" \
        -t 8 \
        -a assembly.bam \
        -s assemble,call \
        -o "${gridss_path}/${bamfile}/${bamfile}.gridss.vcf" \
        "${input_path}/${bamfile}.WGS.bam"

    # 4. Filter PASS variants
    bcftools view \
        -f PASS \
        -o "${gridss_path}/${bamfile}/${bamfile}.filtered.gridss.vcf" \
        "${gridss_path}/${bamfile}/${bamfile}.gridss.vcf"

    # 5. Compress and filter by chromosomes
    bgzip -f "${gridss_path}/${bamfile}/${bamfile}.filtered.gridss.vcf"
    tabix "${gridss_path}/${bamfile}/${bamfile}.filtered.gridss.vcf.gz"
    
    # Filter by specified chromosomes
    bcftools view \
        -r "$chromosomes" \
        "${gridss_path}/${bamfile}/${bamfile}.filtered.gridss.vcf.gz" \
        -Oz \
        -o "${output_path}/${bamfile}/${bamfile}.gridss.vcf.gz"

    tabix "${output_path}/${bamfile}/${bamfile}.gridss.vcf.gz"

    echo "[$(date '+%T')] Successfully processed: $bamfile"
    return 0
}

# ========== Main loop processing each sample ==========

# Read sample list into array
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Process each sample using a for loop (serial processing)
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile"
done