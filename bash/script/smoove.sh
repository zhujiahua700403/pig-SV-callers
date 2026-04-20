#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate smoove  # Activate conda environment

input_dir="${BASE_DIR}/WGS/bam"
output_dir="${BASE_DIR}/result/WGS-SVs/Truvari"
smoove_dir="${BASE_DIR}/result/WGS-SVs/smoove"
software_path="${BASE_DIR}/software"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
bamfile_list="${BASE_DIR}/result/WGS-SVs/1.list"
log_dir="${BASE_DIR}/error"

# ========== Initialization ==========
mkdir -p "$log_dir" "$smoove_dir" "$output_dir" "$smoove_dir/results-genotyped"

# ========== Sample processing function ==========
process_sample() {
    local bamfile="$1"
    
    # Create sample-specific directories
    mkdir -p "${smoove_dir}/${bamfile}" "${output_dir}/${bamfile}"

    # Check input file
    local bam_path="${input_dir}/${bamfile}.WGS.bam"

    if [[ ! -f "$bam_path" ]]; then
        echo "Error: BAM file not found: $bam_path" >> "$log_dir/error.log"
        return 1
    fi

    # Step 1: smoove call
    echo "Calling variants for $bamfile..." >> "$log_dir/error.log"
    if ! $software_path/smoove call \
        --outdir "${smoove_dir}/${bamfile}" \
        --name "$bamfile" \
        --fasta "$reference_genome" \
        -p 3 \
        --genotype "$bam_path" >> "$log_dir/process.log" 2>&1; then
        echo "Error: smoove call failed for $bamfile" >> "$log_dir/smoove.log"
        return 1
    fi

    # Step 2: smoove genotype
    echo "Genotyping variants for $bamfile..." >> "$log_dir/error.log"
    $software_path/smoove genotype -d -x \
        -p 3 \
        --name "${bamfile}.smoove" \
        --outdir "$smoove_dir/results-genotyped" \
        --fasta "$reference_genome" \
        --vcf "${smoove_dir}/${bamfile}/${bamfile}-smoove.genotyped.vcf.gz" \
        "$bam_path"

    # Step 3: Filter chromosomes
    input_vcf="${smoove_dir}/results-genotyped/${bamfile}.smoove-smoove.genotyped.vcf.gz"
    
    echo "Filtering chromosomes for $bamfile..." >> "$log_dir/error.log"
    if ! bcftools view \
        -r "$chromosomes" \
        "$input_vcf" \
        -Oz -o "${output_dir}/${bamfile}/${bamfile}.smoove.vcf.gz" >> "$log_dir/error.log" 2>&1; then
        echo "Error: bcftools view failed for $bamfile" >> "$log_dir/smoove.log"
        return 1
    fi

    tabix "${output_dir}/${bamfile}/${bamfile}.smoove.vcf.gz"
}

# Read sample list into array
if [ ! -f "$bamfile_list" ]; then
    echo "smoove: bamfiles.list not found at $bamfile_list" >> "$log_dir/smoove.log"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Process each sample in parallel
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile" &
done
wait