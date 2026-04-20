#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_path="${BASE_DIR}/WGS/bam"
manta_path="${BASE_DIR}/result/WGS-SVs/manta"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari"
log_dir="${BASE_DIR}/error"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"

# ========== Initialization ==========
mkdir -p "$log_dir" "$manta_path" "$output_path"

# Check dependencies
check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed." >> "$log_dir/error.log"
        exit 1
    fi
}

# Check if required tools are available
check_dependency "samtools"
check_dependency "bcftools"
check_dependency "python"

# Activate conda environment
source activate manta

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list" >> "$log_dir/error.log"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Check reference genome index
if [[ ! -f "${reference_genome}.fai" ]]; then
    samtools faidx "$reference_genome" || exit 1
fi

# ========== Sample processing function ==========
process_sample() {
    local bamfile="$1"    
    
    # Check input file
    local bam_file="${input_path}/${bamfile}.WGS.bam"
    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found for $bamfile" >> "$log_dir/error.log"
        return 1
    fi

    # Create sample-specific directories
    local manta_sample_dir="${manta_path}/${bamfile}"
    rm -rf "${manta_sample_dir}"
    mkdir -p "$manta_sample_dir" "${output_path}/${bamfile}"

    # Step 1: Configure Manta
    configManta.py \
        --bam "$bam_file" \
        --referenceFasta "$reference_genome" \
        --runDir "$manta_sample_dir"

    # Step 2: Run Manta workflow
    python "${manta_sample_dir}/runWorkflow.py" -j 4

    # Step 3: Process output
    local manta_output_vcf="${manta_sample_dir}/results/variants/diploidSV.vcf.gz"
    local manta_out_vcf="${manta_sample_dir}/results/variants/manta.vcf"

    # BND2INV conversion
    ${BASE_DIR}/software/convertInversion.py \
        ${BASE_DIR}/anaconda3/envs/manta/libexec/samtools \
        "$reference_genome" \
        "$manta_output_vcf" > \
        "$manta_out_vcf"

    # Create index
    bgzip "$manta_out_vcf"
    tabix "$manta_out_vcf.gz"
   
    # Filter and compress
    bcftools view -r "$chromosomes" \
        "$manta_out_vcf.gz" \
        -Oz -o "${output_path}/${bamfile}/${bamfile}.manta.vcf.gz"

    # Create index
    tabix -p vcf "${output_path}/${bamfile}/${bamfile}.manta.vcf.gz"

    return 0
}

# ========== Main loop processing each sample ==========

# Read sample list and process in serial (no background &)
for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile"
done

# Wait for all remaining tasks to finish (though none in background)
wait