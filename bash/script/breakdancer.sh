#!/bin/bash
set -e  # Enable error exit mechanism, script exits if any command fails

# ========== Configuration Parameters ==========
input_paths="/path/to/your/folder/bam"
bamfile_list="/path/to/your/folder/bamfiles.list"
breakdancer_path="/path/to/your/folder/breakdancer"
reference_genome="/path/to/your/folder/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_path="/path/to/your/folder/Truvari"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"  # Chromosome list

log_dir="/home/zhujiahua/error"  # Log directory
mkdir -p "$log_dir"

# Activate BreakDancer environment
source activate breakdancer

# Define function to process a single sample
process_bamfile() {
    local bamfile="$1"

    # Create output directories
    mkdir -p "${breakdancer_path}/${bamfile}" "${output_path}/${bamfile}"

    # Generate config file
    bam2cfg.pl "${input_paths}/${bamfile}.WGS.bam" > "${breakdancer_path}/${bamfile}/${bamfile}.config.cfg"

    # Run BreakDancer
    breakdancer-max \
        -d sv.reads \
        "${breakdancer_path}/${bamfile}/${bamfile}.config.cfg" \
        > "${breakdancer_path}/${bamfile}/${bamfile}.out"

    # Convert to VCF format
    python /home/zhujiahua/script/WGS/breakdancer2vcf.py \
        -i "${breakdancer_path}/${bamfile}/${bamfile}.out" \
        -o "${breakdancer_path}/${bamfile}/${bamfile}.breakdancer.vcf" \
        -s "$bamfile"

    echo "Processing completed, output file: ${breakdancer_path}/${bamfile}/${bamfile}.breakdancer.vcf"

    # Sort, compress, and index VCF file
    bcftools sort \
        -o "${breakdancer_path}/${bamfile}/${bamfile}.breakdancer.vcf.gz" \
        -O z \
        "${breakdancer_path}/${bamfile}/${bamfile}.breakdancer.vcf"
    tabix "${breakdancer_path}/${bamfile}/${bamfile}.breakdancer.vcf.gz"

    # Filter by specified chromosomes
    bcftools view \
        -r "${chromosomes}" \
        "${breakdancer_path}/${bamfile}/${bamfile}.breakdancer.vcf.gz" \
        -Oz \
        -o "${output_path}/${bamfile}/${bamfile}.breakdancer.vcf.gz"
    tabix "${output_path}/${bamfile}/${bamfile}.breakdancer.vcf.gz"
}

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Iterate over all samples and process in parallel
for bamfile in "${bamfiles[@]}"; do
    process_bamfile "$bamfile"
done

# Wait for all background processes to finish
wait