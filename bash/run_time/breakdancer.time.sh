#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
breakdancer_path="${BASE_DIR}/result/WGS-SVs/bamfile/breakdancer"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"  # Chromosome list
n_values="25"  # Different n values
run_stats_dir="${BASE_DIR}/result/WGS-SVs/run/breakdancer"  # Runtime statistics directory

# Create necessary directories
mkdir -p "$breakdancer_path" "$run_stats_dir" "$output_base"

# Activate BreakDancer environment
source activate breakdancer || {
    echo "Error: Failed to activate breakdancer environment"
    exit 1
}

# Define function to process a single sample
process_bamfile() {
    local bamfile="$1"
    local n="$2"
    
    local start_time=$(date +%s)
    local sample_id="${bamfile}.p${n}"
    local stats_file="${run_stats_dir}/${sample_id}.stats.txt"
    
    # Input/output paths
    local input_bam="${input_base}/bam_p${n}/${sample_id}.bam"
    local output_dir="${output_base}/p${n}"
    local work_dir="${breakdancer_path}/${sample_id}"
    
    mkdir -p "$work_dir" "$output_dir"

    # Check input file
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file does not exist: $input_bam" >&2
        return 1
    fi

    # Resource monitoring and command execution
    /usr/bin/time -v -o "$stats_file" \
    bash -c "
    set -euo pipefail
        
    # 1. Generate config file
    bam2cfg.pl \"$input_bam\" > \"${work_dir}/${bamfile}.config.cfg\" || exit 1
        
    # 2. Run BreakDancer (unlimited threads)
    breakdancer-max \
        -d \"${work_dir}/sv.reads\" \
        \"${work_dir}/${bamfile}.config.cfg\" \
        > \"${work_dir}/${bamfile}.out\" || exit 1
        
    # 3. Convert to VCF format
    python ${BASE_DIR}/script/WGS/breakdancer2vcf.py \
        -i \"${work_dir}/${bamfile}.out\" \
        -o \"${work_dir}/${bamfile}.breakdancer.vcf\" \
        -s \"$bamfile\" || exit 1
        
    # 4. Process VCF file
    bcftools sort \"${work_dir}/${bamfile}.breakdancer.vcf\" \
        -Oz -o \"${work_dir}/${bamfile}.breakdancer.vcf.gz\" || exit 1
    tabix -p vcf \"${work_dir}/${bamfile}.breakdancer.vcf.gz\" || exit 1
        
    # 5. Filter by specified chromosomes
    bcftools view -r \"$chromosomes\" \
        \"${work_dir}/${bamfile}.breakdancer.vcf.gz\" \
        -Oz -o \"${output_dir}/${bamfile}.breakdancer.vcf.gz\" || exit 1
    tabix -p vcf \"${output_dir}/${bamfile}.breakdancer.vcf.gz\" || exit 1
    "

    # Calculate runtime and memory usage
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    local max_mem=$(grep "Maximum resident set size" "$stats_file" | awk '{printf "%.2f", $6/1024/1024}')
    
    # Record statistics
    echo -e "\n===== Total runtime statistics =====" >> "$stats_file"
    echo "Total runtime: ${elapsed} seconds" >> "$stats_file"
    echo "Maximum memory usage: ${max_mem} GB" >> "$stats_file"
    
    echo "Sample ${sample_id} completed, time: ${elapsed} seconds, max memory: ${max_mem} GB"
}

# Main processing flow
echo "===== Starting sample processing ====="

# Read and process samples
while IFS= read -r bamfile || [[ -n "$bamfile" ]]; do
    bamfile=$(echo "$bamfile" | xargs)  # Remove whitespace
    [[ -z "$bamfile" ]] && continue  # Skip empty lines
    
    for n in $n_values; do
        echo "Processing: $bamfile.p${n}"
        if process_bamfile "$bamfile" "$n"; then
            echo "Successfully processed: $bamfile.p${n}"
        else
            echo "Processing failed: $bamfile.p${n}" >&2
        fi
    done
done < "$bamfile_list"

echo "===== All samples processed ====="