#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
ref_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="25"  # Different n values

insurveyor_dir="${BASE_DIR}/result/WGS-SVs/bamfile/insurveyor"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
run_stats_dir="${BASE_DIR}/result/WGS-SVs/run/insurveyor"  # Runtime statistics directory

# ========== Initialization ==========
mkdir -p "$run_stats_dir" "$insurveyor_dir" "$output_base"

# Activate InSurveyor environment
source activate insurveyor

# ========== Sample processing function ==========
process_sample() {
    local sample="$1"
    local n="$2"
    
    local start_time=$(date +%s)
    local stats_file="${run_stats_dir}/${sample}.p${n}.stats.txt"
    
    # Create output directories
    local input_bam="${input_base}/bam_p${n}/${sample}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    local sample_dir="${insurveyor_dir}/${sample}.p${n}"
    mkdir -p "$sample_dir" "$output_dir"

    # Check if input file exists
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file not found: $input_bam"
        return 1
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start processing sample: ${sample}.p${n}"

    # Monitor resource usage using GNU time (limit 1 thread)
    /usr/bin/time -v -o "$stats_file" \
    bash -c "
    # Step 1: InSurveyor detection
    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Running InSurveyor detection...\"
    insurveyor.py \
        \"$input_bam\" \
        \"$sample_dir\" \
        \"$ref_genome\" \
        --min-insertion-size 50 \
        --threads 1

    # Step 2: Filter by chromosomes
    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Filtering chromosomes...\"
    bcftools view -r \"$chromosomes\" \
        \"$sample_dir/out.pass.vcf.gz\" \
        -Oz -o \"${output_dir}/${sample}.insurveyor.vcf.gz\"
    tabix \"${output_dir}/${sample}.insurveyor.vcf.gz\"
    "

    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    
    # Extract maximum memory usage from stats file
    local max_mem=$(grep "Maximum resident set size" "$stats_file" | awk '{print $6/1024/1024}')
    
    # Append runtime and memory usage to stats file
    echo "===== Total runtime statistics =====" >> "$stats_file"
    echo "Total runtime: ${elapsed} seconds" >> "$stats_file"
    echo "Maximum memory usage: ${max_mem} GB" >> "$stats_file"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample ${sample}.p${n} completed, time: ${elapsed} seconds, max memory: ${max_mem} GB"
    
    # Clean up temporary files
    rm -rf "$sample_dir"
    return 0
}

# ========== Main Process ==========

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi

# Read sample list
readarray -t samples < "$bamfile_list"

# Process each sample for each n value
for sample in "${samples[@]}"; do
    # Remove possible whitespace
    sample=$(echo "$sample" | xargs)
    
    # Skip empty lines
    if [[ -z "$sample" ]]; then
        continue
    fi
    
    # Process sample for each n value
    for n in $n_values; do
        echo "Processing sample: $sample using n=$n"
        process_sample "$sample" "$n"
    done
done

echo "All samples processed"