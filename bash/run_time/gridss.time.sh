#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="25"  # Different n values

gridss_path="${BASE_DIR}/result/WGS-SVs/bamfile/gridss"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
run_stats_dir="${BASE_DIR}/result/WGS-SVs/run/gridss"  # Runtime statistics directory
software="${BASE_DIR}/software"
gridss_jar="$software/gridss-2.13.2-gridss-jar-with-dependencies.jar"

# ========== Initialization ==========
mkdir -p "$run_stats_dir" "$gridss_path" "$output_base"

# Activate Java environment
source activate java

# ========== Sample processing function ==========
process_sample() {
    local sample="$1"
    local n="$2"
    
    local start_time=$(date +%s)
    local stats_file="${run_stats_dir}/${sample}.p${n}.stats.txt"
    
    # Create output directories
    local input_bam="${input_base}/bam_p${n}/${sample}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    local sample_dir="${gridss_path}/${sample}.p${n}"
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
    # 1. Preprocessing
    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Preprocessing...\"
    cd \"$sample_dir\"
    $software/gridss \
        -r \"$reference_genome\" \
        -j \"$gridss_jar\" \
        -s preprocess \
        -t 1 \
        \"$input_bam\"

    # 2. Distributed assembly
    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Distributed assembly...\"
    for jobindex in {0..2}; do
        $software/gridss \
            -r \"$reference_genome\" \
            -j \"$gridss_jar\" \
            -t 1 \
            -s assemble \
            -a assembly.bam \
            --jobnodes 3 \
            --jobindex \"\$jobindex\" \
            \"$input_bam\"
    done

    # 3. Final assembly and calling
    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Final assembly and calling...\"
    $software/gridss \
        -r \"$reference_genome\" \
        -j \"$gridss_jar\" \
        -t 1 \
        -a assembly.bam \
        -s assemble,call \
        -o \"${sample_dir}/${sample}.gridss.vcf\" \
        \"$input_bam\"

    # 4. Filter PASS variants
    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Filtering PASS variants...\"
    bcftools view \
        -f PASS \
        -o \"${sample_dir}/${sample}.filtered.gridss.vcf\" \
        \"${sample_dir}/${sample}.gridss.vcf\"

    # 5. Compress and filter by chromosomes
    echo \"[$(date '+%Y-%m-%d %H:%M:%S')] Compressing and filtering chromosomes...\"
    bgzip -f \"${sample_dir}/${sample}.filtered.gridss.vcf\"
    tabix \"${sample_dir}/${sample}.filtered.gridss.vcf.gz\"
    
    # Filter by specified chromosomes
    bcftools view \
        -r \"$chromosomes\" \
        \"${sample_dir}/${sample}.filtered.gridss.vcf.gz\" \
        -Oz \
        -o \"${output_dir}/${sample}.gridss.vcf.gz\"
    tabix \"${output_dir}/${sample}.gridss.vcf.gz\"
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