#!/bin/bash

#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="25"  # Different n values

manta_path="${BASE_DIR}/result/WGS-SVs/bamfile/manta"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
run_stats_dir="${BASE_DIR}/result/WGS-SVs/run"  # Runtime statistics directory
samtools_path="${BASE_DIR}/anaconda3/envs/manta/libexec/samtools"
convert_script="${BASE_DIR}/software/convertInversion.py"

# ========== Initialization ==========
mkdir -p "$run_stats_dir" "$manta_path" "$output_base"

# Check dependencies
check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed." >&2
        exit 1
    fi
}

# Check if required tools are available
check_dependency "samtools"
check_dependency "bcftools"
check_dependency "python"
check_dependency "/usr/bin/time"  # Ensure GNU time is available

# Activate Manta environment
source activate manta

# Check reference genome index
if [[ ! -f "${reference_genome}.fai" ]]; then
    samtools faidx "$reference_genome" || exit 1
fi

# ========== Sample processing function ==========
process_sample() {
    local sample="$1"
    local n="$2"
    
    local stats_file="${run_stats_dir}/manta/${sample}.p${n}.stats.txt"
    mkdir -p "$(dirname "$stats_file")"
    
    # Create output directories
    local input_bam="${input_base}/bam_p${n}/${sample}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    local sample_dir="${manta_path}/${sample}.p${n}"
    mkdir -p "$sample_dir" "$output_dir/${sample}"

    # Check if input file exists
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file not found: $input_bam" >> "$stats_file"
        return 1
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start processing sample: ${sample}.p${n}" >> "$stats_file"
    local start_time=$(date +%s)

    # Step 1: Configure Manta
    rm -rf "$sample_dir"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Configuring Manta..." >> "$stats_file"
    configManta.py \
        --bam "$input_bam" \
        --referenceFasta "$reference_genome" \
        --runDir "$sample_dir" >> "$stats_file" 2>&1

    # Step 2: Run Manta workflow (monitor memory using time)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Manta workflow..." >> "$stats_file"
    /usr/bin/time -v -o "${stats_file}.time" python "${sample_dir}/runWorkflow.py" -j 1 >> "$stats_file" 2>&1

    # Step 3: Process output
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing output..." >> "$stats_file"
    local manta_output_vcf="${sample_dir}/results/variants/diploidSV.vcf.gz"
    local manta_out_vcf="${sample_dir}/results/variants/manta.vcf"

    # Check if Manta output exists
    if [ ! -f "$manta_output_vcf" ]; then
        echo "Error: Manta output file not found: $manta_output_vcf" >> "$stats_file"
        return 1
    fi

    # BND2INV conversion
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running BND2INV conversion..." >> "$stats_file"
    python "$convert_script" \
        "$samtools_path" \
        "$reference_genome" \
        "$manta_output_vcf" > \
        "$manta_out_vcf" 2>> "$stats_file"

    # Create index
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating index..." >> "$stats_file"
    bgzip -f "$manta_out_vcf" >> "$stats_file" 2>&1
    tabix "$manta_out_vcf.gz" >> "$stats_file" 2>&1
   
    # Filter and compress
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filtering chromosomes..." >> "$stats_file"
    bcftools view -r "$chromosomes" \
        "$manta_out_vcf.gz" \
        -Oz -o "${output_dir}/${sample}/${sample}.manta.vcf.gz" >> "$stats_file" 2>&1

    # Create final index
    tabix -p vcf "${output_dir}/${sample}/${sample}.manta.vcf.gz" >> "$stats_file" 2>&1

    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    
    # Extract maximum memory usage from time output
    local max_mem=$(grep "Maximum resident set size" "${stats_file}.time" | awk '{print $6/1024}')
    
    # Append runtime and memory usage to stats file
    echo "===== Total runtime statistics =====" >> "$stats_file"
    echo "Total runtime: ${elapsed} seconds" >> "$stats_file"
    echo "Maximum memory usage: ${max_mem} MB" >> "$stats_file"
    
    # Append detailed time output
    echo "===== GNU time detailed statistics =====" >> "$stats_file"
    cat "${stats_file}.time" >> "$stats_file"
    rm "${stats_file}.time"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample ${sample}.p${n} completed, time: ${elapsed} seconds, max memory: ${max_mem} MB" >> "$stats_file"
    return 0
}

# ========== Main Process ==========

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list" >&2
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