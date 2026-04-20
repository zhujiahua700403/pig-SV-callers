#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="1 5 10 15 20 25 30 45 60"  # Different n values

gridss_path="${BASE_DIR}/result/WGS-SVs/bamfile/gridss"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
software="${BASE_DIR}/software"
gridss_jar="$software/gridss-2.13.2-gridss-jar-with-dependencies.jar"

# ========== Initialization ==========
mkdir -p "$gridss_path" "$output_base"

# Activate Java environment
source activate java

# ========== Sample processing function ==========
process_sample() {
    local sample="$1"
    local n="$2"
    
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

    echo "Start processing sample: ${sample}.p${n}"
    
    # 1. Preprocessing
    echo "Preprocessing..."
    cd "$sample_dir"
    $software/gridss \
        -r "$reference_genome" \
        -j "$gridss_jar" \
        -s preprocess \
        -t 8 \
        "$input_bam"

    # 2. Distributed assembly
    echo "Distributed assembly..."
    for jobindex in {0..2}; do
        $software/gridss \
            -r "$reference_genome" \
            -j "$gridss_jar" \
            -t 8 \
            -s assemble \
            -a assembly.bam \
            --jobnodes 3 \
            --jobindex "$jobindex" \
            "$input_bam"
    done

    # 3. Final assembly and calling
    echo "Final assembly and calling..."
    $software/gridss \
        -r "$reference_genome" \
        -j "$gridss_jar" \
        -t 8 \
        -a assembly.bam \
        -s assemble,call \
        -o "${sample_dir}/${sample}.gridss.vcf" \
        "$input_bam"

    # 4. Filter PASS variants
    echo "Filtering PASS variants..."
    bcftools view \
        -f PASS \
        -o "${sample_dir}/${sample}.filtered.gridss.vcf" \
        "${sample_dir}/${sample}.gridss.vcf"

    # 5. Compress and filter by chromosomes
    echo "Compressing and filtering chromosomes..."
    bgzip -f "${sample_dir}/${sample}.filtered.gridss.vcf"
    tabix "${sample_dir}/${sample}.filtered.gridss.vcf.gz"
    
    # Filter by specified chromosomes
    bcftools view \
        -r "$chromosomes" \
        "${sample_dir}/${sample}.filtered.gridss.vcf.gz" \
        -Oz \
        -o "${output_dir}/${sample}.gridss.vcf.gz"
    tabix "${output_dir}/${sample}.gridss.vcf.gz"

    echo "Sample ${sample}.p${n} completed"
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
        if process_sample "$sample" "$n"; then
            echo "Successfully processed: $sample.p${n}"
        else
            echo "Processing failed: $sample.p${n}" >&2
        fi
    done
done

echo "All samples processed"