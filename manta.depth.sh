#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"
n_values="5 10 15 20 25 30 45 60"  # Different n values

manta_path="${BASE_DIR}/result/WGS-SVs/bamfile/manta"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
samtools_path="${BASE_DIR}/anaconda3/envs/manta/libexec/samtools"
convert_script="${BASE_DIR}/software/convertInversion.py"

# ========== Initialization ==========
mkdir -p "$manta_path" "$output_base"

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
    
    # Create output directories
    local input_bam="${input_base}/bam_p${n}/${sample}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    local sample_dir="${manta_path}/${sample}.p${n}"
    mkdir -p "$sample_dir" "$output_dir/${sample}"

    # Check if input file exists
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file not found: $input_bam"
        return 1
    fi

    echo "Start processing sample: ${sample}.p${n}"

    # Step 1: Configure Manta
    rm -rf "$sample_dir"
    echo "Configuring Manta..."
    configManta.py \
        --bam "$input_bam" \
        --referenceFasta "$reference_genome" \
        --runDir "$sample_dir"

    # Step 2: Run Manta workflow
    echo "Running Manta workflow..."
    python "${sample_dir}/runWorkflow.py" -j 8

    # Step 3: Process output
    echo "Processing output..."
    local manta_output_vcf="${sample_dir}/results/variants/diploidSV.vcf.gz"
    local manta_out_vcf="${sample_dir}/results/variants/manta.vcf"

    # Check if Manta output exists
    if [ ! -f "$manta_output_vcf" ]; then
        echo "Error: Manta output file not found: $manta_output_vcf"
        return 1
    fi

    # BND2INV conversion
    echo "Running BND2INV conversion..."
    python "$convert_script" \
        "$samtools_path" \
        "$reference_genome" \
        "$manta_output_vcf" > \
        "$manta_out_vcf"

    # Create index
    echo "Creating index..."
    bgzip -f "$manta_out_vcf"
    tabix "$manta_out_vcf.gz"
   
    # Filter and compress
    echo "Filtering chromosomes..."
    bcftools view -r "$chromosomes" \
        "$manta_out_vcf.gz" \
        -Oz -o "${output_dir}/${sample}.manta.vcf.gz"

    # Create final index
    tabix -p vcf "${output_dir}/${sample}.manta.vcf.gz"

    echo "Sample ${sample}.p${n} completed"
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
        if process_sample "$sample" "$n"; then
            echo "Successfully processed: $sample.p${n}"
        else
            echo "Processing failed: $sample.p${n}" >&2
        fi
    done
done

echo "All samples processed"