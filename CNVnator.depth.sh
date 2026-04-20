#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/WGS/bamfile"
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
ref_fa="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
ref_dir="${BASE_DIR}/genome"
n_values="1 5 10 15 20 25 30 45 60"

cnvnator_dir="${BASE_DIR}/result/WGS-SVs/bamfile/CNVnator"
output_base="${BASE_DIR}/result/WGS-SVs/bamfile"
run_stats_dir="${BASE_DIR}/result/WGS-SVs/run/CNVnator"

# Create necessary directories
mkdir -p "$run_stats_dir" "$cnvnator_dir" "$output_base"

# Activate CNVnator environment
source activate cnvnator

# Check reference genome index
if [[ ! -f "${ref_fa}.fai" ]]; then
    samtools faidx "$ref_fa" || { echo "Error: Failed to index reference genome"; exit 1; }
fi

# Function to determine bin_size based on n value
get_bin_size() {
    local n=$1
    case $n in
        1) echo 1000 ;;
        5|10|15) echo 500 ;;
        20|25|30|45|60) echo 100 ;;
        *) echo 100 ;; # Default value
    esac
}

# Define function to process a single sample
process_sample() {
    local bamfile="$1"
    local n="$2"
    local bin_size=$(get_bin_size $n)
    
    # Input/output paths
    local input_bam="${input_base}/bam_p${n}/${bamfile}.p${n}.bam"
    local output_dir="${output_base}/p${n}"
    local sample_dir="${cnvnator_dir}/${bamfile}.p${n}"
    
    mkdir -p "$sample_dir" "$output_dir"

    # Check input file
    if [ ! -f "$input_bam" ]; then
        echo "Error: Input BAM file not found: $input_bam"
        return 1
    fi

    echo "Start processing sample: ${bamfile}.p${n} (bin_size=${bin_size})"
    
    # Run commands
    set -e
    
    root_file="${sample_dir}/${bamfile}.root"
    call_file="${sample_dir}/${bamfile}.cnv.call.txt"

    # 1. Generate analysis tree
    echo "Generating analysis tree..."
    cnvnator -root "$root_file" -tree "$input_bam" -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

    # 2. Generate histogram
    echo "Generating histogram..."
    cnvnator -root "$root_file" -his $bin_size -fasta "$ref_fa"

    # 3. Statistics calculation
    echo "Calculating statistics..."
    cnvnator -root "$root_file" -stat $bin_size

    # 4. Data partitioning
    echo "Partitioning data..."
    cnvnator -root "$root_file" -partition $bin_size

    # 5. CNV calling
    echo "CNV calling..."
    cnvnator -root "$root_file" -call $bin_size > "$call_file"

    # Check if call file was generated
    if [ ! -s "$call_file" ]; then
        echo "Error: CNV call file is empty or not generated: $call_file"
        return 1
    fi

    # 6. Format conversion
    echo "Converting format..."
    cnvnator2VCF.pl -prefix "${bamfile}.p${n}" -reference "$ref_fa" \
        "$call_file" "$ref_dir" > "${output_dir}/${bamfile}.CNVnator.vcf"

    # Check if VCF file was generated
    if [ ! -s "${output_dir}/${bamfile}.CNVnator.vcf" ]; then
        echo "Error: VCF file is empty or not generated"
        return 1
    fi

    # 7. Compress and index
    echo "Compressing and indexing VCF file..."
    bgzip -f "${output_dir}/${bamfile}.CNVnator.vcf"
    tabix -p vcf "${output_dir}/${bamfile}.CNVnator.vcf.gz"

    echo "Sample ${bamfile}.p${n} completed, bin_size=${bin_size}"
    return 0
}

# Main program
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list" >&2
    exit 1
fi

readarray -t samples < "$bamfile_list"

for sample in "${samples[@]}"; do
    sample=$(echo "$sample" | xargs)
    [[ -z "$sample" ]] && continue
    
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