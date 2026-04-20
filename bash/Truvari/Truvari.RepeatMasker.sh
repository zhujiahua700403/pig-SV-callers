#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate truvari  # Activate conda environment

input_dir="${BASE_DIR}/result/WGS-SVs/region/vcf/RepeatMasker"
output_dir="${BASE_DIR}/result/WGS-SVs/Truvari/region/RepeatMasker/result"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
regions=("no")

# ========== Initialization ==========
mkdir -p "$output_dir"

# Check input files
validate_inputs() {
    [ ! -f "$bamfile_list" ] && { echo "ERROR: bamfiles.list not found"; exit 1; }
    [ ! -f "$reference_genome" ] && { echo "ERROR: Reference genome not found"; exit 1; }
    [ ! -d "$input_dir" ] && { echo "ERROR: Input directory not found"; exit 1; }
}

# ========== Processing functions ==========
process_sample() {
    local bamfile="$1"
    local software="$2"
    local region="$3"
    
    echo "Processing $bamfile - $software - $region"
    
    # Define file paths
    local base_vcf="$input_dir/$bamfile/base/$region/${bamfile}_base_${region}_overlap.vcf.gz"
    local comp_vcf="$input_dir/$bamfile/$software/$region/${bamfile}_${software}_${region}_overlap.vcf.gz"
    local output_subdir="$output_dir/$software/$bamfile"
    rm -rf "$output_subdir/$region"

    # Check input files
    if [ ! -f "$base_vcf" ]; then
        echo "ERROR: Base VCF not found: $base_vcf"
        return 1
    fi
    if [ ! -f "$comp_vcf" ]; then
        echo "ERROR: Comparison VCF not found: $comp_vcf"
        return 1
    fi
    
    # Index VCF files
    if [ ! -f "${base_vcf}.tbi" ]; then
        echo "Indexing base VCF..."
        tabix -p vcf "$base_vcf" || {
            echo "ERROR: Failed to index base VCF"
            return 1
        }
    fi
    
    if [ ! -f "${comp_vcf}.tbi" ]; then
        echo "Indexing comparison VCF..."
        tabix -p vcf "$comp_vcf" || {
            echo "ERROR: Failed to index comparison VCF"
            return 1
        }
    fi
    
    # Create output directory
    mkdir -p "$output_subdir" || {
        echo "ERROR: Failed to create output directory"
        return 1
    }
    
    # Run truvari
    echo "Running truvari comparison..."
    truvari bench \
        -b "$base_vcf" \
        -c "$comp_vcf" \
        -o "$output_subdir/$region" \
        -f "$reference_genome" \
        -r 1000 -p 0 -P 0.5 -O 0.5 \
        --passonly || {
        echo "ERROR: truvari failed"
        return 1
    }
    
    echo "Successfully processed $bamfile - $software - $region"
    return 0
}

# ========== Main Process ==========

echo "===== Starting Truvari region comparison ====="
echo "Start time: $(date)"
echo "Sample count: $(wc -l < "$bamfile_list")"
echo "Software count: ${#softwares[@]}"
echo "Region types: ${#regions[@]}"

# Read sample list
readarray -t bamfiles < "$bamfile_list"

# Process each sample and software, parallelizing across regions
for bamfile in "${bamfiles[@]}"; do
    for software in "${softwares[@]}"; do
        for region in "${regions[@]}"; do
            process_sample "$bamfile" "$software" "$region" &
        done
        wait  # Wait for all regions of current software to finish
    done
done

wait  # Wait for all software to finish

echo "===== Processing completed ====="
echo "End time: $(date)"