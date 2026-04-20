#!/bin/bash
set -euo pipefail  # Enable strict error handling

# Activate conda environment
source activate truvari

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/result/WGS-SVs/bamfile"
output_base="${BASE_DIR}/result/WGS-SVs/Truvari/depth"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "insurveyor" "manta" "smoove" "SurVIndel2" "wham")
bamfile_list="${BASE_DIR}/result/WGS-SVs/depth.list"
n_values=(1 5 10 15 20 25 30 45 60)  # Defined n values array

# ========== Initialization Checks ==========
[ ! -f "$bamfile_list" ] && { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: bamfiles.list not found" >&2; exit 1; }
[ ! -f "$reference_genome" ] && { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Reference genome not found" >&2; exit 1; }
[ ! -d "$input_base" ] && { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Input directory not found" >&2; exit 1; }

# ========== Function Definitions ==========
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

check_vcf() {
    local vcf_file="$1"
    local vcf_type="$2"
    
    [ ! -f "$vcf_file" ] && { log "ERROR: $vcf_type VCF not found: $vcf_file"; return 1; }
    
    if [ ! -f "${vcf_file}.tbi" ]; then
        log "Indexing $vcf_type VCF: $vcf_file"
        tabix -p vcf "$vcf_file" || { log "ERROR: Failed to index $vcf_type VCF"; return 1; }
    fi
    
    return 0
}

run_truvari_for_n() {
    local base_vcf="$1"
    local comp_vcf="$2"
    local output_dir="$3"
    local bamfile="$4"
    local software="$5"
    local n="$6"

    rm -rf "$output_dir/${n}"

    log "Running truvari for $bamfile - $software (p$n)"
    if ! truvari bench \
        -b "$base_vcf" \
        -c "$comp_vcf" \
        -o "$output_dir/${n}" \
        -f "$reference_genome" \
        --sizemax 100000 \
        -r 1000 -p 0 -P 0.5 -O 0.5 \
        --passonly > /dev/null 2>&1; then
        log "ERROR: truvari failed for $bamfile - $software (p$n)"
        return 1
    fi
    
    return 0
}

# ========== Main Process ==========
log "Starting processing pipeline"

# Create main output directory
mkdir -p "$output_base" || { log "ERROR: Failed to create main output directory"; exit 1; }

# Read sample list
mapfile -t bamfiles < "$bamfile_list"

# Process each sample and software sequentially
for bamfile in "${bamfiles[@]}"; do
    for software in "${softwares[@]}"; do
        log "Processing $bamfile with $software"
        
        # Prepare base VCF path
        base_dir="${BASE_DIR}/result/WGS-SVs/Truvari/${bamfile}"
        base_vcf="${base_dir}/${bamfile}.base.vcf.gz"
        
        # Check base VCF once per software
        check_vcf "$base_vcf" "base" || continue
        
        # Create output directory for this software
        output_dir="${output_base}/${bamfile}/${software}"
        mkdir -p "$output_dir" || { log "ERROR: Failed to create output directory for $software"; continue; }
        
        # Process all n_values in parallel
        for n in "${n_values[@]}"; do
            input_dir="${input_base}/p${n}"
            [ ! -d "$input_dir" ] && { log "WARNING: Input directory not found: $input_dir, skipping"; continue; }
            
            comp_vcf="${input_dir}/${bamfile}.${software}.vcf.gz"
            check_vcf "$comp_vcf" "comparison" || continue
            
            # Run truvari in background for parallel processing
            run_truvari_for_n "$base_vcf" "$comp_vcf" "$output_dir" "$bamfile" "$software" "$n" &
        done
        
        # Wait for all parallel jobs to finish before moving to next software
        wait
    done
done

log "All processing completed successfully!"