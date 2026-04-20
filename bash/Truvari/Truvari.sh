#!/bin/bash
set -euo pipefail  # Enable strict error handling

# Activate conda environment
source activate truvari

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_paths="${BASE_DIR}/result/WGS-SVs/Truvari"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari/result"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# ========== Initialization Checks ==========
[ ! -f "$bamfile_list" ] && { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: bamfiles.list not found" >&2; exit 1; }
[ ! -f "$reference_genome" ] && { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Reference genome not found" >&2; exit 1; }
[ ! -d "$input_paths" ] && { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Input directory not found" >&2; exit 1; }

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

prepare_output_dir() {
    local dir="$1"
    
    if [ -d "$dir" ]; then
        log "WARNING: Output directory exists, removing: $dir"
        rm -rf "$dir" || { log "ERROR: Failed to remove existing directory"; return 1; }
    fi
    
    mkdir -p "$dir" || { log "ERROR: Failed to create output directory"; return 1; }
    return 0
}

run_truvari() {
    local base_vcf="$1"
    local comp_vcf="$2"
    local output_dir="$3"
    local bamfile="$4"
    local software="$5"
    
    log "Running truvari for $bamfile - $software"
    local log_file="${output_dir}/truvari.log"
    
    if ! truvari bench \
        -b "$base_vcf" \
        -c "$comp_vcf" \
        -o "$output_dir/${software}" \
        -f "$reference_genome" \
        --sizemax 1000000 \
        -r 1000 -p 0 -P 0.5 -O 0.5 \
        --passonly > "$log_file" 2>&1; then
        log "ERROR: truvari failed for $bamfile - $software (exit code: $?)"
        log "See log file for details: $log_file"
        return 1
    fi
    
    log "SUCCESS: Processed $bamfile - $software"
    return 0
}

process_software() {
    local software="$1"
    
    # Read sample list
    mapfile -t bamfiles < "$bamfile_list"
    
    for bamfile in "${bamfiles[@]}"; do
        local base_vcf="${input_paths}/${bamfile}/${bamfile}.base.vcf.gz"
        local comp_vcf="${input_paths}/${bamfile}/${bamfile}.${software}.vcf.gz"
        local output_dir="${output_path}/${bamfile}/${software}"
        
        check_vcf "$base_vcf" "base" || continue
        check_vcf "$comp_vcf" "comparison" || continue
        
        prepare_output_dir "$output_dir" || continue
        
        run_truvari "$base_vcf" "$comp_vcf" "$output_dir" "$bamfile" "$software" || continue
    done
}

# ========== Main Process ==========
log "Starting processing pipeline with software jobs"

# Create main output directory
mkdir -p "$output_path" || { log "ERROR: Failed to create main output directory"; exit 1; }

# Process each software in parallel
for software in "${softwares[@]}"; do
    process_software "$software" &
done

wait

log "All software processing completed successfully!"