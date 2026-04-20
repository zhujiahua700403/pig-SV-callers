#!/bin/bash
set -euo pipefail  # Enable strict error handling

# Activate conda environment
source activate truvari

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_paths="${BASE_DIR}/result/WGS-SVs/Truvari"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari/size/result"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham")
sizemins=(50 500 1000 10000 1000000)  # Minimum size bins
sizemaxs=(499 999 9999 999999 -1)     # Maximum size bins (-1 means no upper limit)
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
        log "Cleaning existing directory: $dir"
        rm -rf "$dir" || { log "ERROR: Failed to clean directory"; return 1; }
    fi
    
    mkdir -p "$dir" || { log "ERROR: Failed to create output directory"; return 1; }
    return 0
}

run_truvari_size() {
    local base_vcf="$1"
    local comp_vcf="$2"
    local output_dir="$3"
    local bamfile="$4"
    local software="$5"
    local sizemin="$6"
    local sizemax="$7"
    
    local size_label
    if [ "$sizemax" -eq -1 ]; then
        size_label="${sizemin}+"
        log "Running truvari for $bamfile - $software (Size: ¡Ý${sizemin}bp)"
    else
        size_label="${sizemin}_${sizemax}"
        log "Running truvari for $bamfile - $software (Size: ${sizemin}-${sizemax}bp)"
    fi
    
    local size_dir="${output_dir}/${size_label}"
    prepare_output_dir "$size_dir" || return 1
    rm -rf "$size_dir"
    local log_file="${output_dir}/truvari.log"
    
    # Build truvari command
    local truvari_cmd="truvari bench \
        -b \"$base_vcf\" \
        -c \"$comp_vcf\" \
        -o \"$size_dir\" \
        -f \"$reference_genome\" \
        --sizemin \"$sizemin\" \
        -r 1000 -p 0 -P 0.5 -O 0.5 \
        --passonly"
    
    # Add sizemax parameter if not -1
    if [ "$sizemax" -ne -1 ]; then
        truvari_cmd+=" --sizemax \"$sizemax\""
    fi
    
    # Execute command
    if ! eval "$truvari_cmd" > "$log_file" 2>&1; then
        log "ERROR: truvari failed for $bamfile - $software (Size: ${size_label})"
        log "See log file for details: $log_file"
        return 1
    fi
    
    log "SUCCESS: Processed $bamfile - $software (Size: ${size_label})"
    return 0
}

process_sample_size() {
    local bamfile="$1"
    local software="$2"
    local sizemin="$3"
    local sizemax="$4"
    
    local base_vcf="${input_paths}/${bamfile}/${bamfile}.base.vcf.gz"
    local comp_vcf="${input_paths}/${bamfile}/${bamfile}.${software}.vcf.gz"
    local output_dir="${output_path}/${bamfile}/${software}"
    
    check_vcf "$base_vcf" "base" || return 1
    check_vcf "$comp_vcf" "comparison" || return 1
    
    run_truvari_size "$base_vcf" "$comp_vcf" "$output_dir" "$bamfile" "$software" "$sizemin" "$sizemax"
}

# ========== Main Process ==========
log "Starting size-based SV analysis pipeline..."

mapfile -t bamfiles < "$bamfile_list"
log "Found ${#bamfiles[@]} samples to process"

for bamfile in "${bamfiles[@]}"; do
    for software in "${softwares[@]}"; do
        # Process each size bin
        for i in "${!sizemins[@]}"; do
            sizemin="${sizemins[i]}"
            if [ $i -lt ${#sizemaxs[@]} ]; then
                sizemax="${sizemaxs[i]}"
            else
                sizemax="-1"  # No upper limit for the last bin
            fi
            
            # Run in parallel for each software combination
            process_sample_size "$bamfile" "$software" "$sizemin" "$sizemax" &
        done
    done
    wait  # Wait for all software bins of current sample to finish
done

wait  # Wait for all processes to complete
log "All size-based SV analysis completed successfully!"