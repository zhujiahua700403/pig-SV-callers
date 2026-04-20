#!/bin/bash

# Activate conda environment
source activate truvari

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_paths="${BASE_DIR}/result/WGS-SVs/region/vcf/trf"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari/region/trf/result"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# ========== Initialization Checks ==========
if [ ! -f "$bamfile_list" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: bamfiles.list not found at $bamfile_list" >&2
    exit 1
fi

if [ ! -f "$reference_genome" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Reference genome not found at $reference_genome" >&2
    exit 1
fi

if [ ! -d "$input_paths" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Input directory not found at $input_paths" >&2
    exit 1
fi

# ========== Function Definitions ==========
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_and_index_vcf() {
    local vcf_file="$1"
    local vcf_type="$2"
    
    if [ ! -f "$vcf_file" ]; then
        log "ERROR: $vcf_type VCF not found: $vcf_file"
        return 1
    fi
    
    if [ ! -f "${vcf_file}.tbi" ]; then
        log "Indexing $vcf_type VCF: $vcf_file"
        if ! tabix -p vcf "$vcf_file"; then
            log "ERROR: Failed to index $vcf_type VCF: $vcf_file"
            return 1
        fi
    fi
    
    return 0
}

process_bamfile() {
    local bamfile="$1"
    local software="$2"
    
    local base_vcf="${input_paths}/${bamfile}/base/overlap/${bamfile}_base_trf_overlap.vcf.gz"
    local comp_vcf="${input_paths}/${bamfile}/${software}/overlap/${bamfile}_${software}_trf_overlap.vcf.gz"
    local output_dir="${output_path}/${bamfile}/${software}"
    
    mkdir -p "${output_dir}"
    rm -rf "${output_dir}/${software}/overlap"

    check_and_index_vcf "$base_vcf" "base" || return 1
    check_and_index_vcf "$comp_vcf" "comparison" || return 1
    
    log "Running truvari bench for $bamfile - $software"
    truvari bench \
        -b "$base_vcf" \
        -c "$comp_vcf" \
        -o "$output_dir/${software}/overlap" \
        -f "$reference_genome" \
        -r 1000 -p 0 -P 0.5 -O 0.5 \
        --passonly
    
    local ret=$?
    if [ $ret -ne 0 ]; then
        log "ERROR: truvari bench failed for $bamfile - $software (exit code: $ret)"
        return $ret
    fi
    
    log "Successfully processed $bamfile - $software"
    return 0
}

# Second process_bamfile function (overwrites the first one)
process_bamfile() {
    local bamfile="$1"
    local software="$2"
    
    local base_vcf="${input_paths}/${bamfile}/base/no_overlap/${bamfile}_base_no_trf_overlap.vcf.gz"
    local comp_vcf="${input_paths}/${bamfile}/${software}/no_overlap/${bamfile}_${software}_no_trf_overlap.vcf.gz"
    local output_dir="${output_path}/${bamfile}/${software}"
    
    mkdir -p "${output_dir}"
    rm -rf "${output_dir}/no_overlap"

    check_and_index_vcf "$base_vcf" "base" || return 1
    check_and_index_vcf "$comp_vcf" "comparison" || return 1
    
    log "Running truvari bench for $bamfile - $software"
    truvari bench \
        -b "$base_vcf" \
        -c "$comp_vcf" \
        -o "$output_dir/no_overlap" \
        -f "$reference_genome" \
        -r 1000 -p 0 -P 0.5 -O 0.5 \
        --passonly
    
    local ret=$?
    if [ $ret -ne 0 ]; then
        log "ERROR: truvari bench failed for $bamfile - $software (exit code: $ret)"
        return $ret
    fi
    
    log "Successfully processed $bamfile - $software"
    return 0
}

# ========== Main Process ==========
log "Starting processing pipeline..."

readarray -t bamfiles < "$bamfile_list"
log "Found ${#bamfiles[@]} samples to process"

# Use parallel control
for software in "${softwares[@]}"; do
    for bamfile in "${bamfiles[@]}"; do
        process_bamfile "$bamfile" "$software"
    done
    # Wait for all tasks of the current software to finish
    wait
done

log "All BAM processing completed successfully!"