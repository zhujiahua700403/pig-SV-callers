#!/bin/bash
set -euo pipefail  # Enable strict error handling

# Activate conda environment
source activate truvari

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_paths="${BASE_DIR}/result/WGS-SVs/region/vcf/snp_indel"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari/region/snp_indel/result"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
numbers=("1" "2" "3" "4" "5" "ge6")

# ========== Initialization Checks ==========
[[ -f "$bamfile_list" ]] || { echo "ERROR: bamfiles.list not found" >&2; exit 1; }
[[ -f "$reference_genome" ]] || { echo "ERROR: Reference genome not found" >&2; exit 1; }
[[ -d "$input_paths" ]] || { echo "ERROR: Input directory not found" >&2; exit 1; }

# ========== Function Definitions ==========
check_vcf() {
    local vcf_file="$1"
    [[ -f "$vcf_file" ]] && [[ -f "${vcf_file}.tbi" ]] || { tabix -p vcf "$vcf_file"; }
}

process_software() {
    local software="$1"
    mapfile -t bamfiles < "$bamfile_list"

    for bamfile in "${bamfiles[@]}"; do
        local output_dir="$output_path/$software/$bamfile"
        mkdir -p "$output_dir"

        for number in "${numbers[@]}"; do
            local base_vcf="$input_paths/$bamfile/base/sv_breakpoints_$number.vcf.gz"
            local comp_vcf="$input_paths/$bamfile/$software/sv_breakpoints_$number.vcf.gz"
                        
            rm -rf "$output_dir/$number"

            check_vcf "$base_vcf" || continue
            check_vcf "$comp_vcf" || continue
            
            truvari bench \
                -b "$base_vcf" \
                -c "$comp_vcf" \
                -o "$output_dir/$number" \
                -f "$reference_genome" \
                -r 1000 -p 0 -P 0.5 -O 0.5 --passonly

            echo "Processed $software - $bamfile - $number"
        done
    done
}

# ========== Main Process ==========
mkdir -p "$output_path"

for software in "${softwares[@]}"; do
    process_software "$software"
done

wait

echo "All software processing completed successfully!"