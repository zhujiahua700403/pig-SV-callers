#!/bin/bash

# Activate necessary environment
source activate truvari

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_path="${BASE_DIR}/result/WGS-SVs/Truvari"
output_dir="${BASE_DIR}/result/WGS-SVs/Truvari/svtype"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# Define SV types and software list
sv_types=("INS" "DEL" "DUP" "INV" "BND")
softwares=("base" "breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham")

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Create output directory structure
mkdir -p "$output_dir"

# Define function to process a single software
process_software() {
    local software="$1"

    # Loop over samples
    for bamfile in "${bamfiles[@]}"; do
        echo "[$(date)] Processing $bamfile - $software"
        
        # Check if input VCF file exists
        input_vcf="$input_path/$bamfile/$bamfile.$software.vcf.gz"
        if [ ! -f "$input_vcf" ]; then
            echo "Warning: Input file $input_vcf not found, skipping..."
            continue
        fi

        # Create software output directory
        mkdir -p "$output_dir/$bamfile/$software"
        mkdir -p "$output_dir/result/$bamfile/$software"

        # Loop over SV types and perform filtering
        for sv_type in "${sv_types[@]}"; do
            output_vcf="$output_dir/$bamfile/$software/$bamfile.$software.$sv_type.vcf"
            echo "Processing $bamfile - $software - $sv_type..."
            
            # Filter by specific SV type using bcftools
            bcftools filter \
                -i "SVTYPE=\"$sv_type\"" \
                "$input_vcf" \
                -o "$output_vcf"

            # Compress and index the result
            if [ -f "$output_vcf" ]; then
                bgzip -f "$output_vcf"
                tabix -f "$output_vcf.gz"
            else
                echo "Error: Failed to create $output_vcf"
            fi
        done

        # Run truvari bench for each SV type
        for sv_type in "${sv_types[@]}"; do
            base_vcf="$output_dir/$bamfile/base/$bamfile.base.$sv_type.vcf.gz"
            comp_vcf="$output_dir/$bamfile/$software/$bamfile.$software.$sv_type.vcf.gz"
            rm -rf "$output_dir/result/$bamfile/$software/$sv_type" 
            mkdir -p "$output_dir/result/$bamfile/$software"
            if [ -f "$base_vcf" ] && [ -f "$comp_vcf" ]; then
                echo "Running truvari bench for $sv_type..."
                truvari bench \
                    -b "$base_vcf" \
                    -c "$comp_vcf" \
                    -o "$output_dir/result/$bamfile/$software/$sv_type" \
                    -f "$reference_genome" -r 1000 -p 0 -P 0.5 -O 0.5 \
                    --passonly > "$output_dir/$bamfile/$software/${sv_type}_truvari.log" 2>&1
            else
                echo "Warning: Missing VCF files for $sv_type comparison (base: $base_vcf, comp: $comp_vcf)"
            fi
        done
    done
}

# Iterate over all software and process in parallel
for software in "${softwares[@]}"; do
    process_software "$software" &
done

# Wait for all background tasks to finish
wait

echo "[$(date)] All BAM processing completed!"