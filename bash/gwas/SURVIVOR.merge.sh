#!/bin/bash
set -euo pipefail  # Enable strict error handling

# Activate SURVIVOR environment
source activate SURVIVOR

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_dir="${BASE_DIR}/base/genotyped"               # Input VCF file directory
output_dir="${BASE_DIR}/result/WGS-SVs/geno"         # Output directory
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"    # Bamfile list file
list_file="$output_dir/geno.list"                    # Sample list file
survivor_params="1000 1 1 1 0 50"                    # SURVIVOR merge parameters
merged_vcf="$output_dir/all_sample.merged.vcf"       # Merged VCF file
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"  # Chromosomes to keep

# ========== Initialize directory ==========
mkdir -p "$output_dir"

# ========== Read bamfile list ==========
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi

readarray -t bamfiles < "$bamfile_list"

# ========== Create sample list ==========
> "$list_file"  # Clear list file

for bamfile in "${bamfiles[@]}"; do
    # Remove possible whitespace
    bamfile=$(echo "$bamfile" | xargs)
    [[ -z "$bamfile" ]] && continue
    
    # Extract sample ID from bamfile name (assuming format sample.pn.bam)
    sample_id=$(basename "$bamfile" | cut -d'.' -f1)
    
    # Find corresponding vcf.gz file
    vcf_gz="${input_dir}/${sample_id}.genotyped.vcf.gz"
    
    if [ -f "$vcf_gz" ]; then
        # Decompress to temporary file
        gunzip -c "$vcf_gz" > "$output_dir/${sample_id}.temp.vcf"
        
        # Rename sample
        echo "$sample_id" > "$output_dir/sample_name.txt"
        bcftools reheader -s "$output_dir/sample_name.txt" \
            "$output_dir/${sample_id}.temp.vcf" \
            -o "$output_dir/${sample_id}.vcf"
        
        # Clean temporary files
        rm -f "$output_dir/${sample_id}.temp.vcf" "$output_dir/sample_name.txt"
        
        # Add to merge list
        echo "$output_dir/${sample_id}.vcf" >> "$list_file"
    else
        echo "Warning: $vcf_gz not found"
    fi
done

# Check if sample list has content
if [ ! -s "$list_file" ]; then
    echo "Error: No valid VCF files found"
    exit 1
fi

# ========== Merge VCF files ==========
SURVIVOR merge "$list_file" $survivor_params "$merged_vcf"

# ========== Process merged VCF ==========
# 1. Sort and validate using bcftools
sorted_vcf="$output_dir/all_sample.sorted.vcf"
filtered_vcf="$output_dir/all_sample.filtered.vcf"

# Sort then filter invalid records
bcftools sort "$merged_vcf" | \
awk 'BEGIN {FS=OFS="\t"} 
     /^#/ {print; next} 
     $1 ~ /^([0-9]+|X|Y|MT)$/ && $2 > 0 && $2 < 4294967295 {print}' \
> "$filtered_vcf"
bgzip -f "$filtered_vcf"
tabix -f "$filtered_vcf.gz"

# 2. Filter by specified chromosomes and compress
compressed_vcf="$output_dir/all_sample.chrT.merged.vcf.gz"
bcftools view -r "$chromosomes" "$filtered_vcf.gz" -Oz -o "$compressed_vcf"
tabix -f "$compressed_vcf"

# ========== Clean temporary files ==========
rm -f "$merged_vcf" "$sorted_vcf" "$output_dir"/*.vcf  # Keep compressed .vcf.gz files

echo "All samples processed!"
echo "Final merged VCF file: $compressed_vcf"