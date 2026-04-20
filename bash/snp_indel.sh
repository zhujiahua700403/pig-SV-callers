#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

sv_path="${BASE_DIR}/result/WGS-SVs/Truvari/INS_DEL"
snp_indel_path="${BASE_DIR}/gvcf"
output_dir="${BASE_DIR}/result/WGS-SVs/region/result"
vcf_output_dir="${BASE_DIR}/result/WGS-SVs/region/vcf/snp_indel"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham" "base")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list" >&2
    exit 1
fi
mapfile -t bamfiles < "$bamfile_list"

# ========== Find compatible VCF file ==========
find_compatible_vcf() {
    local base_dir="$1"
    local sample_name="$2"
    
    # Try three naming patterns
    local vcf_patterns=(
        "${base_dir}/${sample_name}.noMT.h2.g.vcf.gz"
        "${base_dir}/${sample_name}.NOXYMT.h2.g.vcf.gz"
        "${base_dir}/${sample_name}.h2.g.vcf.gz"
    )
    
    for pattern in "${vcf_patterns[@]}"; do
        if [ -f "$pattern" ]; then
            echo "$pattern"
            return 0
        fi
    done
    
    echo "Error: Cannot find compatible VCF file (sample: $sample_name)" >&2
    return 1
}

# ========== Generate BED file function ==========
generate_bed() {
    local input_vcf=$1
    local output_bed=$2
    local variant_type=$3
    local tmp_bed=$(mktemp)
    
    echo "Generating BED file for ${variant_type}: $output_bed"
    
    # Check if input file exists and is readable
    if [ ! -r "$input_vcf" ]; then
        echo "Error: Cannot read VCF file $input_vcf" >&2
        return 1
    fi

    # Process different variant types
    if [[ "$variant_type" == "SV" ]]; then
        if ! bcftools view "$input_vcf" 2>/dev/null | \
            awk 'BEGIN {OFS="\t"} !/^#/ {if($2 > 0) print $1, ($2-1), $2, $3}' > "$tmp_bed"; then
            echo "Error: Failed to process SV VCF file" >&2
            rm -f "$tmp_bed"
            return 1
        fi
    else
        if ! bcftools view -v "$variant_type" "$input_vcf" 2>/dev/null | \
            awk 'BEGIN {OFS="\t"} !/^#/ {if($2 > 0) print $1, ($2-1), $2}' > "$tmp_bed"; then
            echo "Error: Failed to process ${variant_type} VCF file" >&2
            rm -f "$tmp_bed"
            return 1
        fi
    fi

    # Check if output is valid
    if [ ! -s "$tmp_bed" ]; then
        echo "Warning: Generated ${variant_type} BED file is empty" >&2
        rm -f "$tmp_bed"
        return 1
    fi

    # Ensure directory exists
    mkdir -p "$(dirname "$output_bed")"
    mv "$tmp_bed" "$output_bed" || {
        echo "Error: Failed to move temporary file to $output_bed" >&2
        return 1
    }

    echo "Successfully generated ${variant_type} BED file with $(wc -l < "$output_bed") records"
    return 0
}

# ========== Extract VCF entries matching BED ==========
extract_vcf() {
    local input_vcf=$1
    local input_bed=$2
    local output_vcf=$3
    local tmp_bed=$(mktemp)
    
    echo "Extracting VCF entries matching BED: $output_vcf"
    
    # Check input files
    if [ ! -r "$input_vcf" ] || [ ! -r "$input_bed" ]; then
        echo "Error: Input files not readable" >&2
        return 1
    fi

    # Convert BED format to bcftools-compatible format
    awk 'BEGIN {OFS="\t"} {print $1, $2+1, $3}' "$input_bed" > "$tmp_bed" || {
        echo "Error: BED file conversion failed" >&2
        rm -f "$tmp_bed"
        return 1
    }

    # Ensure output directory exists
    mkdir -p "$(dirname "$output_vcf")"

    # Extract VCF entries
    if ! bcftools view -R "$tmp_bed" "$input_vcf" -Oz -o "$output_vcf"; then
        echo "Error: Failed to extract VCF entries" >&2
        rm -f "$tmp_bed" "$output_vcf"
        return 1
    fi

    # Index VCF file
    if ! bcftools index -t "$output_vcf"; then
        echo "Error: Failed to index VCF file" >&2
        rm -f "$output_vcf"
        return 1
    fi

    rm -f "$tmp_bed"
    echo "Successfully generated VCF: $output_vcf ($(zgrep -vc '^#' "$output_vcf" || echo 0) variants)"
    return 0
}

# ========== Process single sample function ==========
process_sample() {
    local bamfile=$1
    local software=$2
    
    echo "Starting processing $bamfile - $software"
    
    # Find compatible VCF file
    snp_vcf=$(find_compatible_vcf "$snp_indel_path" "$bamfile") || return 1
    indel_vcf="$snp_vcf"  # SNP and Indel use the same file
    sv_vcf="$sv_path/$bamfile/$bamfile.$software.vcf.gz"
    
    # Check if SV VCF exists
    if [ ! -f "$sv_vcf" ]; then
        echo "Error: SV VCF file does not exist: $sv_vcf" >&2
        return 1
    fi

    # Create output directories
    mkdir -p "$output_dir/$bamfile/$software" "$vcf_output_dir/$bamfile/$software" || {
        echo "Error: Failed to create output directories" >&2
        return 1
    }

    # 1. Generate BED files for each variant type
    generate_bed "$sv_vcf" "$output_dir/$bamfile/$software/sv.bed" "SV" || return 1
    generate_bed "$snp_vcf" "$output_dir/$bamfile/$software/snps.bed" "snps" || return 1
    generate_bed "$indel_vcf" "$output_dir/$bamfile/$software/indels.bed" "indels" || return 1

    # 2. Process SV breakpoints
    echo "Processing SV breakpoints and extending by 50bp..."
    awk 'BEGIN {OFS="\t"} {
        chr = $1;
        start = int($2);
        end = int($3);
        sv_id = $4;
        
        start_ext = (start - 50) > 0 ? (start - 50) : 0;
        end_ext = end + 50;
        
        print chr, start_ext, start, sv_id "_start";
        print chr, end, end_ext, sv_id "_end";
    }' "$output_dir/$bamfile/$software/sv.bed" > "$output_dir/$bamfile/$software/sv_breakpoints_extended.bed" || {
        echo "Error: Failed to generate extended breakpoint BED" >&2
        return 1
    }

    # 3. Merge SNP and Indel BED files
    echo "Merging SNP and Indel BED files..."
    cat "$output_dir/$bamfile/$software/snps.bed" "$output_dir/$bamfile/$software/indels.bed" > "$output_dir/$bamfile/$software/snp_indel_combined.bed" || {
        echo "Error: Failed to merge BED files" >&2
        return 1
    }
    
    # 4. Compare with merged SNP/Indel file
    echo "Comparing with merged SNP/Indel file..."
    bedtools intersect -a "$output_dir/$bamfile/$software/sv_breakpoints_extended.bed" \
                       -b "$output_dir/$bamfile/$software/snp_indel_combined.bed" -c > \
                       "$output_dir/$bamfile/$software/sv_breakpoints_snp_indel_count.bed" || {
        echo "Error: bedtools intersect failed" >&2
        return 1
    }

    # 5. Add classification information
    echo "Adding classification information..."
    awk 'BEGIN {OFS="\t"} {
        total = $5;
        category = (total <= 5) ? total : "¡Ý6";
        print $1, $2, $3, $4, $5, total, category;
    }' "$output_dir/$bamfile/$software/sv_breakpoints_snp_indel_count.bed" > \
       "$output_dir/$bamfile/$software/sv_breakpoints_snp_indel_count_classified.bed" || {
        echo "Error: Failed to add classification information" >&2
        return 1
    }

    # 6. Output by count category and generate VCF
    echo "Outputting by variant count category and generating VCF..."
    local categories=("1" "2" "3" "4" "5" "¡Ý6")
    for category in "${categories[@]}"; do
        local bed_file output_vcf
        
        if [[ "$category" == "¡Ý6" ]]; then
            bed_file="$output_dir/$bamfile/$software/sv_breakpoints_snp_indel_ge6.bed"
            awk '$7 == "¡Ý6"' "$output_dir/$bamfile/$software/sv_breakpoints_snp_indel_count_classified.bed" > "$bed_file"
            output_vcf="$vcf_output_dir/$bamfile/$software/sv_breakpoints_ge6.vcf.gz"
        else
            bed_file="$output_dir/$bamfile/$software/sv_breakpoints_snp_indel_${category}.bed"
            awk -v cat="$category" '$7 == cat' "$output_dir/$bamfile/$software/sv_breakpoints_snp_indel_count_classified.bed" > "$bed_file"
            output_vcf="$vcf_output_dir/$bamfile/$software/sv_breakpoints_${category}.vcf.gz"
        fi
        
        if [ -s "$bed_file" ]; then
            extract_vcf "$sv_vcf" "$bed_file" "$output_vcf" || {
                echo "Warning: Failed to generate VCF for count=${category}" >&2
                continue
            }
            echo "Successfully generated VCF for count=${category}: $output_vcf"
        else
            echo "Warning: BED file for count=${category} is empty" >&2
        fi
    done
    
    # 7. Generate merged VCF of all count categories
    echo "Generating merged VCF for all count categories..."
    local vcf_list=()
    for category in "${categories[@]}"; do
        local vcf
        if [[ "$category" == "¡Ý6" ]]; then
            vcf="$vcf_output_dir/$bamfile/$software/sv_breakpoints_ge6.vcf.gz"
        else
            vcf="$vcf_output_dir/$bamfile/$software/sv_breakpoints_${category}.vcf.gz"
        fi
        [ -f "$vcf" ] && vcf_list+=("$vcf")
    done
    
    if [ ${#vcf_list[@]} -gt 0 ]; then
        local merged_vcf="$vcf_output_dir/$bamfile/$software/sv_breakpoints_all.vcf.gz"
        if bcftools concat -a "${vcf_list[@]}" -Oz -o "$merged_vcf" && \
           bcftools index -t "$merged_vcf"; then
            echo "Successfully generated merged VCF: $merged_vcf"
        else
            echo "Warning: Failed to merge VCF files" >&2
        fi
    else
        echo "Warning: No VCF files available for merging" >&2
    fi
    
    echo "Completed processing $bamfile - $software"
    return 0
}

# ========== Main pipeline ==========
echo "Starting analysis pipeline..."

# Main processing loop
for software in "${softwares[@]}"; do
    echo "Starting processing for software: $software"
    
    for bamfile in "${bamfiles[@]}"; do
        process_sample "$bamfile" "$software" || {
            echo "Warning: Processing $bamfile - $software failed, continuing to next sample" >&2
            continue
        }
    done
    
    echo "Completed processing for software: $software"
done

echo "All samples and software analysis completed!"
echo "Results saved in: $output_dir and $vcf_output_dir"