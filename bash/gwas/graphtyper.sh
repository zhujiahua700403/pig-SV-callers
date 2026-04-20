#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

# Input parameters
input_base="${BASE_DIR}/result/WGS-SVs/gwas/manta"
bam_path="${BASE_DIR}/LGS/bam/WGS"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
software_path="${BASE_DIR}/software"
bam_list="${BASE_DIR}/result/WGS-SVs/gwas/list/all_bam_list"
bamfile_list="${BASE_DIR}/result/WGS-SVs/gwas/1.list"
threads=18
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"

# Output paths
graphtyper_output="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper"
final_output_path="${BASE_DIR}/result/WGS-SVs/gwas/graphtyper"
region_file="${BASE_DIR}/result/WGS-SVs/graphtyper/regions.txt"
svimmer_output="${BASE_DIR}/result/WGS-SVs/gwas/merge/svimmer_merged"
vcf_list_dir="${BASE_DIR}/result/WGS-SVs/gwas/list"
temp_dir="${BASE_DIR}/result/WGS-SVs/gwas/merge/temp"
genotyped_vcfs_dir="${BASE_DIR}/result/WGS-SVs/gwas/merge/genotyped_vcfs"

# ========== Initialization ==========
# Create all necessary directories
mkdir -p "$graphtyper_output" "$final_output_path" "$svimmer_output" "$vcf_list_dir" \
         "${graphtyper_output}/merged" "$temp_dir" "$genotyped_vcfs_dir" || {
    echo "Error: Failed to create output directories"
    exit 1
}

# Check dependencies
check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed."
        exit 1
    fi
}

for tool in bcftools tabix bgzip python; do
    check_dependency "$tool"
done

# ========== Phase 1: Individual sample genotyping ==========
echo "[$(date)] Starting individual sample genotyping pipeline"

# 1. Generate BAM file list (based on gwas.list)
echo "[$(date)] Generating BAM file list based on gwas.list..."
> "$bam_list"  # Clear file

# Check if bamfile list exists
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfile list not found at $bamfile_list"
    exit 1
fi

# Declare associative array
declare -A bamfile_to_bam_path

# Read bamfile list and find matching BAM files
readarray -t bamfiles < "$bamfile_list"

for bamfile in "${bamfiles[@]}"; do
    # Remove possible whitespace
    bamfile=$(echo "$bamfile" | xargs)
    
    # Skip empty lines
    if [ -z "$bamfile" ]; then
        continue
    fi
    
    echo "Looking for BAM files matching pattern: $bamfile.*.bam"
    
    # Find matching BAM files (format: $bamfile.*.bam)
    found_bams=$(find "$bam_path" -name "${bamfile}.*.bam" -type f)
    
    if [[ -n "$found_bams" ]]; then
        # If multiple matches, take the first one
        first_bam=$(echo "$found_bams" | head -n 1)
        echo "$first_bam" >> "$bam_list"
        bamfile_to_bam_path["$bamfile"]="$first_bam"
        echo "Found BAM for $bamfile: $first_bam"
    else
        echo "Warning: No BAM file found matching pattern $bamfile.*.bam"
        # Try other possible BAM file patterns
        alt_bams=$(find "$bam_path" -name "*${bamfile}*.bam" -type f)
        if [[ -n "$alt_bams" ]]; then
            first_alt_bam=$(echo "$alt_bams" | head -n 1)
            echo "$first_alt_bam" >> "$bam_list"
            bamfile_to_bam_path["$bamfile"]="$first_alt_bam"
            echo "Found alternative BAM for $bamfile: $first_alt_bam"
        else
            echo "Error: No BAM file found for $bamfile with any pattern"
        fi
    fi
done

# Check BAM list file
if [[ ! -s "$bam_list" ]]; then
    echo "Error: No BAM files found in the list"
    exit 1
fi

# 2. Individual sample genotyping
echo "[$(date)] Starting individual sample genotyping..."
genotyped_vcfs_list="${vcf_list_dir}/genotyped_vcfs.list"
> "$genotyped_vcfs_list"

# Read bamfile list and process each sample
readarray -t bamfiles < "$bamfile_list"
for bamfile in "${bamfiles[@]}"; do
    # Remove possible whitespace
    bamfile=$(echo "$bamfile" | xargs)
    
    # Skip empty lines
    if [ -z "$bamfile" ]; then
        continue
    fi
    
    echo "[$(date)] Processing sample: $bamfile"
    
    # Possible VCF file paths
    vcf_files=(
        "${input_base}/${bamfile}/${bamfile}.vcf.gz"
        "${input_base}/${bamfile}/${bamfile}.deduped.vcf.gz"
        "${input_base}/${bamfile}/${bamfile}.rmdup.vcf.gz"
    )
    
    # Find the first existing VCF file
    found_vcf=""
    for vcf in "${vcf_files[@]}"; do
        if [[ -f "$vcf" ]]; then
            found_vcf="$vcf"
            break
        fi
    done
    
    if [[ -z "$found_vcf" ]]; then
        echo "Warning: No VCF file found for $bamfile, skipping"
        continue
    fi
    
    # Get corresponding BAM file path (with safety check)
    if [[ -z "${bamfile_to_bam_path[$bamfile]+exists}" ]]; then
        echo "Warning: No BAM file mapped for $bamfile, skipping"
        continue
    fi
    
    bam_path_for_sample="${bamfile_to_bam_path[$bamfile]}"
    if [[ -z "$bam_path_for_sample" ]]; then
        echo "Warning: BAM path is empty for $bamfile, skipping"
        continue
    fi
    
    # Create sample-specific BAM list file (containing only the current sample)
    sample_bam_list="${temp_dir}/${bamfile}_bam.list"
    echo "$bam_path_for_sample" > "$sample_bam_list"
    
    # Filter VCF to keep only sites with POS > 0
    filtered_vcf="${temp_dir}/${bamfile}_filtered.vcf.gz"
    echo "Filtering VCF: $found_vcf -> $filtered_vcf"
    
    bcftools view -i 'POS > 0' "$found_vcf" | \
    awk 'BEGIN {OFS="\t"} 
         /^#/ {print; next} 
         $1 != "" && $1 != " " {print}' | \
    bgzip -c > "$filtered_vcf" || {
        echo "Error: Failed to filter VCF $found_vcf"
        continue
    }
    
    tabix -p vcf "$filtered_vcf" || {
        echo "Error: Failed to index filtered VCF $filtered_vcf"
        continue
    }
    
    # Create output directory for current sample
    sample_output_dir="${genotyped_vcfs_dir}/${bamfile}"
    mkdir -p "$sample_output_dir"
    
    # Run graphtyper genotype_sv for current sample
    echo "[$(date)] Running graphtyper genotype_sv for sample: $bamfile"
    
    if [[ -f "$region_file" ]]; then
        "$software_path"/graphtyper genotype_sv \
            "$reference_genome" \
            "$filtered_vcf" \
            --sams="$sample_bam_list" \
            --region_file="$region_file" \
            --threads=4 \
            --output="$sample_output_dir" || {
            echo "Error: graphtyper genotype_sv failed for sample $bamfile"
            continue
        }
    else
        "$software_path"/graphtyper genotype_sv \
            "$reference_genome" \
            "$filtered_vcf" \
            --sams="$sample_bam_list" \
            --threads=4 \
            --output="$sample_output_dir" || {
            echo "Error: graphtyper genotype_sv failed for sample $bamfile"
            continue
        }
    fi
    
    # Merge chromosome VCFs for this sample
    sample_merged_vcf="${genotyped_vcfs_dir}/${bamfile}.genotyped.vcf.gz"
    sample_chrom_vcf_list="${temp_dir}/${bamfile}_chrom_vcfs.list"
    
    # Find all chromosome VCFs for this sample
    find "$sample_output_dir" -name "*.vcf.gz" | sort > "$sample_chrom_vcf_list"
    
    if [[ -s "$sample_chrom_vcf_list" ]]; then
        bcftools concat \
            --naive \
            --file-list "$sample_chrom_vcf_list" \
            -Oz -o "$sample_merged_vcf" || {
            echo "Error: Failed to merge chromosome VCFs for sample $bamfile"
            continue
        }
        
        tabix -p vcf "$sample_merged_vcf" || {
            echo "Error: Failed to index merged VCF for sample $bamfile"
            continue
        }
        
        # Add genotyped VCF to list
        echo "$sample_merged_vcf" >> "$genotyped_vcfs_list"
        echo "Successfully genotyped and merged sample: $bamfile"
    else
        echo "Warning: No VCF files found for sample $bamfile after genotyping"
    fi
    
    # Clean sample temporary files
    rm -f "$sample_bam_list" "$sample_chrom_vcf_list"
done

# Check genotyped VCF list file
if [[ ! -s "$genotyped_vcfs_list" ]]; then
    echo "Error: No genotyped VCF files found"
    exit 1
fi

# ========== Phase 2: Merge genotyped VCFs using svimmer ==========
echo "[$(date)] Starting svimmer merging of genotyped VCFs..."

# Ensure svimmer output directory exists
mkdir -p "$(dirname "$svimmer_output")"

# Build chromosome arguments for svimmer
chromosome_args=""
for chrom in {1..18}; do
    chromosome_args="$chromosome_args $chrom"
done

# Merge all genotyped VCFs using svimmer
merged_vcf="${svimmer_output}/all.genotyped.merged.vcf.gz"

echo "[$(date)] Merging genotyped VCFs using svimmer..."
python3 "${software_path}/svimmer-0.1/svimmer" \
    "$genotyped_vcfs_list" \
    $chromosome_args \
    --threads "$threads" \
    --ids \
    | bgzip -c > "$merged_vcf" || {
    echo "Error: svimmer merging of genotyped VCFs failed"
    exit 1
}

tabix -p vcf "$merged_vcf" || {
    echo "Error: Failed to index merged genotyped VCF"
    exit 1
}

# ========== Phase 3: Final filtering and processing ==========
echo "[$(date)] Starting final filtering and processing..."

# 1. Sort merged VCF
echo "[$(date)] Sorting merged VCF..."
sorted_vcf="${svimmer_output}/all.genotyped.merged.sorted.vcf.gz"
bcftools sort \
    -Oz -o "$sorted_vcf" \
    "$merged_vcf" || {
    echo "Error: Failed to sort merged VCF"
    exit 1
}
tabix -p vcf "$sorted_vcf"

# 2. Final filtering step
echo "[$(date)] Filtering final VCF..."

# Build complex filter expression combining all conditions
FILTER_EXPR='(SVTYPE="INS" || SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="DUP") && '
FILTER_EXPR+='((SVTYPE="DEL" && QD > 12 && (ABHet > 0.30 || ABHet < 0) && (AC / NUM_MERGED_SVS) < 25 && PASS_AC > 0 && PASS_ratio > 0.1) || '
FILTER_EXPR+='(SVTYPE="INS" && PASS_AC > 0 && (AC / NUM_MERGED_SVS) < 25 && PASS_ratio > 0.1 && (ABHet > 0.25 || ABHet < 0) && MaxAAS > 4)) && '
FILTER_EXPR+='FILTER="PASS"'

# Add SV size filter
FILTER_EXPR+=' && ((SVTYPE="INS" && SVLEN <= 1000000) || '
FILTER_EXPR+='(SVTYPE="DEL" && SVLEN >= -1000000))'

echo "Applying complex filter expression: $FILTER_EXPR"

# Apply all filters in one step
final_vcf="${final_output_path}/graphtyper.all.genotyped.merged.filtered.vcf.gz"
bcftools view -i "$FILTER_EXPR" "$sorted_vcf" -Oz -o "$final_vcf"
if [ $? -ne 0 ]; then
    echo "Error: Failed to apply complex filtering expression."
    exit 1
fi

# Create index
echo "[$(date)] Indexing final VCF..."
bcftools index -t -f "$final_vcf"

# ========== Post-processing ==========
echo "[$(date)] Starting post-processing pipeline"

# 1. Create reference allele file
ref_allele_file="${final_output_path}/ref_alleles.txt"
echo "[$(date)] Creating reference allele file: $ref_allele_file"
zcat "$final_vcf" | \
awk '!/^#/ && $4 != "." {print $3"\t"$4}' > "$ref_allele_file" || {
    echo "Error: Failed to create reference allele file"
    exit 1
}

# 2. Process VCF: remove '+' characters from ALT fields
processed_vcf="${final_output_path}/graphtyper.all.genotyped.merged.filtered.processed.vcf.gz"
echo "[$(date)] Processing VCF file and removing '+' in ALT fields"
zcat "$final_vcf" | \
awk 'BEGIN {OFS="\t"} 
     /^#/ {print; next} 
     {gsub(/\+/, "", $5); print}' | \
bgzip -c > "$processed_vcf" || {
    echo "Error: Failed to process VCF file"
    exit 1
}

# Create index for processed VCF
tabix -p vcf "$processed_vcf" || {
    echo "Error: Failed to index processed VCF"
    exit 1
}

# 3. Replace original file
echo "[$(date)] Replacing original VCF with processed version"
mv "$processed_vcf" "$final_vcf"
mv "${processed_vcf}.tbi" "${final_vcf}.tbi"

# ========== Cleanup intermediate files ==========
echo "[$(date)] Cleaning up intermediate files..."

# Delete temporary directory
if [[ -d "$temp_dir" ]]; then
    rm -rf "$temp_dir"
    echo "Removed temporary directory: $temp_dir"
fi

# Delete genotyped individual sample VCF directory (optional, keep if needed)
if [[ -d "$genotyped_vcfs_dir" ]]; then
    echo "Note: Genotyped individual VCFs are preserved in: $genotyped_vcfs_dir"
    # To clean, uncomment the following line:
    # rm -rf "$genotyped_vcfs_dir"
fi

# Delete VCF list file
if [[ -f "$genotyped_vcfs_list" ]]; then
    rm -f "$genotyped_vcfs_list"
    echo "Removed genotyped VCF list file: $genotyped_vcfs_list"
fi

# Delete BAM list file
if [[ -f "$bam_list" ]]; then
    rm -f "$bam_list"
    echo "Removed BAM list file: $bam_list"
fi

echo "[$(date)] Pipeline completed successfully."
echo "Final merged and filtered VCF: $final_vcf"
echo "Individual genotyped VCFs are preserved in: $genotyped_vcfs_dir"