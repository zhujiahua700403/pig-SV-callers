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
bamfile_list="${BASE_DIR}/result/WGS-SVs/gwas/gwas.list"
threads=18
chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18"

# Output paths
graphtyper_output="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper"
final_output_path="${BASE_DIR}/result/WGS-SVs/gwas/graphtyper"
region_file="${BASE_DIR}/result/WGS-SVs/graphtyper/regions.txt"
svimmer_output="${BASE_DIR}/result/WGS-SVs/gwas/merge/svimmer_merged"
vcf_list_dir="${BASE_DIR}/result/WGS-SVs/gwas/list"

# ========== Initialization ==========
# Create all necessary directories
mkdir -p "$graphtyper_output" "$final_output_path" "$svimmer_output" "$vcf_list_dir" \
         "${graphtyper_output}/merged" || {  # Ensure merged directory exists
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

# ========== SV merging stage ==========
echo "[$(date)] Starting SV merging pipeline"

# 1. Generate BAM file list
echo "[$(date)] Generating BAM file list..."
find "${bam_path}" -name "*.bam" | sort > "$bam_list" || {
    echo "Error: Failed to generate BAM list"
    exit 1
}

# 2. Generate VCF file list (based on bamfiles in gwas.list)
echo "[$(date)] Generating VCF file list..."
vcf_list="${vcf_list_dir}/vcfs.list"
> "$vcf_list"  # Clear file

# Check if bamfile list exists
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfile list not found at $bamfile_list"
    exit 1
fi

# Read bamfile list and find corresponding VCF files
readarray -t bamfiles < "$bamfile_list"
for bamfile in "${bamfiles[@]}"; do
    # Remove possible whitespace
    bamfile=$(echo "$bamfile" | xargs)
    
    # Skip empty lines
    if [ -z "$bamfile" ]; then
        continue
    fi
    
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
    
    # If VCF file found, add to list
    if [[ -n "$found_vcf" ]]; then
        echo "$found_vcf" >> "$vcf_list"
        echo "Found VCF for $bamfile: $found_vcf"
    else
        echo "Warning: No VCF file found for $bamfile"
    fi
done

# 3. Check VCF list file
if [[ ! -s "$vcf_list" ]]; then
    echo "Error: No VCF files found in the list"
    exit 1
fi

# 4. Merge VCFs using svimmer
echo "[$(date)] Merging VCFs using svimmer..."
merged_vcf="${svimmer_output}/all.merged.vcf.gz"
rm -rf $merged_vcf $merged_vcf.tbi
# Ensure svimmer output directory exists
mkdir -p "$(dirname "$merged_vcf")"

python3 "${software_path}/svimmer-0.1/svimmer" \
    "$vcf_list" \
    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 \
    --threads "$threads" \
    --ids \
    | bgzip -c > "$merged_vcf" || {
    echo "Error: svimmer merging failed"
    exit 1
}

tabix -p vcf "$merged_vcf" || {
    echo "Error: Failed to index merged VCF"
    exit 1
}

# ========== VCF comparison and filtering ==========
echo "[$(date)] Starting VCF comparison and filtering pipeline"

# 1. Define output path
comparison_dir="${BASE_DIR}/result/WGS-SVs/gwas"
mkdir -p "$comparison_dir" || {
    echo "Error: Failed to create comparison directory"
    exit 1
}

# 2. Extract SV information from long-read VCF (CHROM, POS, END, SVTYPE)
lrs_vcf="${BASE_DIR}/result/LSV/all_sample.merge.vcf.gz"
lrs_bed="${comparison_dir}/LRS.bed"

echo "[$(date)] Extracting SV information from long-read VCF: $lrs_vcf"
zcat "$lrs_vcf" | \
awk 'BEGIN {OFS="\t"} 
     !/^#/ {
         n = split($8, info, ";");
         svtype = "";
         end = $2;
         for (i = 1; i <= n; i++) {
             if (info[i] ~ /^SVTYPE=/) {
                 svtype = substr(info[i], 8);
             }
             if (info[i] ~ /^END=/) {
                 end = substr(info[i], 5);
             }
         }
         if (svtype != "") {
             print $1, $2, end, svtype
         }
     }' > "$lrs_bed" || {
    echo "Error: Failed to extract SV info from long-read VCF"
    exit 1
}

# 3. Extract SV information from short-read VCF (CHROM, POS, END, SVTYPE)
srs_bed="${comparison_dir}/SRS.bed"

echo "[$(date)] Extracting SV information from short-read VCF: $merged_vcf"
zcat "$merged_vcf" | \
awk 'BEGIN {OFS="\t"} 
     !/^#/ {
         n = split($8, info, ";");
         svtype = "";
         end = $2;
         for (i = 1; i <= n; i++) {
             if (info[i] ~ /^SVTYPE=/) {
                 svtype = substr(info[i], 8);
             }
             if (info[i] ~ /^END=/) {
                 end = substr(info[i], 5);
             }
         }
         if (svtype != "") {
             print $1, $2, end, svtype
         }
     }' > "$srs_bed" || {
    echo "Error: Failed to extract SV info from short-read VCF"
    exit 1
}

# Fix chromosome ordering
echo "[$(date)] Sorting BED files to ensure consistent chromosome ordering..."

# Sort original BED files
sort -k1,1V -k2,2n "$lrs_bed" > "${lrs_bed}.sorted" && mv "${lrs_bed}.sorted" "$lrs_bed"
sort -k1,1V -k2,2n "$srs_bed" > "${srs_bed}.sorted" && mv "${srs_bed}.sorted" "$srs_bed"

# 4. Compare overlaps with different ratios per SV type
merged_bed="${comparison_dir}/merge.bed"
filtered_vcf="${comparison_dir}/filtered_merged.vcf.gz"

# Initialize merged file
> "$merged_bed"

# Define minimum overlap ratio for each SV type
declare -A min_overlap=(
    ["INS"]=0.90
    ["DEL"]=0.50
    ["INV"]=0.80
    ["DUP"]=0.80
)

# Process each SV type separately
for svtype in INS DEL INV DUP; do
    echo "[$(date)] Processing $svtype variants with min overlap: ${min_overlap[$svtype]}"
    
    # Extract variants of current type
    awk -v t="$svtype" '$4 == t' "$lrs_bed" > "${comparison_dir}/lrs_${svtype}.bed"
    awk -v t="$svtype" '$4 == t' "$srs_bed" > "${comparison_dir}/srs_${svtype}.bed"
    
    # Sort extracted BED files
    sort -k1,1V -k2,2n "${comparison_dir}/lrs_${svtype}.bed" > "${comparison_dir}/lrs_${svtype}.sorted.bed"
    sort -k1,1V -k2,2n "${comparison_dir}/srs_${svtype}.bed" > "${comparison_dir}/srs_${svtype}.sorted.bed"
    
    # Skip empty files
    if [ ! -s "${comparison_dir}/lrs_${svtype}.sorted.bed" ] || [ ! -s "${comparison_dir}/srs_${svtype}.sorted.bed" ]; then
        echo "Warning: No $svtype variants found in one of the files"
        continue
    fi
    
    # Use bedtools for overlap comparison
    bedtools intersect \
        -a "${comparison_dir}/srs_${svtype}.sorted.bed" \
        -b "${comparison_dir}/lrs_${svtype}.sorted.bed" \
        -f ${min_overlap[$svtype]} \
        -wa \
        -sorted \
        >> "$merged_bed" || {
        echo "Error: bedtools intersect failed for $svtype"
        exit 1
    }
done

# 5. Filter VCF based on matched CHROM, POS, END, SVTYPE
echo "[$(date)] Filtering VCF based on matched CHROM,POS,END,SVTYPE"

# Create temporary match file
matched_variants="${comparison_dir}/matched_variants.txt"
awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' "$merged_bed" > "$matched_variants"

# Filter VCF using AWK
zcat "$merged_vcf" | \
awk -v OFS="\t" '
    BEGIN {
        # Load matched variant information into hash table
        while (getline < "'"$matched_variants"'" > 0) {
            key = $1 "|" $2 "|" $3 "|" $4;
            matched[key] = 1;
        }
    }
    /^#/ { print; next }   # Print header
    {
        # Extract coordinates and type of current variant
        chrom = $1;
        pos = $2;
        end = "";
        svtype = "";
        
        # Extract END and SVTYPE from INFO field
        n = split($8, info, ";");
        for (i = 1; i <= n; i++) {
            if (info[i] ~ /^END=/) {
                split(info[i], arr, "=");
                end = arr[2];
            }
            if (info[i] ~ /^SVTYPE=/) {
                split(info[i], arr, "=");
                svtype = arr[2];
            }
        }
        
        # Check if match
        if (end != "" && svtype != "") {
            key = chrom "|" pos "|" end "|" svtype;
            if (key in matched) {
                print;
            }
        }
    }' | bgzip -c > "$filtered_vcf" || {
    echo "Error: Failed to filter VCF by coordinates and type"
    exit 1
}

# 6. Index filtered VCF
echo "[$(date)] Indexing filtered VCF: $filtered_vcf"
tabix -p vcf "$filtered_vcf" || {
    echo "Error: Failed to index filtered VCF"
    exit 1
}

# 7. Clean temporary files
echo "[$(date)] Cleaning temporary comparison files"
rm -f "${comparison_dir}/lrs_"*".bed" "${comparison_dir}/srs_"*".bed" "$matched_variants"

echo "[$(date)] VCF comparison and filtering completed. Final VCF: $filtered_vcf"