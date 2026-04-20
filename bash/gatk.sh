#!/bin/bash
set -euo pipefail

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

reference="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
input_dir="${BASE_DIR}/gvcf"                # Directory containing per-sample gVCF files
temp_dir="${BASE_DIR}/gvcf/temp"            # Temporary directory for filtered gVCFs
output_gvcf="${BASE_DIR}/gvcf/all.g.vcf.gz" # Merged gVCF output
final_vcf_dir="${BASE_DIR}/gvcf"            # Directory for final SNP/INDEL VCFs
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
java_opts="-Xmx10g -Djava.io.tmpdir=/tmp"

# ========== Initialization ==========
mkdir -p "$temp_dir" "$final_vcf_dir"

# Check tool availability
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not installed"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "Error: tabix not installed"; exit 1; }

# ========== 1. Prepare per?sample gVCFs (filter to autosomes 1-18) ==========
# Read sample list
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

for bamfile in "${bamfiles[@]}"; do
    bamfile=$(echo "$bamfile" | tr -d '\r\n')
    
    # Possible input file names
    input1="${input_dir}/${bamfile}.noMT.h2.g.vcf.gz"
    input2="${input_dir}/${bamfile}.NOXYMT.h2.g.vcf.gz"
    output="${temp_dir}/${bamfile}.1_18.h2.g.vcf.gz"
    
    if [ -f "$input1" ]; then
        input_file="$input1"
    elif [ -f "$input2" ]; then
        input_file="$input2"
    else
        echo "Warning: No gVCF found for $bamfile (tried $input1 and $input2)"
        continue
    fi
    
    echo "Processing $bamfile: filtering to chromosomes 1-18..."
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 \
        -Oz -o "$output" "$input_file"
    tabix "$output"
    
    if [ ! -f "$output" ]; then
        echo "Error: Failed to process $bamfile"
        exit 1
    fi
done

# ========== 2. Combine gVCFs using GATK ==========
gvcf_list="${temp_dir}/gvcf.list"
find "$temp_dir" -name "*.1_18.h2.g.vcf.gz" | sort > "$gvcf_list"

if [ ! -s "$gvcf_list" ]; then
    echo "Error: No valid gVCF files found in $temp_dir"
    exit 1
fi

echo "Activating GATK environment..."
source activate GATK || { echo "Error: Failed to activate gatk environment"; exit 1; }

echo "Combining GVCFs..."
gatk CombineGVCFs \
    -R "$reference" \
    --variant "$gvcf_list" \
    -O "$output_gvcf"

if [ ! -f "$output_gvcf" ]; then
    echo "Error: Failed to combine GVCFs"
    exit 1
fi
tabix -p vcf "$output_gvcf"
echo "Successfully created merged gVCF: $output_gvcf"

# ========== 3. GenotypeGVCFs (population variant calling) ==========
echo "Running GenotypeGVCFs..."
gatk --java-options "${java_opts}" GenotypeGVCFs \
    -R "${reference}" \
    -V "${output_gvcf}" \
    -O "${final_vcf_dir}/all.merge_raw.vcf" || {
    echo "Error: GenotypeGVCFs failed"
    exit 1
}

# ========== 4. SNP processing ==========
echo "Processing SNPs..."

gatk --java-options "${java_opts}" SelectVariants \
    -R "${reference}" \
    -V "${final_vcf_dir}/all.merge_raw.vcf" \
    --select-type SNP \
    -O "${final_vcf_dir}/all.raw.snp.vcf" || {
    echo "Error: SNP selection failed"
    exit 1
}

gatk --java-options "${java_opts}" VariantFiltration \
    -R "${reference}" \
    -V "${final_vcf_dir}/all.raw.snp.vcf" \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name 'SNP_filter' \
    -O "${final_vcf_dir}/all.filter.snp.vcf" || {
    echo "Error: SNP filtration failed"
    exit 1
}

gatk --java-options "${java_opts}" SelectVariants \
    -R "${reference}" \
    -V "${final_vcf_dir}/all.filter.snp.vcf" \
    --exclude-filtered \
    -O "${final_vcf_dir}/all.filtered.snp.vcf.gz" || {
    echo "Error: Filtered SNP selection failed"
    exit 1
}

# ========== 5. INDEL processing ==========
echo "Processing INDELs..."

gatk --java-options "${java_opts}" SelectVariants \
    -R "${reference}" \
    -V "${final_vcf_dir}/all.merge_raw.vcf" \
    --select-type INDEL \
    -O "${final_vcf_dir}/all.raw.indel.vcf" || {
    echo "Error: INDEL selection failed"
    exit 1
}

gatk --java-options "${java_opts}" VariantFiltration \
    -R "${reference}" \
    -V "${final_vcf_dir}/all.raw.indel.vcf" \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name 'INDEL_filter' \
    -O "${final_vcf_dir}/all.filter.indel.vcf" || {
    echo "Error: INDEL filtration failed"
    exit 1
}

gatk --java-options "${java_opts}" SelectVariants \
    -R "${reference}" \
    -V "${final_vcf_dir}/all.filter.indel.vcf" \
    --exclude-filtered \
    -O "${final_vcf_dir}/all.filtered.indel.vcf.gz" || {
    echo "Error: Filtered INDEL selection failed"
    exit 1
}

# ========== 6. Index final files ==========
echo "Creating index files..."
tabix -p vcf "${final_vcf_dir}/all.filtered.snp.vcf.gz"
tabix -p vcf "${final_vcf_dir}/all.filtered.indel.vcf.gz"

# ========== 7. Cleanup intermediate files ==========
echo "Cleaning up intermediate files..."
rm -f "${final_vcf_dir}/all.merge_raw.vcf" \
      "${final_vcf_dir}/all.raw.snp.vcf" \
      "${final_vcf_dir}/all.filter.snp.vcf" \
      "${final_vcf_dir}/all.raw.indel.vcf" \
      "${final_vcf_dir}/all.filter.indel.vcf"

echo "Pipeline completed successfully!"