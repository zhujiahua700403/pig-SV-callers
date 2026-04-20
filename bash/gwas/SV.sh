#!/bin/bash
# complete_sv_filtering.sh

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

# Input parameters
INPUT_VCF="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper2/merged/all.sorted.genotyped.vcf.gz"
OUTPUT_VCF="${BASE_DIR}/result/WGS-SVs/gwas/graphtyper2/filtered_svs.vcf.gz"

# Temporary files
TEMP0="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper2/merged/temp0_transformed.vcf.gz"
TEMP1="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper2/merged/temp1_filtered.vcf.gz"
TEMP2="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper2/merged/temp2_filtered.vcf.gz"
TEMP3="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper2/merged/temp3_filtered.vcf.gz"
KEPT_VARIANTS="${BASE_DIR}/result/WGS-SVs/gwas/merge/graphtyper2/merged/kept_variants.txt"

# Step 0: Check if input file exists and is BGZF compressed
echo "Step 0: Checking input file..."
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "Error: Input file $INPUT_VCF does not exist."
    exit 1
fi

# Step 0.5: Transform VCF format (modify REF/ALT and chromosome names)
echo "Step 0.5: Transforming VCF format (modifying REF/ALT and chromosome names)..."
# Extract header
bcftools view -h "$INPUT_VCF" > header.txt
# Process body
bcftools view -H "$INPUT_VCF" | awk 'BEGIN {OFS = "\t"} 
    $8 ~ /SVTYPE/ { 
        $4 = "A";
        $5 = "T";
        gsub(/Gm0?/, "", $1);
        print
    } !($8 ~ /SVTYPE/) {
        print
    }' | cat header.txt - | bgzip -c > "$TEMP0"

# Create index
tabix -p vcf -f "$TEMP0"
rm -f header.txt

# Update input file to transformed file
INPUT_VCF="$TEMP0"
echo "VCF transformation completed successfully! Using transformed file as input for further processing."

# Step 1: Keep only INS, DEL, INV, DUP SV types with size <= 1000000bp
echo "Step 1: Filtering for INS/DEL/INV/DUP with size <= 1000000bp..."
# First extract all INS, DEL, INV, DUP variants
bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="DUP"' "$INPUT_VCF" -Oz -o "$TEMP1"
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract INS/DEL/INV/DUP variants."
    exit 1
fi
bcftools index -t -f "$TEMP1"

# Get variant IDs that meet size requirement
bcftools query -f '%ID\t%SVTYPE\t%SVLEN\n' "$TEMP1" | awk '
    ($2 == "INS" || $2 == "DEL" || $2 == "INV" || $2 == "DUP") {
        len = $3;
        if (len == "") len = 0;
        # Take absolute value for size comparison
        abs_len = (len < 0) ? -len : len;
        if (abs_len <= 1000000) {
            print $1;
        }
    }' > "$KEPT_VARIANTS"

if [[ ! -s "$KEPT_VARIANTS" ]]; then
    echo "Warning: No INS/DEL/INV/DUP variants with size <= 1000000bp found. Exiting."
    exit 1
else
    # Keep only variants meeting size requirement
    bcftools view -i "ID=@${KEPT_VARIANTS}" "$TEMP1" -Oz -o "$TEMP2"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter variants by size."
        exit 1
    fi
fi
bcftools index -t -f "$TEMP2"

# Step 2: Apply complex filtering expression using bcftools view
echo "Step 2: Applying complex filtering expression with bcftools view..."
FILTER_EXPR='(SVTYPE = "DEL" && QD > 12 && (ABHet > 0.30 || ABHet < 0) && (AC / NUM_MERGED_SVS) < 25 && PASS_AC > 0 && PASS_ratio > 0.1) || (SVTYPE = "INS" && PASS_AC > 0 && (AC / NUM_MERGED_SVS) < 25 && PASS_ratio > 0.1 && (ABHet > 0.25 || ABHet < 0) && MaxAAS > 4)'
if ! bcftools view -i "$FILTER_EXPR" "$TEMP2" -Oz -o "$TEMP3"; then
    echo "Error: Step 2 (applying complex filtering expression) failed."
    exit 1
fi
bcftools index -t -f "$TEMP3"

# Step 3: Filter sites with FILTER not PASS
echo "Step 3: Filtering sites with FILTER not PASS..."
bcftools view -i 'FILTER="PASS"' "$TEMP3" -Oz -o "$OUTPUT_VCF"
if [ $? -ne 0 ]; then
    echo "Error: Failed to filter sites by FILTER=PASS."
    exit 1
fi

# Index final file
echo "Indexing final VCF..."
bcftools index -t -f "$OUTPUT_VCF"  # Force overwrite existing index

# Step 4: Create reference allele file
echo "Step 4: Creating reference allele file..."
REF_ALLELE_FILE="${OUTPUT_VCF%.vcf.gz}.ref_alleles.txt"
bcftools query -f '%ID\t%REF\n' "$OUTPUT_VCF" > "$REF_ALLELE_FILE" || {
    echo "Error: Failed to create reference allele file"
    exit 1
}

# Step 5: Process VCF to remove '+' characters from ALT fields
echo "Step 5: Removing '+' characters from ALT fields..."
PROCESSED_VCF="${OUTPUT_VCF%.vcf.gz}.processed.vcf.gz"
bcftools view -h "$OUTPUT_VCF" | head -n -1 > header.txt
bcftools view -h "$OUTPUT_VCF" | tail -n 1 >> header.txt
bcftools view -H "$OUTPUT_VCF" | awk 'BEGIN {OFS="\t"} {gsub(/\+/, "", $5); print}' | \
cat header.txt - | bgzip -c > "$PROCESSED_VCF" || {
    echo "Error: Failed to process VCF file"
    rm -f header.txt
    exit 1
}
rm -f header.txt

# Create index for processed VCF
tabix -p vcf -f "$PROCESSED_VCF" || {
    echo "Error: Failed to index processed VCF"
    exit 1
}

# Replace original file
echo "Step 6: Replacing original VCF with processed version..."
mv -f "$PROCESSED_VCF" "$OUTPUT_VCF"
mv -f "${PROCESSED_VCF}.tbi" "${OUTPUT_VCF}.tbi"

# Clean temporary files
echo "Cleaning up temporary files..."
rm -f "$TEMP0" "$TEMP0.tbi" "$TEMP1" "$TEMP1.tbi" "$TEMP2" "$TEMP2.tbi" "$TEMP3" "$TEMP3.tbi" "$KEPT_VARIANTS"

echo "SV filtering and post-processing completed successfully! Output file: $OUTPUT_VCF"
echo "Reference allele file: $REF_ALLELE_FILE"