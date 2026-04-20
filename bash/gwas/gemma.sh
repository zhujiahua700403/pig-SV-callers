#!/bin/bash
# SV-GWAS analysis pipeline integration script

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

WORK_DIR="${BASE_DIR}/result/WGS-SVs/gwas/result/1"
VCF_FILE="${BASE_DIR}/result/WGS-SVs/gwas/graphtyper2/filtered_svs.vcf.gz"
PLINK="${BASE_DIR}/software/plink/plink"

# Change to working directory
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

echo "Starting SV-GWAS analysis pipeline..."

# Restore original phenotype file
if [ -f "${WORK_DIR}/pheno.txt.backup" ]; then
    echo "Restoring original phenotype file..."
    cp ${WORK_DIR}/pheno.txt.backup ${WORK_DIR}/pheno.txt
fi

# Step 1: Convert VCF to PLINK format
echo "Step 1: Converting VCF to PLINK format..."
${PLINK} --vcf ${VCF_FILE} --recode --out ${WORK_DIR}/pig --make-bed --allow-extra-chr

# Calculate missingness
${PLINK} --bfile ${WORK_DIR}/pig --missing --out ${WORK_DIR}/pig.plink_missing

# Step 2: Quality control filtering
echo "Step 2: Performing MAF and missingness filtering..."
${PLINK} --bfile ${WORK_DIR}/pig --maf 0.05 --geno 0.1 --make-bed --mind 0.9 --out pig_filtered

# Calculate missingness after filtering
${PLINK} --bfile ${WORK_DIR}/pig_filtered --missing --out ${WORK_DIR}/pig_filtered.plink_missing

# Additional step: Compare missingness files before and after filtering, update phenotype file
echo "Additional step: Comparing missingness files before/after filtering, updating phenotype file..."

# Check if files exist
if [ ! -f "${WORK_DIR}/pig_filtered.plink_missing.imiss" ]; then
    echo "Error: Filtered missingness file does not exist"
    exit 1
fi

if [ ! -f "${WORK_DIR}/pheno.txt" ]; then
    echo "Error: Phenotype file does not exist"
    exit 1
fi

# Extract sample IDs retained after filtering (second column, sample ID)
awk 'NR>1 {print $2}' ${WORK_DIR}/pig_filtered.plink_missing.imiss > ${WORK_DIR}/kept_samples.txt

# Check number of retained samples
n_kept_samples=$(wc -l < ${WORK_DIR}/kept_samples.txt)
echo "Number of retained samples after filtering: $n_kept_samples"

# Backup original phenotype file
cp ${WORK_DIR}/pheno.txt ${WORK_DIR}/pheno.txt.backup

# Update phenotype file, keeping only retained samples
awk '
BEGIN {
    # Read retained sample IDs
    while (getline < "'${WORK_DIR}'/kept_samples.txt") {
        kept[$1] = 1
    }
}
{
    # Check if the first column (sample ID) is in the kept set
    if ($1 in kept) {
        print $0
    }
}' ${WORK_DIR}/pheno.txt > ${WORK_DIR}/pheno_filtered.txt

# Check number of lines in updated phenotype file
n_pheno=$(wc -l < ${WORK_DIR}/pheno_filtered.txt)
echo "Number of lines in updated phenotype file: $n_pheno"

# Check number of samples after filtering
n_samples=$(wc -l < ${WORK_DIR}/pig_filtered.fam)
echo "Number of samples after filtering: $n_samples"

if [ $n_pheno -eq 0 ]; then
    echo "Error: Updated phenotype file is empty"
    echo "Check sample ID format consistency:"
    echo "First 5 lines of phenotype file:"
    head -n 5 ${WORK_DIR}/pheno.txt
    echo "First 5 lines of retained samples file:"
    head -n 5 ${WORK_DIR}/kept_samples.txt
    exit 1
fi

if [ $n_samples -ne $n_pheno ]; then
    echo "Warning: Number of lines in phenotype file ($n_pheno) does not match filtered sample count ($n_samples)"
    echo "Attempting to update phenotype file using samples from FAM file..."
    
    # Use samples from FAM file to update phenotype file
    awk '{print $2}' ${WORK_DIR}/pig_filtered.fam > ${WORK_DIR}/fam_samples.txt
    
    awk '
    BEGIN {
        # Read sample IDs from FAM file
        while (getline < "'${WORK_DIR}'/fam_samples.txt") {
            fam[$1] = 1
        }
    }
    NR==FNR {
        # First pass: read phenotype file, build index
        pheno[$1] = $0
        next
    }
    ($1 in fam) {
        # If sample exists in FAM file, output phenotype
        if ($1 in pheno) {
            print pheno[$1]
        } else {
            # If no phenotype data, output default values
            # Output appropriate number of columns based on phenotype file structure
            printf "0\t%s", $1
            # Output enough -9 (assuming 16 phenotypes)
            for(i=1; i<=16; i++) {
                printf "\t-9"
            }
            printf "\n"
        }
    }' ${WORK_DIR}/pheno.txt ${WORK_DIR}/fam_samples.txt > ${WORK_DIR}/pheno_filtered2.txt
    
    # Use the newly created phenotype file
    mv ${WORK_DIR}/pheno_filtered2.txt ${WORK_DIR}/pheno.txt
    n_pheno=$(wc -l < ${WORK_DIR}/pheno.txt)
    echo "New phenotype file line count: $n_pheno"
else
    # Use the updated phenotype file
    mv ${WORK_DIR}/pheno_filtered.txt ${WORK_DIR}/pheno.txt
fi

echo "Phenotype file update completed, retained $n_pheno samples"

# Modification step: Change the 6th column of the FAM file before PCA calculation
echo "Modification step: Changing the 6th column (phenotype value) of FAM file from -9 to 1 before PCA calculation..."
# Backup original FAM file
cp ${WORK_DIR}/pig_filtered.fam ${WORK_DIR}/pig_filtered.fam.backup

# Change the 6th column from -9 to 1
awk '{
    if ($6 == -9) {
        $6 = 1
    }
    print $0
}' ${WORK_DIR}/pig_filtered.fam > ${WORK_DIR}/pig_filtered.fam.tmp

# Replace original file
mv ${WORK_DIR}/pig_filtered.fam.tmp ${WORK_DIR}/pig_filtered.fam

echo "FAM file modified, changed -9 to 1 in the 6th column"

# Step 3: Calculate PCA
echo "Step 3: Calculating PCA..."
${PLINK} --allow-extra-chr --threads 2 --bfile ${WORK_DIR}/pig_filtered --pca 3 --out ${WORK_DIR}/pig_filtered

# Additional step: Compare sample order between pig_filtered.eigenvec and pig_filtered.fam
echo "Additional step: Comparing sample order between pig_filtered.eigenvec and pig_filtered.fam..."

# Extract sample IDs from eigenvec file (combine first two columns as FID_IID)
awk '{print $1 "_" $2}' ${WORK_DIR}/pig_filtered.eigenvec > ${WORK_DIR}/eigenvec_samples.tmp

# Extract sample IDs from fam file (combine first two columns as FID_IID)
awk '{print $1 "_" $2}' ${WORK_DIR}/pig_filtered.fam > ${WORK_DIR}/fam_samples.tmp

# Compare the two files
if diff ${WORK_DIR}/eigenvec_samples.tmp ${WORK_DIR}/fam_samples.tmp > /dev/null; then
    echo "Sample order consistent: eigenvec and fam files have the same sample order"
else
    echo "Warning: Sample order between eigenvec and fam files is inconsistent"
    echo "Differences:"
    diff ${WORK_DIR}/eigenvec_samples.tmp ${WORK_DIR}/fam_samples.tmp
fi

# Clean temporary files
rm ${WORK_DIR}/eigenvec_samples.tmp ${WORK_DIR}/fam_samples.tmp

# Step 4: Convert PCA file format and add sex covariate
echo "Step 4: Converting PCA file format and adding sex covariate..."

# Check if sex file exists
if [ ! -f "${WORK_DIR}/sex.txt" ]; then
    echo "Error: Sex file sex.txt does not exist"
    exit 1
fi

# Create a temporary file merging sex information with PCA
# First create a temporary file with sex information
awk '
BEGIN {
    # Read sex file, map sample ID to sex
    while (getline < "'${WORK_DIR}'/sex.txt") {
        sex[$1] = $2
    }
}
{
    # Extract sample ID (second column of eigenvec)
    sample_id = $2
    # Look up corresponding sex
    if (sample_id in sex) {
        gender = sex[sample_id]
    } else {
        # If no sex information, use NA
        gender = "NA"
    }
    # Output sample ID and sex
    print sample_id, gender
}' ${WORK_DIR}/pig_filtered.eigenvec > ${WORK_DIR}/gender_info.tmp

# Merge sex information with PCA
paste -d " " ${WORK_DIR}/gender_info.tmp ${WORK_DIR}/pig_filtered.eigenvec | awk '
{
    # Output format: 1 sex PCA1 PCA2 PCA3
    printf "1 %s", $2  # First column is intercept 1, second column is sex
    # Add PCA values (starting from column 5)
    for(i=5; i<=NF; i++) {
        printf " %s", $i
    }
    printf "\n"
}' > pca.txt

# Clean temporary file
rm ${WORK_DIR}/gender_info.tmp

echo "PCA file (with sex covariate) conversion completed"

# Step 5: Activate conda environment and calculate kinship matrix
echo "Step 5: Calculating kinship matrix..."
source activate gemma
gemma -bfile ${WORK_DIR}/pig_filtered -gk 2 -o pig.kinship -p pheno.txt

# Step 6: Association analysis (loop over phenotypes 1-16)
echo "Step 6: Starting association analysis (phenotypes 1-16)..."

# Get actual name of kinship matrix file
KINSHIP_FILE=$(ls ${WORK_DIR}/output/pig.kinship.s*.txt | head -n 1)

for PHENO_NUM in {1..16}; do
    echo "Analyzing phenotype ${PHENO_NUM}..."
    gemma \
        -bfile pig_filtered \
        -k ${KINSHIP_FILE} \
        -lmm 4 \
        -n ${PHENO_NUM} \
        -c pca.txt \
        -p pheno.txt \
        -o gwas_results_pheno${PHENO_NUM}
done

echo "All analyses completed!"
echo "Result files are saved in: ${WORK_DIR}/output/"