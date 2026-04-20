#!/bin/bash
# GWAS analysis pipeline script - includes VCF conversion, PCA analysis, population structure analysis and association analysis

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

# === Step 1: Convert VCF file to PLINK format ===
# Use vcftools to convert compressed VCF file to PLINK format
vcftools \
  --gzvcf ${BASE_DIR}/gvcf/all.filtered.snp.vcf.gz \    # Input compressed VCF file
  --plink \                                              # Specify output as PLINK format
  --out ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp       # Output file prefix

# === Step 2: Create PLINK binary files ===
# Convert PLINK text files to binary format (.bed/.bim/.fam) for efficiency
${BASE_DIR}/software/plink/plink \
  --file ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp \    # Input files from previous step
  --make-bed \                                           # Generate binary files
  --out ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp       # Output file prefix

# === Step 3: Principal Component Analysis (PCA) ===
# Perform PCA analysis using PLINK (for population stratification correction)
${BASE_DIR}/software/plink/plink \
  --threads 8 \                                          # Use 8 threads to speed up computation
  --bfile ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp \   # Input binary files
  --pca 20 \                                             # Compute first 20 principal components
  --out ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp.pca   # Output prefix

# === Step 4: Population structure analysis (Admixture) ===
# Perform Admixture analysis for population structure (K=2 to 10)
conda activate admixture  # Activate admixture environment

for K in {2..10}; do      # Loop over K=2 to 10 population structures
  admixture \
    --cv \                                                # Perform cross-validation
    ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp.bed \      # Input PLINK binary file
    $K | tee ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp.log${K}.out  # Output log file
done

conda deactivate          # Deactivate environment

# === Step 5: GEMMA association analysis ===
conda activate gemma      # Activate gemma environment

# 5.1 Calculate kinship matrix
gemma -bfile ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp \
      -gk 1 \                                             # Generate standardized kinship matrix
      -o snp_kinship \                                    # Output file prefix
      -outdir ${BASE_DIR}/result/WGS-SVs/gwas/manta       # Output directory

# 5.2 Perform association analysis using linear mixed model
gemma \
  -bfile ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp \
  -k ${BASE_DIR}/result/WGS-SVs/gwas/manta/snp_kinship.cXX.txt \  # Input kinship matrix
  -lmm 1 \                                                # LMM with Wald test
  -o snp_gwas_results \                                   # Output file prefix
  -outdir ${BASE_DIR}/result/WGS-SVs/gwas/manta          # Output directory

conda deactivate          # Deactivate environment

echo "Analysis pipeline completed!"