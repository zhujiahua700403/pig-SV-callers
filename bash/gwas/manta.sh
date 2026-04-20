#!/bin/bash

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_base="${BASE_DIR}/LGS/bam/WGS"
manta_path="${BASE_DIR}/result/LSV/WGS"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_path="${BASE_DIR}/result/LSV/WGS/vcf/manta"
chromosomes="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"
list_path="${BASE_DIR}/result/LSV/list/WGS"
bamfile_list="${BASE_DIR}/result/WGS-SVs/gwas/1.list"

# ========== Initialization ==========
mkdir -p "$manta_path" "$output_path" "$list_path" || {
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

# Check if required tools are available
for tool in samtools bcftools python bgzip tabix; do
    check_dependency "$tool"
done

# Activate conda environment
source activate manta || {
    echo "Error: Failed to activate manta environment"
    exit 1
}

# Check reference genome index
if [[ ! -f "${reference_genome}.fai" ]]; then
    echo "Indexing reference genome..."
    samtools faidx "$reference_genome" || {
        echo "Error: Failed to index reference genome"
        exit 1
    }
fi

# Check BAM file list
if [[ ! -f "$bamfile_list" ]]; then
    echo "Error: BAM file list not found at $bamfile_list"
    exit 1
fi

# ========== Function to find BAM file ==========
find_bam_file() {
    local bamfile="$1"
    local found_bam=""
    
    # Recursively search for matching BAM file under input_base
    found_bam=$(find "$input_base" -name "${bamfile}.*.bam" -type f | head -n1)
    
    if [[ -z "$found_bam" ]]; then
        # If not found with suffix, try directly bamfile.bam
        found_bam=$(find "$input_base" -name "${bamfile}.bam" -type f | head -n1)
    fi
    
    echo "$found_bam"
}

# ========== Main processing function ==========
process_sample() {
    local bamfile="$1"
    local bam_path=$(find_bam_file "$bamfile")
    
    if [[ -z "$bam_path" ]]; then
        echo "Error: Could not find BAM file for $bamfile"
        return 1
    fi
    
    local manta_sample_dir="${manta_path}/${bamfile}"
    local sample_output_dir="${output_path}"
    local final_vcf_gz="${sample_output_dir}/${bamfile}.vcf.gz"
    
    # Check if final output file already exists (skip without deleting)
    if [[ -f "$final_vcf_gz" && -s "$final_vcf_gz" ]]; then
        echo "Final VCF already exists: $final_vcf_gz. Skipping processing for $bamfile."
        return 0
    fi
    
    echo "Processing sample: $bamfile"
    echo "Found BAM file: $bam_path"
    
    # 1. Create sample-specific directory
    mkdir -p "$manta_sample_dir" "$sample_output_dir" || {
        echo "Error: Failed to create directories for $bamfile"
        return 1
    }
    
    # 2. Configure Manta
    echo "Configuring Manta for $bamfile..."
    configManta.py \
        --bam "$bam_path" \
        --referenceFasta "$reference_genome" \
        --outputContig \
        --runDir "$manta_sample_dir" || {
        echo "Error: Manta configuration failed for $bamfile"
        return 1
    }
    
    # 3. Run Manta workflow
    echo "Running Manta workflow for $bamfile..."
    python "${manta_sample_dir}/runWorkflow.py" -j 8 || {
        echo "Error: Manta workflow failed for $bamfile"
        return 1
    }
    
    # 4. Check Manta output
    local manta_output_vcf="${manta_sample_dir}/results/variants/diploidSV.vcf.gz"
    if [[ ! -f "$manta_output_vcf" ]]; then
        echo "Error: Manta output not found at $manta_output_vcf"
        return 1
    fi
    
    # 5. Convert VCF format
    echo "Converting VCF format for $bamfile..."
    local converted_vcf="${sample_output_dir}/${bamfile}.converted.vcf"
    bcftools view "$manta_output_vcf" > "$converted_vcf" || {
        echo "Error: VCF conversion failed for $bamfile"
        return 1
    }
    
    # 6. BND to INV conversion
    echo "Running BND to INV conversion for $bamfile..."
    local final_vcf="${sample_output_dir}/${bamfile}.vcf"
    ${BASE_DIR}/software/convertInversion.py \
        ${BASE_DIR}/anaconda3/envs/manta/libexec/samtools \
        "$reference_genome" \
        "$converted_vcf" > \
        "$final_vcf" || {
        echo "Error: BND to INV conversion failed for $bamfile"
        return 1
    }
    
    # 7. Compress and index
    echo "Compressing and indexing final VCF for $bamfile..."
    bgzip -f "$final_vcf" || {
        echo "Error: Failed to compress final VCF"
        return 1
    }
    tabix -p vcf "${final_vcf}.gz" || {
        echo "Error: Failed to index final VCF"
        return 1
    }
    
    # 8. Clean up temporary files
    rm -f "$converted_vcf"
    echo "Successfully processed sample: $bamfile"
}

# ========== Main processing loop ==========
echo "Starting variant calling pipeline"

# Read bamfile list and process
total_samples=0
processed_samples=0
skipped_samples=0

while IFS= read -r bamfile; do
    # Skip empty lines and comment lines
    [[ -z "$bamfile" ]] && continue
    [[ "$bamfile" =~ ^# ]] && continue
    
    # Remove possible whitespace
    bamfile=$(echo "$bamfile" | xargs)
    
    ((total_samples++))
    
    if process_sample "$bamfile"; then
        ((processed_samples++))
    else
        echo "Warning: Processing failed for sample $bamfile"
    fi
    
done < "$bamfile_list"

echo "Variant calling pipeline completed successfully"
echo "Summary: $processed_samples processed, $((total_samples - processed_samples)) failed out of $total_samples total samples"