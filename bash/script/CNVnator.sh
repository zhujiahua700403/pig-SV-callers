#!/bin/bash
set -e  # Enable error exit mechanism, script exits if any command fails

source activate cnvnator

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

input_paths="${BASE_DIR}/WGS/bam"
output_base="${BASE_DIR}/result/WGS-SVs"
ref_fa="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
bin_size=100
log_dir="${BASE_DIR}/error"
mkdir -p "$log_dir"

bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

cnvnator_dir="${output_base}/CNVnator"
output_dir="${output_base}/Truvari"

process_sample() {
    local bamfile="$1"
    local log_file="$log_dir/${bamfile}.log"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start processing sample: ${bamfile}" | tee -a "$log_file"
    
    mkdir -p "${cnvnator_dir}/${bamfile}" "${output_dir}/${bamfile}"
    root_file="${cnvnator_dir}/${bamfile}/${bamfile}.root"

    bam_path="${input_paths}/${bamfile}.WGS.bam"

    echo "BAM file path: ${bam_path}" | tee -a "$log_file"

    # 1. Generate analysis tree
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating analysis tree..." | tee -a "$log_file"
    cnvnator -root "$root_file" -tree "$bam_path" -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 2>> "$log_file" || {
        echo "ERROR: Failed to generate analysis tree" | tee -a "$log_file"
        return 1
    }

    # 2. Generate histogram
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating histogram..." | tee -a "$log_file"
    cnvnator -root "$root_file" -his $bin_size -d $(dirname "$ref_fa") 2>> "$log_file" || {
        echo "ERROR: Failed to generate histogram" | tee -a "$log_file"
        return 1
    }

    # 3. Statistics calculation
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Calculating statistics..." | tee -a "$log_file"
    cnvnator -root "$root_file" -stat $bin_size 2>> "$log_file" || {
        echo "ERROR: Statistics calculation failed" | tee -a "$log_file"
        return 1
    }

    # 4. Data partitioning
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Partitioning data..." | tee -a "$log_file"
    cnvnator -root "$root_file" -partition $bin_size 2>> "$log_file" || {
        echo "ERROR: Data partitioning failed" | tee -a "$log_file"
        return 1
    }

    # 5. CNV calling
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] CNV calling..." | tee -a "$log_file"
    cnvnator -root "$root_file" -call $bin_size > "${cnvnator_dir}/${bamfile}/${bamfile}.cnv.call.txt" 2>> "$log_file" || {
        echo "ERROR: CNV calling failed" | tee -a "$log_file"
        return 1
    }

    # 6. Format conversion
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converting format..." | tee -a "$log_file"
    cnvnator2VCF.pl -prefix "$bamfile" -reference "$ref_fa" \
        "${cnvnator_dir}/${bamfile}/${bamfile}.cnv.call.txt" $(dirname "$ref_fa") \
        > "${output_dir}/${bamfile}/${bamfile}.CNVnator.vcf" 2>> "$log_file" || {
        echo "ERROR: VCF conversion failed" | tee -a "$log_file"
        return 1
    }

    # 7. Compress and index
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Compressing and indexing VCF file..." | tee -a "$log_file"
    bgzip -f "${output_dir}/${bamfile}/${bamfile}.CNVnator.vcf" 2>> "$log_file" || {
        echo "ERROR: VCF compression failed" | tee -a "$log_file"
        return 1
    }
    tabix -f "${output_dir}/${bamfile}/${bamfile}.CNVnator.vcf.gz" 2>> "$log_file" || {
        echo "ERROR: VCF indexing failed" | tee -a "$log_file"
        return 1
    }

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Successfully completed sample: ${bamfile}" | tee -a "$log_file"
    return 0
}

# ========== Main pipeline ==========

for bamfile in "${bamfiles[@]}"; do
    process_sample "$bamfile" &
done

# Wait for all remaining tasks to finish
wait