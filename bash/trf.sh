#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

output_dir="${BASE_DIR}/result/WGS-SVs/region/trf"
sv_path="${BASE_DIR}/result/WGS-SVs/Truvari/INS_DEL"
vcf_output_dir="${BASE_DIR}/result/WGS-SVs/region/vcf/trf"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham" "base")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
trf_bed="${BASE_DIR}/genome/susScr11.trf.cleaned.bed"

# ========== Initialization ==========
mkdir -p "$output_dir" "$vcf_output_dir"

# Check input files
if [ ! -f "$trf_bed" ]; then
    echo "Error: TRF BED file does not exist: $trf_bed" >&2
    exit 1
fi
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list does not exist: $bamfile_list" >&2
    exit 1
fi

# ========== Processing functions ==========
generate_filtered_bed() {
    local input_vcf=$1
    local output_bed=$2
    
    if [ ! -f "$input_vcf" ]; then
        echo "Warning: VCF file does not exist: $input_vcf" >&2
        return 1
    fi
    
    # Check if VCF file is valid
    if ! bcftools view -h "$input_vcf" &>/dev/null; then
        echo "Warning: VCF file may be corrupted: $input_vcf" >&2
        return 1
    fi
    
    # Generate filtered BED file
    # Use bcftools to extract position information and convert to 0-based coordinates
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "$input_vcf" 2>/dev/null | \
    sed 's/^chr//' | \
    awk 'BEGIN {OFS="\t"} {
        # Handle missing END values
        if ($3 == ".") $3 = $2 + 1;
        
        # Ensure coordinates are valid
        if ($2 < 0) $2 = 0;
        if ($3 < 0) $3 = 0;
        if ($2 >= $3) $3 = $2 + 1;
        
        print $1, $2, $3, $4
    }' > "$output_bed" || {
        echo "Error: Failed to process SV VCF file: $input_vcf" >&2
        return 1
    }
    
    [ -s "$output_bed" ] || {
        echo "Warning: Filtered BED file is empty: $output_bed" >&2
        return 1
    }
    
    return 0
}

process_sample_software() {
    local bamfile=$1
    local software=$2
    
    echo "===== Processing $bamfile - $software ====="
    
    # Define file paths
    sv_vcf="$sv_path/$bamfile/$bamfile.$software.vcf.gz"
    sample_dir="$output_dir/$bamfile/$software"
    vcf_sample_dir="$vcf_output_dir/$bamfile/$software"
    
    # Create output directories
    mkdir -p "$sample_dir" "$vcf_sample_dir" || {
        echo "Error: Failed to create directories" >&2
        return 1
    }

    # Step 1: Generate filtered BED file
    echo "[1/5] Generating filtered BED file..."
    generate_filtered_bed "$sv_vcf" "$sample_dir/sv.bed" || return 1

    # Step 2: Calculate overlap with TRF regions
    echo "[2/5] Calculating TRF overlap..."
    overlap_bed="$sample_dir/sv_trf_overlap80.bed"
    bedtools intersect -a "$sample_dir/sv.bed" -b "$trf_bed" \
        -f 0.8 -wa -wb > "$overlap_bed" || {
        echo "Error: bedtools intersect failed" >&2
        return 1
    }

    # Step 3: Statistics
    echo "[3/5] Collecting statistics..."
    total_svs=$(wc -l < "$sample_dir/sv.bed")
    overlapping_svs=$(wc -l < "$overlap_bed")
    
    if [ "$total_svs" -gt 0 ]; then
        overlap_percent=$(awk "BEGIN {printf \"%.2f\", $overlapping_svs/$total_svs*100}")
        echo "Statistics: total SV=$total_svs, overlapping TRF=$overlapping_svs ($overlap_percent%)"
    else
        echo "Warning: No valid SVs" >&2
        return 1
    fi

    # Step 4: Extract overlapping VCF
    echo "[4/5] Extracting overlapping VCF..."
    mkdir -p "$vcf_sample_dir/overlap" "$vcf_sample_dir/no_overlap"
    
    # Extract overlapping SV IDs
    if [ -s "$overlap_bed" ]; then
        awk '{print $4}' "$overlap_bed" | sort -u > "$sample_dir/trf_overlap_sv_ids.txt"
        
        # Extract VCF records for overlapping SVs
        bcftools view -i "ID=@$sample_dir/trf_overlap_sv_ids.txt" "$sv_vcf" \
            -Oz -o "$vcf_sample_dir/overlap/${bamfile}_${software}_trf_overlap.vcf.gz" && \
        tabix -p vcf "$vcf_sample_dir/overlap/${bamfile}_${software}_trf_overlap.vcf.gz" || {
            echo "Error: Failed to extract overlapping VCF" >&2
            return 1
        }
    fi

    # Step 5: Extract non-overlapping VCF
    echo "[5/5] Extracting non-overlapping VCF..."
    # Get all SV IDs
    awk '{print $4}' "$sample_dir/sv.bed" > "$sample_dir/all_sv_ids.txt"
    
    # Get non-overlapping SV IDs (using comm command)
    comm -23 <(sort "$sample_dir/all_sv_ids.txt") \
            <(sort "$sample_dir/trf_overlap_sv_ids.txt" 2>/dev/null) \
            > "$sample_dir/trf_no_overlap_sv_ids.txt"
    
    no_overlap_count=$(wc -l < "$sample_dir/trf_no_overlap_sv_ids.txt")
    echo "Non-overlapping SV count: $no_overlap_count"
    
    if [ "$no_overlap_count" -gt 0 ]; then
        # Extract VCF records for non-overlapping SVs
        bcftools view -i "ID=@$sample_dir/trf_no_overlap_sv_ids.txt" "$sv_vcf" \
            -Oz -o "$vcf_sample_dir/no_overlap/${bamfile}_${software}_no_trf_overlap.vcf.gz" && \
        tabix -p vcf "$vcf_sample_dir/no_overlap/${bamfile}_${software}_no_trf_overlap.vcf.gz" || {
            echo "Error: Failed to extract non-overlapping VCF" >&2
            return 1
        }
    else
        echo "Warning: No non-overlapping SVs"
    fi

    echo "Completed processing $bamfile - $software"
    return 0
}

# ========== Main pipeline ==========
echo "===== Starting TRF region analysis ====="
echo "Start time: $(date)"

# Read sample list
mapfile -t bamfiles < "$bamfile_list"
echo "Number of samples: ${#bamfiles[@]}"
echo "Number of software tools: ${#softwares[@]}"

# Main processing loop
for bamfile in "${bamfiles[@]}"; do
    for software in "${softwares[@]}"; do
        process_sample_software "$bamfile" "$software"
    done
done

echo "===== Analysis completed ====="
echo "End time: $(date)"
exit 0