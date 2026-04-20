#!/bin/bash

# Set logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

# Decompress input file
log "Starting to decompress RepeatMasker result file..."
input_file="${BASE_DIR}/genome/susScr11.fa.out"
if [ ! -f "${input_file}.gz" ]; then
    log "Error: Input file ${input_file}.gz does not exist"
    exit 1
fi
gunzip -c "${input_file}.gz" > "$input_file" || {
    log "Error: Failed to decompress susScr11.fa.out.gz"
    exit 1
}

# Define output files
output_dir="${BASE_DIR}/result/WGS-SVs/region/RepeatMasker"
sv_path="${BASE_DIR}/result/WGS-SVs/Truvari/INS_DEL"
vcf_output_dir="${BASE_DIR}/result/WGS-SVs/region/vcf/RepeatMasker"
softwares=("breakdancer" "CNVnator" "delly" "dysgu" "gridss" "insurveyor" "manta" "smoove" "SurVIndel2" "wham" "base")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"
regions=("DNA" "LINE" "SINE" "LTR" "Low_complexity" "Simple_repeat" "Other")

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    log "Error: bamfiles.list does not exist at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Create output directories
mkdir -p "$output_dir" "$vcf_output_dir" || {
    log "Error: Failed to create output directories"
    exit 1
}

# Define repeat class BED files
declare -A repeat_beds=(
    ["DNA"]="$output_dir/DNA.bed"
    ["LINE"]="$output_dir/LINE.bed"
    ["SINE"]="$output_dir/SINE.bed"
    ["LTR"]="$output_dir/LTR.bed"
    ["Low_complexity"]="$output_dir/Low_complexity.bed"
    ["Simple_repeat"]="$output_dir/Simple_repeat.bed"
    ["Other"]="$output_dir/Other.bed"
)

# Clean old BED files
for region in "${regions[@]}"; do
    rm -f "${repeat_beds[$region]}"
done

# Process RepeatMasker output and generate classified BED files
log "Processing RepeatMasker output file..."
if [ ! -s "$input_file" ]; then
    log "Error: RepeatMasker output file is empty or does not exist"
    exit 1
fi

awk 'BEGIN{OFS="\t"} 
!/^#/ && !/^$/ && NF >= 10 {
    # Extract chromosome and normalize
    chrom = $5;
    sub(/^chr/, "", chrom);
    
    # Extract coordinates and convert to 0-based
    start = $6 - 1;  # Convert to 0-based
    end = $7;
    
    # Classify repeat type
    class_family = $11;
    if (class_family ~ /DNA/) {type="DNA"}
    else if (class_family ~ /LINE/) {type="LINE"}
    else if (class_family ~ /SINE/) {type="SINE"}
    else if (class_family ~ /LTR/) {type="LTR"}
    else if (class_family ~ /Low_complexity/) {type="Low_complexity"}
    else if (class_family ~ /Simple_repeat/) {type="Simple_repeat"}
    else {type="Other"}
    
    # Output only valid coordinates
    if (start >= 0 && end > start) {
        print chrom, start, end > "'"$output_dir/"'" type ".bed"
    } else {
        print "Warning: Invalid coordinates - chrom:" chrom ", position:" start+1 "-" end | "cat >&2"
    }
}' "$input_file" || {
    log "Error: Failed to process RepeatMasker output file"
    exit 1
}

# Check generated BED files
for region in "${regions[@]}"; do
    bed_file="${repeat_beds[$region]}"
    if [ ! -s "$bed_file" ]; then
        log "Warning: ${region} BED file is empty or not generated"
        # Create empty file to avoid later errors
        touch "$bed_file"
    else
        # Validate BED file format
        if ! head -n 1 "$bed_file" | awk '{exit ($2 < 0 || $3 < 0 || $2 >= $3 ? 1 : 0)}'; then
            log "Fixing ${region} BED file format issues..."
            awk 'BEGIN{OFS="\t"} {
                if ($2 < 0) $2 = 0;
                if ($3 < 0) $3 = 0;
                if ($2 > $3) {tmp = $2; $2 = $3; $3 = tmp};
                if ($2 < $3) print $0
            }' "$bed_file" > "${bed_file}.tmp" && mv "${bed_file}.tmp" "$bed_file"
        fi
    fi
done

# Function to generate and filter SV BED file
generate_filtered_bed() {
    local input_vcf=$1
    local output_bed=$2
    local temp_bed="${output_bed}.tmp"
    
    if [ ! -f "$input_vcf" ]; then
        log "Warning: VCF file does not exist $input_vcf"
        return 1
    fi
    
    # Check if VCF file is valid
    if ! bcftools view -h "$input_vcf" &>/dev/null; then
        log "Warning: VCF file may be corrupted or incorrectly formatted $input_vcf"
        return 1
    fi
    
    # Use bcftools query to extract information and filter negative values
    log "Generating and filtering SV BED file..."
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "$input_vcf" 2>/dev/null | \
    sed 's/^chr//' | \
    awk 'BEGIN {OFS="\t"} {
        # Handle missing END values
        if ($3 == ".") $3 = $2 + 1;
        
        # Check if coordinates are negative
        if ($2 < 0) $2 = 0;
        if ($3 < 0) $3 = 0;
        
        # Ensure start <= end
        if ($2 > $3) {
            tmp = $2;
            $2 = $3;
            $3 = tmp;
        }
        
        # Ensure minimum length
        if ($3 - $2 < 1) $3 = $2 + 1;
        
        print $1, $2, $3, $4
    }' > "$temp_bed" || {
        log "Error: Failed to process SV VCF file $input_vcf"
        return 1
    }
    
    # Check if generated BED file is empty
    if [ ! -s "$temp_bed" ]; then
        log "Warning: Generated BED file is empty $temp_bed"
        rm -f "$temp_bed"
        return 1
    fi
    
    # Move temporary file to final location
    mv "$temp_bed" "$output_bed" || {
        log "Error: Failed to rename temporary BED file"
        return 1
    }
    
    return 0
}

# Main pipeline
log "Starting processing of samples and software..."
for bamfile in "${bamfiles[@]}"; do
    for software in "${softwares[@]}"; do
        log "Processing sample $bamfile with software $software"
        
        # Define file paths
        sv_vcf="$sv_path/$bamfile/$bamfile.$software.vcf.gz"
        sample_dir="$output_dir/$bamfile/$software"
        vcf_sample_dir="$vcf_output_dir/$bamfile/$software"
        
        # Create output directories
        mkdir -p "$sample_dir" "$vcf_sample_dir" "$vcf_sample_dir/no" || {
            log "Error: Failed to create directories $sample_dir or $vcf_sample_dir"
            continue
        }

        # Check if VCF file exists
        if [ ! -f "$sv_vcf" ]; then
            log "Warning: VCF file does not exist $sv_vcf"
            continue
        fi

        # Generate and filter SV BED file
        if ! generate_filtered_bed "$sv_vcf" "$sample_dir/sv.bed"; then
            continue
        fi
        
        # Create empty file to record all overlapping SV IDs
        all_overlap_ids_file="$sample_dir/all_overlap_ids.txt"
        > "$all_overlap_ids_file"  # Clear file

        # Compare with genomic repeat regions
        log "Comparing with genomic repeat regions..."
        for region in "${regions[@]}"; do
            region_bed="${repeat_beds[$region]}"
            overlap_bed="$sample_dir/sv_${region}_overlap80.bed"
            count_file="$sample_dir/sv_${region}_count.txt"
            region_vcf_dir="$vcf_output_dir/$bamfile/$software/$region"
            
            if [ ! -s "$region_bed" ]; then
                log "Warning: ${region} BED file is empty $region_bed"
                continue
            fi
            
            # Calculate 80% overlap regions
            bedtools intersect -a "$sample_dir/sv.bed" \
                               -b "$region_bed" -f 0.8 -wa -wb > "$overlap_bed" 2>/dev/null || {
                log "Error: bedtools intersect failed, possibly due to empty input files"
                continue
            }
            
            # Count overlaps
            total_svs=$(wc -l < "$sample_dir/sv.bed" 2>/dev/null || echo 0)
            overlapping_svs=$(wc -l < "$overlap_bed" 2>/dev/null || echo 0)
            
            if [ "$total_svs" -gt 0 ]; then
                overlap_percent=$(awk -v total="$total_svs" -v overlap="$overlapping_svs" \
                    'BEGIN {printf "%.2f", overlap/total*100}')
                echo -e "Total SV count\t${total_svs}\nOverlap with ${region} region SV count\t${overlapping_svs}\nOverlap percentage\t${overlap_percent}%" > "$count_file"
                log "Sample $bamfile ($software) ${region}: total SV=$total_svs, overlapping SV=$overlapping_svs ($overlap_percent%)"
            else
                log "Warning: Sample $bamfile ($software) has no valid SVs detected"
                continue
            fi
            
            # If there are overlapping regions, extract corresponding VCF records
            if [ -s "$overlap_bed" ]; then
                mkdir -p "$region_vcf_dir" || {
                    log "Error: Failed to create directory $region_vcf_dir"
                    continue
                }
                
                # Extract overlapping SV IDs
                awk '{print $4}' "$overlap_bed" | sort -u > "$sample_dir/${region}_overlap_sv_ids.txt"
                
                # Add to all overlap IDs file
                cat "$sample_dir/${region}_overlap_sv_ids.txt" >> "$all_overlap_ids_file"
                
                # Extract corresponding VCF records
                bcftools view -i "ID=@$sample_dir/${region}_overlap_sv_ids.txt" "$sv_vcf" \
                    -Oz -o "$region_vcf_dir/${bamfile}_${software}_${region}_overlap.vcf.gz" && \
                tabix -p vcf "$region_vcf_dir/${bamfile}_${software}_${region}_overlap.vcf.gz" 2>/dev/null || {
                    log "Error: Failed to extract and index overlapping VCF file"
                    continue
                }
                
                log "Successfully generated ${region} overlapping VCF file: $region_vcf_dir/${bamfile}_${software}_${region}_overlap.vcf.gz"
            else
                log "Warning: No SVs overlapping with ${region} region found"
            fi
        done
        
        # ==================================================================
        # New section: Extract SVs not in any repeat sequence
        # ==================================================================
        log "Extracting SVs not in any repeat sequence..."
        
        # 1. Get all SV IDs
        awk '{print $4}' "$sample_dir/sv.bed" > "$sample_dir/all_sv_ids.txt"
        
        # 2. Deduplicate all overlapping IDs
        if [ -s "$all_overlap_ids_file" ]; then
            sort -u "$all_overlap_ids_file" > "$sample_dir/all_overlap_sv_ids.unique.txt"
        else
            > "$sample_dir/all_overlap_sv_ids.unique.txt"  # Create empty file
        fi
        
        # 3. Find SV IDs not in repeat sequences
        comm -23 <(sort "$sample_dir/all_sv_ids.txt") <(sort "$sample_dir/all_overlap_sv_ids.unique.txt") \
            > "$sample_dir/no_overlap_sv_ids.txt"
        
        # 4. Count number
        no_overlap_count=$(wc -l < "$sample_dir/no_overlap_sv_ids.txt")
        log "Number of SVs not in any repeat sequence: $no_overlap_count"
        
        # 5. Extract corresponding VCF records
        if [ "$no_overlap_count" -gt 0 ]; then
            bcftools view -i "ID=@$sample_dir/no_overlap_sv_ids.txt" "$sv_vcf" \
                -Oz -o "$vcf_sample_dir/no/${bamfile}_${software}_no_overlap.vcf.gz" && \
            tabix -p vcf "$vcf_sample_dir/no/${bamfile}_${software}_no_overlap.vcf.gz" 2>/dev/null || {
                log "Error: Failed to extract and index non-overlapping VCF file"
                continue
            }
            log "Successfully generated non-overlapping VCF file: $vcf_sample_dir/no/${bamfile}_${software}_no_overlap.vcf.gz"
        else
            log "Warning: No SVs outside repeat sequences found"
        fi
    done
done

log "All samples and software analysis completed! Results saved in: $output_dir and $vcf_output_dir"