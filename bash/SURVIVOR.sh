#!/bin/bash

# Activate SURVIVOR environment
source activate SURVIVOR

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

# File path configuration
vcf_path="${BASE_DIR}/result/WGS-SVs/Truvari"
input_paths="${BASE_DIR}/result/WGS-SVs/Truvari/SURVIVOR/vcf"
softwares=("manta" "wham" "dysgu" "SurVIndel2" "smoove" "delly")
SURVIVOR_path="${BASE_DIR}/result/WGS-SVs/SURVIVOR"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari"
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# Read sample list from file
if [ ! -f "$bamfile_list" ]; then
    echo "Error: bamfiles.list not found at $bamfile_list"
    exit 1
fi
readarray -t bamfiles < "$bamfile_list"

# Ensure directories exist
mkdir -p "$input_paths"

# Function to process a single sample
process_bamfile() {
    local bamfile="$1"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing $bamfile..."

    # Create sample output directories
    sample_dir="${SURVIVOR_path}/${bamfile}"
    mkdir -p "${sample_dir}/2to1" "${sample_dir}/3to2" "${sample_dir}/4to2" "${sample_dir}/5to3"
    mkdir -p "${output_path}/SURVIVOR/2to1/${bamfile}" "${output_path}/SURVIVOR/3to2/${bamfile}" "${output_path}/SURVIVOR/4to2/${bamfile}" "${output_path}/SURVIVOR/5to3/${bamfile}"
    mkdir -p "$input_paths/${bamfile}"

    # Copy VCF files to SURVIVOR directory and decompress
    for software in "${softwares[@]}"; do
        vcf_file="${vcf_path}/${bamfile}/${bamfile}.${software}.vcf.gz"
        if [[ -f "$vcf_file" ]]; then
            cp "$vcf_file" "$input_paths/${bamfile}"
            bgzip -d "$input_paths/${bamfile}/${bamfile}.${software}.vcf.gz"
        fi
    done

    # Generate all software combinations
    total_softwares=${#softwares[@]}

    # Function to fix VCF header
    fix_vcf_header() {
        local input_vcf="$1"
        local output_vcf="$2"
        local temp_header=$(mktemp)
        
        # Extract and modify header
        bcftools view -h "$input_vcf" | grep "^##" | \
            sed -e 's/##INFO=<ID=CIEND,Number=2,Type=String/##INFO=<ID=CIEND,Number=2,Type=Integer/' \
                -e 's/##INFO=<ID=CIPOS,Number=2,Type=String/##INFO=<ID=CIPOS,Number=2,Type=Integer/' > "$temp_header"
        
        # Modify VCF using bcftools
        bcftools annotate \
            --remove INFO/CIEND,INFO/CIPOS \
            --output "$output_vcf" \
            --header-lines "$temp_header" \
            "$input_vcf" && \
            bcftools index "$output_vcf"
        
        rm "$temp_header"
    }

    # 2-software combinations
    for ((i=0; i<total_softwares-1; i++)); do
        for ((j=i+1; j<total_softwares; j++)); do
            combo1="${softwares[i]}"
            combo2="${softwares[j]}"
            output_file="${sample_dir}/2to1/${bamfile}.${combo1}_${combo2}.list"
            merged_vcf="${sample_dir}/2to1/${bamfile}.${combo1}_${combo2}_merged.vcf"
            sorted_vcf="${output_path}/SURVIVOR/2to1/${bamfile}/${bamfile}_${combo1}_${combo2}_merged_sorted.vcf"
            fixed_vcf="${output_path}/SURVIVOR/2to1/${bamfile}/${bamfile}_${combo1}_${combo2}_merged_sorted_fixed.vcf.gz"
            filtered_vcf="${output_path}/SURVIVOR/2to1/${bamfile}/${bamfile}_${combo1}_${combo2}_merged_filtered.vcf.gz"

            # Clear old file
            > "$output_file"

            # Build VCF paths for current sample
            file1="${input_paths}/${bamfile}/${bamfile}.${combo1}.vcf"
            file2="${input_paths}/${bamfile}/${bamfile}.${combo2}.vcf"

            # Write file paths (only if files exist)
            [[ -f "$file1" ]] && echo "$file1" >> "$output_file"
            [[ -f "$file2" ]] && echo "$file2" >> "$output_file"

            # Check if file list is valid
            if [[ $(wc -l < "$output_file") -ge 2 ]]; then
                # Merge VCFs using SURVIVOR (parameters: 1000 1 0 0 0 50)
                SURVIVOR merge "$output_file" 1000 1 0 0 0 50 "$merged_vcf"

                # Sort and compress VCF
                bcftools sort "${merged_vcf}" -O z -o "$sorted_vcf.gz"
                rm -f "$merged_vcf"

                # Fix VCF header
                fix_vcf_header "$sorted_vcf.gz" "$fixed_vcf"
                
                # Keep only INS and DEL types
                bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL"' "$fixed_vcf" -O z -o "$filtered_vcf"
                tabix "$filtered_vcf"
                
                # Remove intermediate files
                rm -f "$sorted_vcf.gz" "$fixed_vcf" "${fixed_vcf}.tbi"
            else
                rm -f "$output_file" "$sorted_vcf.gz"
            fi
        done
    done

    # 3-software combinations
    for ((i=0; i<total_softwares-2; i++)); do
        for ((j=i+1; j<total_softwares-1; j++)); do
            for ((k=j+1; k<total_softwares; k++)); do
                combo1="${softwares[i]}"
                combo2="${softwares[j]}"
                combo3="${softwares[k]}"
                output_file="${sample_dir}/3to2/${bamfile}.${combo1}_${combo2}_${combo3}.list"
                merged_vcf="${sample_dir}/3to2/${bamfile}.${combo1}_${combo2}_${combo3}_merged.vcf"
                sorted_vcf="${output_path}/SURVIVOR/3to2/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_merged_sorted.vcf"
                fixed_vcf="${output_path}/SURVIVOR/3to2/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_merged_sorted_fixed.vcf.gz"
                filtered_vcf="${output_path}/SURVIVOR/3to2/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_merged_filtered.vcf.gz"

                > "$output_file"

                file1="${input_paths}/${bamfile}/${bamfile}.${combo1}.vcf"
                file2="${input_paths}/${bamfile}/${bamfile}.${combo2}.vcf"
                file3="${input_paths}/${bamfile}/${bamfile}.${combo3}.vcf"

                [[ -f "$file1" ]] && echo "$file1" >> "$output_file"
                [[ -f "$file2" ]] && echo "$file2" >> "$output_file"
                [[ -f "$file3" ]] && echo "$file3" >> "$output_file"

                if [[ $(wc -l < "$output_file") -ge 2 ]]; then
                    SURVIVOR merge "$output_file" 1000 2 1 1 0 50 "$merged_vcf"
                    bcftools sort "${merged_vcf}" -O z -o "$sorted_vcf.gz"
                    rm -f "$merged_vcf"
                    fix_vcf_header "$sorted_vcf.gz" "$fixed_vcf"
                    
                    bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL"' "$fixed_vcf" -O z -o "$filtered_vcf"
                    tabix "$filtered_vcf"
                    
                    rm -f "$sorted_vcf.gz" "$fixed_vcf" "${fixed_vcf}.tbi"
                else
                    rm -f "$output_file" "$sorted_vcf.gz"
                fi
            done
        done
    done

    # 4-software combinations
    for ((i=0; i<total_softwares-3; i++)); do
        for ((j=i+1; j<total_softwares-2; j++)); do
            for ((k=j+1; k<total_softwares-1; k++)); do
                for ((l=k+1; l<total_softwares; l++)); do
                    combo1="${softwares[i]}"
                    combo2="${softwares[j]}"
                    combo3="${softwares[k]}"
                    combo4="${softwares[l]}"
                    output_file="${sample_dir}/4to2/${bamfile}.${combo1}_${combo2}_${combo3}_${combo4}.list"
                    merged_vcf="${sample_dir}/4to2/${bamfile}.${combo1}_${combo2}_${combo3}_${combo4}_merged.vcf"
                    sorted_vcf="${output_path}/SURVIVOR/4to2/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_${combo4}_merged_sorted.vcf"
                    fixed_vcf="${output_path}/SURVIVOR/4to2/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_${combo4}_merged_sorted_fixed.vcf.gz"
                    filtered_vcf="${output_path}/SURVIVOR/4to2/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_${combo4}_merged_filtered.vcf.gz"

                    > "$output_file"

                    file1="${input_paths}/${bamfile}/${bamfile}.${combo1}.vcf"
                    file2="${input_paths}/${bamfile}/${bamfile}.${combo2}.vcf"
                    file3="${input_paths}/${bamfile}/${bamfile}.${combo3}.vcf"
                    file4="${input_paths}/${bamfile}/${bamfile}.${combo4}.vcf"

                    [[ -f "$file1" ]] && echo "$file1" >> "$output_file"
                    [[ -f "$file2" ]] && echo "$file2" >> "$output_file"
                    [[ -f "$file3" ]] && echo "$file3" >> "$output_file"
                    [[ -f "$file4" ]] && echo "$file4" >> "$output_file"

                    if [[ $(wc -l < "$output_file") -ge 2 ]]; then
                        SURVIVOR merge "$output_file" 1000 2 1 1 0 50 "$merged_vcf"
                        bcftools sort "${merged_vcf}" -O z -o "$sorted_vcf.gz"
                        rm -f "$merged_vcf"
                        fix_vcf_header "$sorted_vcf.gz" "$fixed_vcf"
                        
                        bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL"' "$fixed_vcf" -O z -o "$filtered_vcf"
                        tabix "$filtered_vcf"
                        
                        rm -f "$sorted_vcf.gz" "$fixed_vcf" "${fixed_vcf}.tbi"
                    else
                        rm -f "$output_file" "$sorted_vcf.gz"
                    fi
                done
            done
        done
    done

    # 5-software combinations
    for ((i=0; i<total_softwares-4; i++)); do
        for ((j=i+1; j<total_softwares-3; j++)); do
            for ((k=j+1; k<total_softwares-2; k++)); do
                for ((l=k+1; l<total_softwares-1; l++)); do
                    for ((m=l+1; m<total_softwares; m++)); do
                        combo1="${softwares[i]}"
                        combo2="${softwares[j]}"
                        combo3="${softwares[k]}"
                        combo4="${softwares[l]}"
                        combo5="${softwares[m]}"
                        output_file="${sample_dir}/5to3/${bamfile}.${combo1}_${combo2}_${combo3}_${combo4}_${combo5}.list"
                        merged_vcf="${sample_dir}/5to3/${bamfile}.${combo1}_${combo2}_${combo3}_${combo4}_${combo5}_merged.vcf"
                        sorted_vcf="${output_path}/SURVIVOR/5to3/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_${combo4}_${combo5}_merged_sorted.vcf"
                        fixed_vcf="${output_path}/SURVIVOR/5to3/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_${combo4}_${combo5}_merged_sorted_fixed.vcf.gz"
                        filtered_vcf="${output_path}/SURVIVOR/5to3/${bamfile}/${bamfile}_${combo1}_${combo2}_${combo3}_${combo4}_${combo5}_merged_filtered.vcf.gz"

                        > "$output_file"

                        file1="${input_paths}/${bamfile}/${bamfile}.${combo1}.vcf"
                        file2="${input_paths}/${bamfile}/${bamfile}.${combo2}.vcf"
                        file3="${input_paths}/${bamfile}/${bamfile}.${combo3}.vcf"
                        file4="${input_paths}/${bamfile}/${bamfile}.${combo4}.vcf"
                        file5="${input_paths}/${bamfile}/${bamfile}.${combo5}.vcf"

                        [[ -f "$file1" ]] && echo "$file1" >> "$output_file"
                        [[ -f "$file2" ]] && echo "$file2" >> "$output_file"
                        [[ -f "$file3" ]] && echo "$file3" >> "$output_file"
                        [[ -f "$file4" ]] && echo "$file4" >> "$output_file"
                        [[ -f "$file5" ]] && echo "$file5" >> "$output_file"

                        if [[ $(wc -l < "$output_file") -ge 2 ]]; then
                            SURVIVOR merge "$output_file" 1000 3 1 1 0 50 "$merged_vcf"
                            bcftools sort "${merged_vcf}" -O z -o "$sorted_vcf.gz"
                            rm -f "$merged_vcf"
                            fix_vcf_header "$sorted_vcf.gz" "$fixed_vcf"
                            
                            bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL"' "$fixed_vcf" -O z -o "$filtered_vcf"
                            tabix "$filtered_vcf"
                            
                            rm -f "$sorted_vcf.gz" "$fixed_vcf" "${fixed_vcf}.tbi"
                        else
                            rm -f "$output_file" "$sorted_vcf.gz"
                        fi
                    done
                done
            done
        done
    done
}

# Process all samples
for bamfile in "${bamfiles[@]}"; do
    process_bamfile "$bamfile"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] All BAM processing completed!"