#!/bin/bash
set -euo pipefail  # Enable strict error handling

# ========== Configuration Parameters ==========
BASE_DIR="/path/to/your/folder"

source activate truvari  # Activate conda environment

input_path="${BASE_DIR}/result/WGS-SVs/Truvari"
output_path="${BASE_DIR}/result/WGS-SVs/Truvari/SURVIVOR/result2"
reference_genome="${BASE_DIR}/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa"
combines=("2to1" "3to2" "4to2" "5to3")
svtypes=("BND" "INV" "DUP")
bamfile_list="${BASE_DIR}/result/WGS-SVs/bamfiles.list"

# ========== Initialization ==========
mkdir -p "$output_path"

# ========== Sample processing function ==========
process_combination() {
    local bamfile="$1"
    local combine="$2"
    
    echo "===== Processing $bamfile - $combine ====="
    
    # Get VCF file list
    vcf_dir="$input_path/SURVIVOR/vcf2/$combine/$bamfile"
    if [ ! -d "$vcf_dir" ]; then
        echo "Warning: Directory does not exist: $vcf_dir"
        return 1
    fi

    # Extract software combination list
    software_list=()
    while IFS= read -r -d '' file; do
        filename=$(basename "$file")
        # Extract software part (remove sample prefix and SV type suffix)
        software_part=${filename#"${bamfile}_"}
        software_part=${software_part%"_merged_"*".vcf.gz"}
        software_list+=("$software_part")
    done < <(find "$vcf_dir" -name "${bamfile}_*_merged_*.vcf.gz" -print0)

    # Remove duplicates
    readarray -t unique_software_list < <(printf "%s\n" "${software_list[@]}" | sort -u)
    
    echo "Found ${#unique_software_list[@]} software combinations"

    # Process each software combination
    for software in "${unique_software_list[@]}"; do
        echo "Processing combination: $software"
        
        # Process each SV type
        for svtype in "${svtypes[@]}"; do
            echo "Processing SV type: $svtype"
            
            # Build input/output paths
            input_vcf="$vcf_dir/${bamfile}_${software}_merged_${svtype}.vcf.gz"
            output_subdir="$output_path/$bamfile/$combine/$svtype"
            base_vcf="${BASE_DIR}/result/WGS-SVs/Truvari/svtype/$bamfile/base/$bamfile.base.${svtype}.vcf.gz"
            
            # Check input files
            if [ ! -f "$input_vcf" ]; then
                echo "Warning: Input VCF does not exist: $input_vcf"
                continue
            fi
            
            if [ ! -f "$base_vcf" ]; then
                echo "Error: Base VCF does not exist: $base_vcf"
                continue
            fi

            # Create output directory
            mkdir -p "$output_subdir"
            rm -rf "$output_subdir/$software"

            # Run truvari
            echo "Running truvari:"
            echo "  Base VCF: $base_vcf"
            echo "  Input VCF: $input_vcf"
            echo "  Output directory: $output_subdir/$software"
            
            if truvari bench \
                -b "$base_vcf" \
                -c "$input_vcf" \
                -o "$output_subdir/$software" \
                -f "$reference_genome" \
                --sizemax 100000 \
                -r 1000 -p 0 -P 0.5 -O 0.5 \
                --passonly > "$output_subdir/${software}_truvari.log" 2>&1
            then
                echo "Successfully completed: $software - $svtype"
            else
                echo "Error: truvari failed: $software - $svtype"
            fi
        done
    done
    
    echo "Finished processing: $bamfile - $combine"
}

# ========== Main Process ==========
echo "===== Starting processing ====="
echo "Start time: $(date)"

# Read sample list
mapfile -t bamfiles < "$bamfile_list"
echo "Number of samples: ${#bamfiles[@]}"
echo "Number of combinations: ${#combines[@]}"
echo "Number of SV types: ${#svtypes[@]}"

# Main processing loop
for bamfile in "${bamfiles[@]}"; do
    echo "===== Starting sample $bamfile ====="
    
    # Start background tasks for all combinations
    for combine in "${combines[@]}"; do
        process_combination "$bamfile" "$combine"
    done
        
    echo "===== Finished all combinations for sample $bamfile ====="
done

echo "===== Processing completed ====="
echo "End time: $(date)"
exit 0