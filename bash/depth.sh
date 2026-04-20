#!/bin/bash
set -euo pipefail

# 配置参数
input_path="/data/zhujiahua/WGS/bam/deduped"
output_path="/data/zhujiahua/WGS/bamfile"
software="/home/zhujiahua/software"
depths=(5 10 15 20 25 30 45 60)
bamfile_list="/data/zhujiahua/result/WGS-SVs/bamfile.list"
threads=8  # 总并行任务数 = 样本数 × 深度数

# 检查依赖
command -v parallel >/dev/null 2>&1 || { echo >&2 "需要安装GNU parallel"; exit 1; }
[ ! -f "$bamfile_list" ] && { echo "样本列表不存在: $bamfile_list" >&2; exit 1; }

# 读取样本列表
readarray -t bamfiles < "$bamfile_list"

# 创建输出目录
mkdir -p "$output_path/depth_reports"
for depth in "${depths[@]}"; do
    mkdir -p "$output_path/bam_${depth}x"
done

# 处理单个样本的单个深度
process_depth() {
    local sample="$1"
    local target_depth="$2"
    
    sample=$(echo "$sample" | tr -d '\r\n')
    input_bam="$input_path/${sample}.deduped.bam"
    output_bam="$output_path/bam_${target_depth}x/${sample}.${target_depth}x.bam"
    
    echo "[$(date '+%T')] 处理 $sample 到 ${target_depth}x..."
    
    # 计算原始深度（如果尚未计算）
    if [ ! -f "$output_path/depth_reports/${sample}.original.depth" ]; then
        $software/mosdepth_d4 -n -t 3 --by 500 \
            "$output_path/depth_reports/${sample}.original" \
            "$input_bam"
        tail -1 "$output_path/depth_reports/${sample}.original.mosdepth.summary.txt" | awk '{print $4}' \
            > "$output_path/depth_reports/${sample}.original.depth"
    fi
    
    original_depth=$(cat "$output_path/depth_reports/${sample}.original.depth")
    
    # 检查目标深度是否高于原始深度
    if (( $(echo "$target_depth >= $original_depth" | bc -l) )); then
        echo "跳过: $sample 原始深度 ${original_depth}x 低于目标深度 ${target_depth}x"
        return 0
    fi
    
    # 计算下采样比例
    fraction=$(echo "scale=6; ${target_depth}.03 / $original_depth" | bc)
    
    # 执行下采样
    samtools view -@ 8 -b -s "$fraction" "$input_bam" > "$output_bam" &&
    samtools index -@ 8 "$output_bam"
    
    # 验证结果
    $software/mosdepth_d4 -n -t 3 --by 500 \
        "$output_path/depth_reports/${sample}.${target_depth}x" \
        "$output_bam"
    
    final_depth=$(tail -1 "$output_path/depth_reports/${sample}.${target_depth}x.mosdepth.summary.txt" | awk '{print $4}')
    
    echo "$sample ${target_depth}x: 原始=${original_depth}x, 实际=${final_depth}x"
}

export -f process_depth
export input_path output_path software

# 生成任务列表
task_list=()
for sample in "${bamfiles[@]}"; do
    for depth in "${depths[@]}"; do
        task_list+=("$sample $depth")
    done
done

# 并行执行
printf "%s\n" "${task_list[@]}" | \
parallel -j "$threads" \
    "process_depth {1} {2}" 2>&1

echo "所有任务完成! 结果保存在: $output_path/bam_*x/"