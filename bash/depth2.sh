#!/bin/bash

# 配置部分（用户需修改这些变量）
input_dir="/data/dingrongrong/PGRPv2/bam"    # 输入BAM文件目录
output_dir="/data/zhujiahua/LGS/bam/WGS/Berkshire"           # 输出目录
software_dir="/home/zhujiahua/software"        # 软件安装目录
target_depth=25.03                               # 目标测序深度
sample_list="/data/zhujiahua/LGS/bam/WGS/bamfiles.list"  # 样本列表文件

# 检查必要文件和目录
[ ! -f "$sample_list" ] && { echo "错误：样本列表文件 $sample_list 不存在"; exit 1; }
[ ! -d "$input_dir" ] && { echo "错误：输入目录 $input_dir 不存在"; exit 1; }
mkdir -p "${output_dir}/depth_report" "${output_dir}/depth"          # 创建深度报告目录

# 检查软件是否存在
[ ! -f "${software_dir}/mosdepth_d4" ] && { echo "错误：mosdepth_d4 未找到"; exit 1; }

# 主处理流程
while IFS= read -r sample; do
    # 清理样本名（去除特殊字符）
    sample=$(echo "$sample" | tr -d '\n\r')
    input_bam="${input_dir}/${sample}.deduped.bam"
    
    # 检查输入文件是否存在
    [ ! -f "$input_bam" ] && { echo "警告：${sample} 的BAM文件不存在，跳过..."; continue; }
    
    echo "正在处理样本: $sample"
    
    # 1. 计算原始深度
    echo "  ▶ 计算原始测序深度..."
    ${software_dir}/mosdepth_d4 -n -t 3 --by 500 \
        "${output_dir}/depth_report/${sample}.original" \
        "$input_bam" || { echo "错误: mosdepth 计算深度失败"; continue; }
    original_depth=$(awk 'END{print $4}' "${output_dir}/depth_report/${sample}.original.mosdepth.summary.txt")
    echo "  ✅ 原始深度: ${original_depth}x"
    
    # 2. 下采样到目标深度
    echo "  ▶ 下采样至 ${target_depth}x..."
    fraction=$(echo "scale=5; $target_depth / $original_depth" | bc)
    output_bam="${output_dir}/${sample}.deduped.bam"
    
    samtools view -@ 6 -b -s "$fraction" "$input_bam" > "$output_bam" && \
    samtools index "$output_bam" || { echo "错误: samtools下采样失败"; continue; }
    
    # 3. 计算最终深度
    echo "  ▶ 计算下采样后深度..."
    ${software_dir}/mosdepth_d4 -n -t 3 --by 500 \
        "${output_dir}/depth_report/${sample}.downsampled" \
        "$output_bam"
    final_depth=$(awk 'END{print $4}' "${output_dir}/depth_report/${sample}.downsampled.mosdepth.summary.txt")
    
    # 4. 生成处理报告
    echo "样本 ${sample} 处理完成: 原始深度=${original_depth}x, 目标深度=${target_depth}x, 实际深度=${final_depth}x"
    echo "输出文件: ${output_bam}"
    echo "--------------------------------------------------"
done < "$sample_list"

echo "所有样本处理完成！结果保存在: ${output_dir}"