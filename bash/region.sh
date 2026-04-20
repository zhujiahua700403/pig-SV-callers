#!/bin/bash

# 基因组文件
input_path="/home/zhujiahua/genome"
output_path="/home/zhujiahua/test"
temp_path="/home/zhujiahua/temp"
bamfile="Sus_scrofa.Sscrofa11.1"

# 提取染色体信息
faidx \
  $input_path/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa \
-i chromsizes | grep -v 'M' \
  > $temp_path/$bamfile.size.genome

# 生成窗口文件
bedtools makewindows \
  -g $temp_path/$bamfile.size.genome \
  -w 50 > \
  $temp_path/$bamfile.windows.bed || { echo "Error generating windows"; exit 1; }

# 检查窗口文件内容
echo "生成的窗口文件:"
cat $temp_path/$bamfile.windows.bed

# 计算 GC 含量
bedtools nuc \
  -fi $input_path/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa \
  -bed $temp_path/$bamfile.windows.bed > \
  $temp_path/$bamfile.Sus_scrofa.Sscrofa11.1.gcstat.txt || { echo "Error calculating GC content"; exit 1; }

# 检查 GC 含量结果
echo "GC 含量统计结果:"
cat $temp_path/$bamfile.gcstat.txt

# 提取 GC 含量高于 80% 的区域
awk 'NR==1 || $5 > 0.80' \
  $temp_path/$bamfile.Sus_scrofa.Sscrofa11.1.gcstat.txt > \
  $temp_path/$bamfile.high_gc_regions.bed



# 定义重复序列分类文件
declare -A repeat_beds=(
    ["DNA"]="$output_dir/DNA.bed"
    ["LINE"]="$output_dir/LINE.bed"
    ["SINE"]="$output_dir/SINE.bed"
    ["LTR"]="$output_dir/LTR.bed"
    ["Low_complexity"]="$output_dir/Low_complexity.bed"
    ["Simple_repeat"]="$output_dir/Simple_repeat.bed"
    ["Other"]="$output_dir/Other.bed"
)

# 清空输出文件
for bed_file in "${repeat_beds[@]}"; do
    > "$bed_file" || {
        log "错误: 无法清空文件 $bed_file"
        exit 1
    }
done

# 处理RepeatMasker文件并分类
log "开始处理RepeatMasker文件..."
while read -r line; do
    # 跳过注释行和空行
    [[ "$line" =~ ^# ]] && continue
    [[ -z "$line" ]] && continue
    
    # 提取各列并处理染色体名称（去掉"chr"前缀）
    chr=$(echo "$line" | awk '{print $5}' | sed 's/^chr//')
    start=$(echo "$line" | awk '{print $6}')
    end=$(echo "$line" | awk '{print $7}')
    full_type=$(echo "$line" | awk '{print $11}')

    # 根据full_type设置type
    case "$full_type" in
        SINE/*) type="SINE" ;;
        LINE/*) type="LINE" ;;
        LTR/*) type="LTR" ;;
        DNA/*) type="DNA" ;;
        "Simple_repeat") type="Simple_repeat" ;;
        "Low_complexity") type="Low_complexity" ;;
        *) type="Other" ;;
    esac

    # 写入对应的bed文件
    echo -e "$chr\t$start\t$end\t$type" >> "${repeat_beds[$type]}" || {
        log "错误: 无法写入 ${repeat_beds[$type]}"
        exit 1
    }
done < /home/zhujiahua/genome/susScr11.fa.out

# 删除临时解压文件
rm /home/zhujiahua/genome/susScr11.fa.out || {
    log "警告: 无法删除临时文件 /home/zhujiahua/genome/susScr11.fa.out"
}

log "RepeatMasker分类完成，结果保存在 $output_dir 目录下"
