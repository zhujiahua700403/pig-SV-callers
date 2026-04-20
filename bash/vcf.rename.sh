#!/bin/bash

# 定义文件路径
input_dir="/home/zhujiahua/base/genotyped"

# 进入指定目录
cd "$input_dir" || { echo "目录不存在"; exit 1; }

# 遍历所有 *.sniffles.vcf.gz 文件
for file in *.sniffles.vcf.gz; do
    # 检查文件是否存在
    if [[ -e "$file" ]]; then
        # 构造新的文件名
        new_file="${file/.sniffles.vcf.gz/.genotyped.vcf.gz}"

        # 重命名文件
        mv "$file" "$new_file"
        echo "重命名: $file -> $new_file"

        # 使用 tabix 重新建立索引
        tabix -p vcf "$new_file"
        echo "索引建立: $new_file"
    else
        echo "未找到匹配的文件"
    fi
done
rm *.sniffles.vcf.gz.tbi