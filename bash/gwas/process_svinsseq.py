#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import gzip
import re
import sys
import argparse

def process_vcf(input_vcf, output_vcf):
    """处理单个VCF文件"""
    # 判断文件类型并打开
    open_func = gzip.open if input_vcf.endswith('.gz') else open
    mode = 'rt' if input_vcf.endswith('.gz') else 'r'
    
    with open_func(input_vcf, mode) as fin, gzip.open(output_vcf, 'wt') as fout:
        svinsseq_header_added = False
        
        for line in fin:
            # 写入头信息
            if line.startswith('#'):
                if line.startswith('##INFO=') and 'ID=SVINSSEQ' in line:
                    continue  # 跳过已有的SVINSSEQ定义
                
                # 在##INFO行之前添加SVINSSEQ定义
                if not svinsseq_header_added and line.startswith('##INFO='):
                    fout.write('##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">\n')
                    svinsseq_header_added = True
                
                fout.write(line)
                continue
            
            # 处理数据行
            fields = line.strip().split('\t')
            if len(fields) < 8:  # 确保有足够的列
                fout.write(line)
                continue
                
            info_field = fields[7]
            
            # 解析INFO字段
            info_dict = {}
            for item in info_field.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = True
            
            # 判断并处理SVINSSEQ
            if 'LEFT_SVINSSEQ' in info_dict and 'RIGHT_SVINSSEQ' in info_dict:
                # 情况1：合并LEFT和RIGHT序列，中间插入100个N
                left_seq = info_dict['LEFT_SVINSSEQ']
                right_seq = info_dict['RIGHT_SVINSSEQ']
                svinsseq = left_seq + 'N' * 100 + right_seq
            else:
                # 情况2：使用ALT列作为SVINSSEQ
                svinsseq = fields[4]  # ALT列
            
            # 添加SVINSSEQ到INFO字段
            if 'SVINSSEQ' in info_dict:
                # 替换现有的SVINSSEQ
                info_field = re.sub(r'SVINSSEQ=[^;]+', f'SVINSSEQ={svinsseq}', info_field)
            else:
                # 添加新的SVINSSEQ
                if info_field and info_field != '.':
                    info_field += f';SVINSSEQ={svinsseq}'
                else:
                    info_field = f'SVINSSEQ={svinsseq}'
            
            # 更新INFO字段
            fields[7] = info_field
            
            # 写入处理后的行
            fout.write('\t'.join(fields) + '\n')
        
        # 如果还没有添加SVINSSEQ的INFO定义，则在文件末尾添加
        if not svinsseq_header_added:
            fout.write('##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Sequence of insertion">\n')

def find_vcf_files_recursive(base_dir):
    """在指定目录下递归查找所有VCF.GZ文件"""
    vcf_files = []
    
    # 使用glob递归查找所有.vcf.gz文件
    pattern = os.path.join(base_dir, "**", "*.vcf.gz")
    vcf_files = glob.glob(pattern, recursive=True)
    
    # 如果没有找到，尝试非递归查找
    if not vcf_files:
        pattern = os.path.join(base_dir, "*.vcf.gz")
        vcf_files = glob.glob(pattern)
    
    return vcf_files

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='处理VCF文件，添加SVINSSEQ信息')
    parser.add_argument('--bamfile-list', required=True, help='BAM文件列表路径')
    parser.add_argument('--manta-dir', required=True, help='Manta VCF文件所在的基础目录')
    parser.add_argument('--output-dir', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    # 设置配置参数
    bamfile_list = args.bamfile_list
    manta_base_dir = args.manta_dir
    output_base_dir = args.output_dir
    
    # 确保输出目录存在
    os.makedirs(output_base_dir, exist_ok=True)
    
    print(f"配置参数:")
    print(f"  BAM文件列表: {bamfile_list}")
    print(f"  Manta目录: {manta_base_dir}")
    print(f"  输出目录: {output_base_dir}")
    print("-" * 50)
    
    # 检查输入目录是否存在
    if not os.path.exists(manta_base_dir):
        print(f"错误: Manta目录不存在: {manta_base_dir}")
        return
    
    # 读取bamfile列表
    try:
        with open(bamfile_list, 'r') as f:
            bamfiles = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"错误: 找不到bamfile列表文件: {bamfile_list}")
        return
    
    print(f"读取到 {len(bamfiles)} 个bamfile")
    
    # 递归查找所有VCF文件
    print("正在递归查找VCF文件...")
    all_vcf_files = find_vcf_files_recursive(manta_base_dir)
    print(f"找到 {len(all_vcf_files)} 个VCF文件")
    
    # 为每个bamfile查找匹配的VCF文件
    processed_count = 0
    for bamfile in bamfiles:
        print(f"处理bamfile: {bamfile}")
        
        # 查找包含bamfile名称的VCF文件
        matching_vcfs = [vcf for vcf in all_vcf_files if bamfile in os.path.basename(vcf)]
        
        if not matching_vcfs:
            print(f"  警告: 未找到包含 '{bamfile}' 的VCF文件")
            continue
        
        print(f"  找到 {len(matching_vcfs)} 个匹配的VCF文件")
        
        for input_vcf in matching_vcfs:
            # 构建输出文件名
            basename = os.path.basename(input_vcf)
            output_vcf = os.path.join(output_base_dir, basename)
            
            print(f"  处理文件: {os.path.basename(input_vcf)} -> {os.path.basename(output_vcf)}")
            try:
                process_vcf(input_vcf, output_vcf)
                processed_count += 1
                print(f"  处理成功")
            except Exception as e:
                print(f"  错误: 处理文件 {input_vcf} 时出错: {e}")
    
    print(f"处理完成! 成功处理 {processed_count} 个VCF文件")
    print(f"输出目录: {output_base_dir}")

if __name__ == "__main__":
    main()
if __name__ == "__main__":
    main()