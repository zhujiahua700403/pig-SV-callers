#!/usr/bin/env python
# coding=utf-8

import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Convert TRA format VCF to BND format while preserving other SVs')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    return parser.parse_args()

def enhance_header(original_header):
    new_header = []
    required_headers = {
        'fileformat': '##fileformat=VCFv4.2',
        'alt_bnd': '##ALT=<ID=BND,Description="Breakend">',
        'info_chr2': '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">',
        'info_mateid': '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">'
    }
    
    # 确保文件格式声明正确
    fileformat_line = required_headers['fileformat']
    for line in original_header:
        if line.startswith("##fileformat"):
            fileformat_line = line.strip()
            break
    
    # 重建头信息结构
    new_header = [fileformat_line]
    
    # 添加其他原始头信息
    for line in original_header:
        stripped_line = line.strip()
        if stripped_line.startswith("##fileformat"):
            continue
        if stripped_line:
            new_header.append(stripped_line)
    
    # 检查并添加必要头信息
    headers_to_add = []
    if not any(line.startswith("##ALT=<ID=BND") for line in new_header):
        headers_to_add.append(required_headers['alt_bnd'])
    if not any(line.startswith("##INFO=<ID=CHR2") for line in new_header):
        headers_to_add.append(required_headers['info_chr2'])
    if not any(line.startswith("##INFO=<ID=MATEID") for line in new_header):
        headers_to_add.append(required_headers['info_mateid'])
    
    # 在contig定义前插入新头信息
    contig_index = next((i for i, line in enumerate(new_header) if line.startswith("##contig")), -1)
    if contig_index != -1:
        for header in reversed(headers_to_add):
            new_header.insert(contig_index, header)
    else:
        new_header.extend(headers_to_add)
    
    return new_header

def process_vcf(input_file, output_file):
    records = []
    
    with open(input_file, 'r') as fin:
        # 读取完整头信息（包含#CHROM行）
        original_header = []
        line = fin.readline().strip()
        while line.startswith('#'):
            original_header.append(line)
            line = fin.readline().strip()
        
        # 处理数据行
        data_lines = []
        if line and not line.startswith('#'):
            data_lines.append(line)
        data_lines.extend([l.strip() for l in fin if l.strip()])
        
        # 增强头信息（已包含#CHROM行）
        enhanced_header = enhance_header(original_header)
        
        # 获取列标题结构
        header_columns = [h for h in original_header if h.startswith('#CHROM')]
        column_count = len(header_columns[0].split('\t')) if header_columns else 8

        # 处理数据行（保持原始顺序）
        for line in data_lines:
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 8:
                print(f"Skipping malformed line: {line}", file=sys.stderr)
                continue
            
            # 补齐列数
            full_fields = fields + [''] * (column_count - len(fields))
            info_dict = {}
            for item in full_fields[7].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
            
            sv_type = info_dict.get('SVTYPE', 'UNKNOWN')
            
            if sv_type == 'TRA':
                chrom = full_fields[0]
                original_pos = int(full_fields[1])
                chr2 = info_dict.get('CHR2')
                original_end = int(info_dict.get('END', original_pos))
                adjusted = False

                # 关键修复：当染色体相同时，确保END >= POS
                if chr2 == chrom:
                    if original_end < original_pos:
                        # 交换POS和END的位置
                        original_pos, original_end = original_end, original_pos
                        adjusted = True
                        print(f"ADJUSTED: Swapped POS and END for {chrom}: new_POS={original_pos} new_END={original_end}")

                pos = original_pos
                end = original_end
                mateid_base = f"{chrom}_{pos}_BND"
                
                # 构造BND1记录
                bnd1 = list(full_fields)
                bnd1[1] = str(pos)
                bnd1[2] = f"{mateid_base}_1"
                bnd1[4] = f"N[{chr2}:{end}["
                bnd1_info = [
                    f"SVTYPE=BND",
                    f"MATEID={mateid_base}_2",
                    f"CHR2={chr2}"
                ]
                if adjusted:
                    bnd1_info.append(f"ORIGINAL_POS={full_fields[1]}")
                    bnd1_info.append(f"ORIGINAL_END={info_dict.get('END')}")
                if 'EVENT' in info_dict:
                    bnd1_info.append(f"EVENT={info_dict['EVENT']}")
                bnd1[7] = ";".join(bnd1_info)
                
                # 构造BND2记录
                bnd2 = [
                    chr2,                           # CHROM
                    str(end),                       # POS
                    f"{mateid_base}_2",             # ID
                    "N",                            # REF
                    f"]{chrom}:{pos}]N",            # ALT
                    bnd1[5],                        # QUAL
                    bnd1[6],                        # FILTER
                    ";".join([                      # INFO
                        f"SVTYPE=BND",
                        f"MATEID={mateid_base}_1",
                        f"CHR2={chrom}"
                    ]),
                    *bnd1[8:]                       # 其他列
                ]
                
                records.append("\t".join(map(str, bnd1)))
                records.append("\t".join(map(str, bnd2)))
            else:
                records.append("\t".join(full_fields))
    
    # 直接写入文件
    with open(output_file, 'w') as fout:
        # 写入完整的头信息（已包含#CHROM行）
        fout.write("\n".join(enhanced_header))
        
        # 直接追加数据记录
        if records:
            fout.write("\n" + "\n".join(records))

def main():
    args = parse_args()
    process_vcf(args.input, args.output)

if __name__ == '__main__':
    main()