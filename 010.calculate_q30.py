import argparse
import os
import numpy as np

def calculate_q30(file_name, phred=33):
    # 逐行读取并选择性提取质量行（内存优化）
    quality_chunks = []
    with open(file_name, 'r') as f:
        while True:
            lines = [f.readline() for _ in range(4)]  # 批量读取4行
            if not lines[0]:  # 文件结束检查
                break
            quality_chunks.append(lines[3].strip())  # 直接提取质量行
    
    # 向量化处理质量分数（速度提升关键）
    qual_bytes = ''.join(quality_chunks).encode()
    q_scores = np.frombuffer(qual_bytes, dtype=np.uint8) - phred
    
    # 快速统计指标
    q30_bases = np.sum(q_scores >= 30)
    total_bases = q_scores.size
    q30_percent = round((q30_bases / total_bases) * 100, 3)
    
    # 写入结果
    output_file = f"{os.path.basename(file_name)}.q30.txt"
    with open(output_file, 'w') as f:
        f.write("File\tQ30_Percent\tQ30_Bases\tTotal_Bases\tRead_Count\n")
        f.write(f"{os.path.basename(file_name)}\t{q30_percent}\t{q30_bases}\t{total_bases}\t{len(quality_chunks)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate Q30 metrics from FASTQ file')
    parser.add_argument('file_name', type=str, help='Path to the FASTQ file')
    parser.add_argument('--phred', type=int, default=33, help='Phred offset (default: 33)')
    args = parser.parse_args()
    calculate_q30(args.file_name, args.phred)
