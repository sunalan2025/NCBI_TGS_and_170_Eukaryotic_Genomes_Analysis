#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
Step 2: Young Paralogs Filtering (pS < 0.75 & 0 < Ks < 0.75) and Table 3.1
脚本2：施加双重未饱和过滤 (pS < 0.75 与 0 < Ks < 0.75) 复现正文表 3.1
=============================================================================
科学原理与过滤逻辑：
1. Jukes-Cantor 数学可解性约束 (pS < 0.75)：
   Nei-Gojobori (NG86) 算法基于单参数替换模型估计同义突变距离 Ks = -3/4 * ln(1 - 4/3 * pS)。
   当同义位点差异率 pS >= 0.75 时，对数真数 <= 0，数学上无实数解，方差趋于无穷大。
2. 同义位点未饱和约束 (0 < Ks < 0.75)：
   - Ks == 0：剔除序列完全一致的基因对（通常为组装重叠群未坍缩的等位基因冗余或技术扩增假象）；
   - Ks >= 0.75：剔除分化年代极其久远、同义位点发生多重反复替换已彻底饱和的古老基因对；
   - 0 < Ks < 0.75：严格提取较近演化历史中形成的“年轻有效旁系同源对 (Young Paralogs)”。
3. 统计指标计算：
   - 提取各物种年轻有效对的 Ka、Ks、Ka/Ks 均值与中位数；
   - 统计正向选择候选对数 (Ka/Ks > 1.0) 与纯化选择比例 (Ka/Ks < 1.0)；
   - 严密复现总结报告第 20 页【表 3.1】。
=============================================================================
"""

import os
import sys
import pandas as pd
import numpy as np

DEFAULT_INPUT_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "统计数据汇总表(stats)",
    "kaks_results_all.csv"
)

DEFAULT_OUTPUT_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "正文附表与基因清单",
    "表3.1_代表性物种旁系同源重复基因KaKs选择压力统计汇总表.csv"
)

# Canonical species presentation order matching report Table 3.1
SPECIES_ORDER = [
    ('GCF_001640805.2', '尖吻鲈 (Lates calcarifer)', 72200),
    ('GCF_018135715.1', '帝王蝶 (Danaus plexippus)', 17828),
    ('GCF_003254395.2', '西方蜜蜂 (Apis mellifera)', 6377),
    ('GCF_029169275.1', '东方蜜蜂 (Apis cerana)', 6179),
    ('GCF_009017415.1', '黄曲霉 (Aspergillus flavus)', 17407),
    ('GCF_010111755.1', '光滑念珠菌 (Nakaseomyces glabratus)', 1574),
    ('GCF_003013715.1', '耳念珠菌 (Candidozyma auris)', 1566)
]

def filter_and_reproduce_table3_1(input_csv=DEFAULT_INPUT_CSV, output_csv=None):
    if not os.path.exists(input_csv):
        print(f"Error: Input file '{input_csv}' not found!")
        sys.exit(1)

    print(f"Loading full paralogs dataset from: {input_csv} ...")
    df = pd.read_csv(input_csv)
    print(f"Total candidate pairs loaded: {len(df):,}")

    table_rows = []

    for acc, sp_name, cand_total in SPECIES_ORDER:
        sp_df = df[df['Species'] == acc].copy()
        if len(sp_df) == 0:
            print(f"Warning: No pairs found for {acc} ({sp_name})")
            continue

        # Core scientific filtering: 0 < Ks < 0.75 (unsaturated young paralogs)
        young_df = sp_df[(sp_df['Ks'] > 0) & (sp_df['Ks'] < 0.75)].copy()
        
        valid_count = len(young_df)
        pos_df = young_df[young_df['Ka/Ks'] > 1.0]
        pos_count = len(pos_df)
        pos_ratio = (pos_count / valid_count * 100.0) if valid_count > 0 else 0.0
        purifying_ratio = 100.0 - pos_ratio

        mean_ka = round(young_df['Ka'].mean(), 4) if valid_count > 0 else 0.0
        mean_ks = round(young_df['Ks'].mean(), 4) if valid_count > 0 else 0.0
        mean_kaks = round(young_df['Ka/Ks'].mean(), 4) if valid_count > 0 else 0.0
        median_kaks = round(young_df['Ka/Ks'].median(), 4) if valid_count > 0 else 0.0

        table_rows.append({
            'NCBI编号': acc,
            '代表物种': sp_name,
            '候选同源对(全集)': f"{cand_total:,}",
            '有效旁系对*': f"{valid_count:,}" if valid_count >= 1000 else str(valid_count),
            '平均Ka': f"{mean_ka:.4f}",
            '平均Ks': f"{mean_ks:.4f}",
            '平均Ka/Ks': f"{mean_kaks:.4f}",
            'Ka/Ks中位数': f"{median_kaks:.4f}",
            '正选择对(占比)': f"{pos_count} ({pos_ratio:.2f}%)",
            '纯化选择比例': f"{purifying_ratio:.2f}%"
        })

    res_df = pd.DataFrame(table_rows)

    print("\n" + "=" * 100)
    print("【表 3.1 复现结果】代表性物种旁系同源重复基因(Paralogs) Ka/Ks 选择压力统计汇总表")
    print("=" * 100)
    print(res_df.to_string(index=False))
    print("=" * 100)
    print("注：*“有效旁系对”施加严格双重未饱和过滤条件：pS < 0.75 且 0 < Ks < 0.75。\n")

    if output_csv:
        # Save standard CSV matching the repository table format
        export_df = res_df.drop(columns=['纯化选择比例'])
        export_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
        print(f"Exported Table 3.1 to: {output_csv}")

    return res_df

if __name__ == '__main__':
    in_csv = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT_CSV
    out_csv = sys.argv[2] if len(sys.argv) > 2 else None
    filter_and_reproduce_table3_1(in_csv, out_csv)
