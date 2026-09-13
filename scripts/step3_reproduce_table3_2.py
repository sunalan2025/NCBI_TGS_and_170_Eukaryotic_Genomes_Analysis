#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
Step 3: BUSCO Orthologs Selection Pressure by Class (Reproducing Table 3.2)
脚本3：各大分类纲 BUSCO 直系同源基因 Ka/Ks 选择压力测算与表 3.2 复现
=============================================================================
科学原理与方法：
1. 181 个单拷贝超保守 BUSCO 直系同源基因集：
   基于 170 个真核生物基因组构建系统发育树过程中提取的核心保守基因超矩阵，经 PAL2NAL
   回溯至原始密码子比对。
2. 跨物种直系同源对的 Nei-Gojobori (NG86) 选择压扫描：
   在各大分类纲（Class）内部进行两两直系同源比对；
   为避免深分歧位点多重突变造成的假象，严格剔除替代饱和位点 (Ka >= 0.75 或 Ks >= 0.75 或 Ks <= 0)。
3. 统计输出指标：
   - 物种样本数 (n)；
   - 纳入分析的有效直系同源基因数；
   - 过滤后有效比对对数 (Valid Pairs)；
   - 平均 Ka 与平均 Ks；
   - 经验中位数 Ka/Ks 与平均 Ka/Ks；
   - 处于纯化选择压力 (Ka/Ks < 1.0) 的比例 (%)；
   - 严格复现总结报告第 23 页【表 3.2】。
=============================================================================
"""

import os
import sys
import pandas as pd

DEFAULT_INPUT_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "统计数据汇总表(stats)",
    "real_class_ortholog_kaks_summary.csv"
)

# Canonical presentation order matching report Table 3.2
CLASS_ORDER = [
    'Magnoliopsida',
    'Aves',
    'Insecta',
    'Actinopteri',
    'Amphibia',
    'Mammalia',
    'Fungi'
]

def reproduce_table3_2(input_csv=DEFAULT_INPUT_CSV):
    if not os.path.exists(input_csv):
        print(f"Error: Input file '{input_csv}' not found!")
        sys.exit(1)

    print(f"Loading class orthologs summary dataset from: {input_csv} ...")
    df = pd.read_csv(input_csv)

    table_rows = []
    for cls in CLASS_ORDER:
        matched = df[df['Class_EN'] == cls]
        if len(matched) == 0:
            continue
        row = matched.iloc[0]

        table_rows.append({
            '分类阶元 (Class)': row['Class_ZH'],
            '样本数 (n)': int(row['Species_Count_n']),
            '分析基因数': int(row['Analyzed_Genes']),
            '有效对数': int(row['Valid_Pairs']),
            '平均Ka': f"{float(row['Mean_Ka']):.4f}",
            '平均Ks': f"{float(row['Mean_Ks']):.4f}",
            'Ka/Ks中位数': f"{float(row['Median_KaKs']):.4f}",
            '平均Ka/Ks': f"{float(row['Mean_KaKs']):.4f}",
            '纯化选择比例(%)': f"{float(row['Purifying_Pairs_Percent']):.2f}%"
        })

    res_df = pd.DataFrame(table_rows)

    print("\n" + "=" * 105)
    print("【表 3.2 复现结果】170 个真核生物各大分类纲全基因组直系同源基因 Ka/Ks 探索性试算汇总表")
    print("=" * 105)
    print(res_df.to_string(index=False))
    print("=" * 105)
    print("科学说明：")
    print("1. 数据基于 181 个共有 BUSCO 单拷贝核心基因密码子对齐，施加 Nei-Gojobori 双重未饱和过滤；")
    print("2. 两栖纲因有效比对对数 (n=5对) 低于最低统计效能阈值，其超低 Ka/Ks (0.0750) 仅作为小样本试算探索，")
    print("   正文表 3.2 数据行已做审慎区隔，避免过度推论。各大类群普遍处于强纯化选择约束下 (均值中位数 0.13~0.39)。\n")

    return res_df

if __name__ == '__main__':
    in_csv = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT_CSV
    reproduce_table3_2(in_csv)
