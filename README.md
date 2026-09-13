# NCBI 三代测序数据库与 170 个真核生物参考基因组演化分析开源数据集

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![GitHub Repository](https://img.shields.io/badge/GitHub-NCBI__TGS__and__170__Eukaryotic__Genomes__Analysis-brightgreen)](https://github.com/sunalan2025/NCBI_TGS_and_170_Eukaryotic_Genomes_Analysis)

本仓库为《暑期科研见习岗总结报告》的官方配套开源数据与图表仓库，收录了总结报告正文引用的全量 12 张高清图表（图 3.1 - 图 3.12）、结构化附表、Top 25 核心基因多序列比对超级矩阵、170 物种系统发育拓扑进化树及全基因组选择压力分析全量数据集。

---

## 📌 研究背景与数据集概述 (Overview)

本开源数据集对应总结报告的两大核心研究内容：

1. **NCBI 第三代测序（TGS）公开基因组数据库规模与发展趋势分析**：
   - 包含基于 NCBI Assembly 数据库 E-utilities API 审计的完整图表与数据表 (`正文图表(Figures)/图3.1 - 图3.5` 及 `统计数据汇总表(stats)/ncbi_tgs_...csv`)，定量展示长读长测序（PacBio HiFi/CLR、Oxford Nanopore ONT）组装数据的历年提交趋势、完成度分层及生态位分布。
2. **170 个代表性真核生物基因组质量评估与系统发育演化分析**：
   - 包含 QUAST 连续性与 BUSCO 完整度评估数据集及对比图表 (`图3.6 - 图3.10`)；
   - 包含用于重构 170 物种系统发育拓扑树的 646 KB 氨基酸多序列比对超级矩阵 (`concatenated_alignment_170.faa`)、涵盖 43 个代表性真核生物的系统发育拓扑骨架树 (`图3.11`) 及 170 物种全景环形进化树 (`species_phylogeny_tree_170_circular.pdf`)；
   - 包含基于 NCBI RefSeq 官方审校 CDS 最长转录本代表序列的代表性物种旁系同源基因 Ka/Ks 选择压力数据集 (`表3.1_代表性物种旁系同源重复基因KaKs选择压力统计汇总表.csv` 及 `kaks_species_summary.csv`, `kaks_results_all.csv`)；
   - 包含 181 个共有 BUSCO 单拷贝核心直系同源基因的密码子对齐 (PAL2NAL) 与基于 Nei-Gojobori (NG86) 双重未饱和过滤 (\( K_a < 0.75, 0 < K_s < 0.75 \)) 的各大分类纲探索性试算数据集 (`图3.12`、`181_BUSCO_Orthologs_KaKs_Summary.csv`、`表3.2_各大纲类群BUSCO直系同源基因KaKs选择压力统计汇总表.csv` 及 `real_class_ortholog_kaks_summary.csv`)。

---

## 📂 仓库目录结构 (Repository Structure)

```text
NCBI_TGS_and_170_Eukaryotic_Genomes_Analysis/
├── README.md
├── 正文图表(Figures)/
│   ├── 图3.1_NCBI三代测序平台与联合组装比例分布.png
│   ├── 图3.2_2013-2026年三代测序基因组历年提交数量趋势.png
│   ├── 图3.3_全量三代基因组组装完成度级别分布.png
│   ├── 图3.4_三代测序基因组生物学域界及真核生物子类群分布.png
│   ├── 图3.5_环境宏基因组（MAGs）生态位分布.png
│   ├── 图3.6_不同测序技术Contig_N50连续度分布对比.png
│   ├── 图3.7_不同测序技术BUSCO完整度分布对比.png
│   ├── 图3.8_170个基因组评估指标主成分分析(PCA)聚类图.png
│   ├── 图3.9_QUAST与BUSCO评估指标相关性热图.png
│   ├── 图3.10_不同分类群的BUSCO重复基因比例分布.png
│   ├── 图3.11_170个真核物种最大似然系统发育拓扑进化树.png
│   └── 图3.12_181个BUSCO直系同源基因KaKs选择压分布与各大类群演化对比热图.png
├── 正文附表与基因清单/
│   ├── 181_BUSCO_Orthologs_KaKs_Summary.csv
│   ├── NCBI_TGS_Yearly_Submission_Stats_Table.csv
│   ├── Top25_BUSCO_Genes_Table.csv
│   ├── real_class_ortholog_kaks_summary.csv
│   ├── 表3.1_代表性物种旁系同源重复基因KaKs选择压力统计汇总表.csv
│   └── 表3.2_各大纲类群BUSCO直系同源基因KaKs选择压力统计汇总表.csv
├── 系统发育树与比对文件/
│   ├── concatenated_alignment_170.faa
│   ├── species_phylogeny_tree_170.pdf
│   ├── species_phylogeny_tree_170.svg
│   ├── species_phylogeny_tree_170_circular.pdf
│   ├── species_phylogeny_tree_backbone.pdf
│   ├── species_phylogeny_tree_backbone.svg
│   └── species_tree_170.nwk
└── 统计数据汇总表(stats)/
    ├── busco_summary_stats.csv
    ├── duplicated_genes_functions.csv
    ├── kaks_results_all.csv
    ├── kaks_species_summary.csv
    ├── ncbi_tgs_assembly_and_domain_stats.csv
    ├── ncbi_tgs_mags_ecological_distribution.csv
    ├── ncbi_tgs_yearly_submission_trends.csv
    ├── quast_summary_stats.csv
    ├── real_class_ortholog_kaks_summary.csv
    ├── repeatmasker_summary_all.csv
    ├── structural_features_all.csv
    ├── taxonomy_summary_stats.csv
    └── tech_summary_stats.csv
```

---

## 🔬 数据来源与统计说明 (Methods & Data Notes)

1. **旁系同源基因过滤标准**：
   - 提取自 NCBI RefSeq 官方蛋白质/CDS 数据库，每个基因位点严格仅保留最长转录本代表序列，剔除可变剪接异构体冗余；
   - 采用 BLASTP 全对全比对 (\( E < 10^{-5} \))，严格剔除自比对 (\( \text{ID}_1 = \text{ID}_2 \)) 与镜像重复，并通过双向比对对齐；
   - 采用 Nei-Gojobori (NG86) 算法计算 \( K_a \)、\( K_s \) 及 \( K_a/K_s \)，并设立 \( K_s < 0.75 \) 严格未饱和阈值以识别近期演化重复对。
2. **BUSCO 直系同源基因选择压力扫描**：
   - 基于 170 物种建树中提取的 181 个单拷贝直系同源基因集；
   - 采用 PAL2NAL 将氨基酸多序列比对回溯至原始密码子比对，计算各分类群物种对的比值；
   - 严格剔除饱和替代位点 (\( K_a \ge 0.75 \) 或 \( K_s \ge 0.75 \) 或 \( K_s = 0 \))，保证各大类群演化保守性评估的稳健性。
