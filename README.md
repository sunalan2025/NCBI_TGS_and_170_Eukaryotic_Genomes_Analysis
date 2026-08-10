# NCBI 三代测序数据库与 170 个真核生物参考基因组演化分析开源数据集

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![GitHub Repository](https://img.shields.io/badge/GitHub-NCBI__TGS__and__170__Eukaryotic__Genomes__Analysis-brightgreen)](https://github.com/sunalan2025/NCBI_TGS_and_170_Eukaryotic_Genomes_Analysis)

本仓库为《暑期科研见习岗总结报告》的官方配套开源数据与图表仓库，收录了总结报告正文引用的全量 12 张高清图表 (图 3.1 - 图 3.12)、结构化附表、Top 25 核心基因多序列比对超级矩阵及 170 物种系统发育进化树文件。

---

## 📌 研究背景与数据集概述 (Overview)

本开源数据集对应总结报告的两大核心研究内容：

1. **NCBI 第三代测序（TGS）公开基因组数据库规模与发展趋势分析**：
   - 包含基于 NCBI Assembly 数据库 E-utilities API 审计的完整图表与数据表 (`正文图表(Figures)/图3.1 - 图3.5` 及 `统计数据汇总表(stats)/ncbi_tgs_...csv`)，定量展示长读长测序（PacBio HiFi/CLR、Oxford Nanopore ONT）组装数据的历年提交趋势、完成度分层及生态位分布。
2. **170 个代表性真核生物基因组质量评估与系统发育演化分析**：
   - 包含 QUAST 连续性与 BUSCO 完整度评估数据集及对比图表 (`图3.6 - 图3.10`)；
   - 包含用于重构 170 物种系统发育拓扑树的 646 KB 氨基酸多序列比对超级矩阵 (`concatenated_alignment_170.faa`) 与最大似然进化树 (`图3.11`)；
   - 包含 181 个共有 BUSCO 直系同源基因的密码子比对 (PAL2NAL) 与 PAML CodeML Branch-Site 选择压 ($Ka/Ks$) 扫描汇总数据集 (`图3.12`)。

---

## 📂 仓库目录结构 (Repository Structure)

```
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
│   └── Top25_BUSCO_Genes_Table.csv
├── 系统发育树与比对文件/
│   ├── concatenated_alignment_170.faa
│   ├── species_phylogeny_tree_170.pdf
│   └── species_tree_170.nwk
└── 统计数据汇总表(stats)/
    ├── busco_summary_stats.csv
    ├── duplicated_genes_functions.csv
    ├── kaks_results_all.csv
    ├── kaks_species_summary.csv
    ├── ncbi_tgs_assembly_and_domain_stats.csv
    ├── ncbi_tgs_mags_ecological_distribution.csv
    ├── ncbi_tgs_yearly_submission_trends.csv
    ├── positive_selection_branch_site_results.csv
    ├── purifying_selection_wgd.csv
    ├── quast_summary_stats.csv
    ├── repeatmasker_summary_all.csv
    ├── structural_features_all.csv
    ├── taxonomy_summary_stats.csv
    └── tech_summary_stats.csv
```