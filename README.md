# NCBI 三代测序数据库与 170 个真核生物参考基因组演化分析开源数据集

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![GitHub Repository](https://img.shields.io/badge/GitHub-NCBI__TGS__and__170__Eukaryotic__Genomes__Analysis-brightgreen)](https://github.com/sunalan2025/NCBI_TGS_and_170_Eukaryotic_Genomes_Analysis)
[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)

本仓库为《暑期科研见习岗总结报告》的官方配套开源数据、图表与复现分析代码仓库。收录了总结报告正文引用的全量 12 张高清图表（图 3.1 - 图 3.12）、结构化附表、Top 25 核心基因多序列比对超级矩阵、170 物种系统发育拓扑进化树，以及用于计算旁系同源（表 3.1）与直系同源（表 3.2）选择压力（\( K_a/K_s \)）的全流程质控过滤脚本与全量数据集。

---

## 📌 研究背景与数据集概述 (Overview)

本开源数据集对应总结报告的两大核心研究内容：

1. **NCBI 第三代测序（TGS）公开基因组数据库规模与发展趋势分析**：
   - 包含基于 NCBI Assembly 数据库 E-utilities API 审计的完整图表与数据表 (`正文图表(Figures)/图3.1 - 图3.5` 及 `统计数据汇总表(stats)/ncbi_tgs_...csv`)，定量展示长读长测序（PacBio HiFi/CLR、Oxford Nanopore ONT）组装数据的历年提交趋势、完成度分层及生态位分布。
2. **170 个代表性真核生物基因组质量评估与系统发育演化分析**：
   - 包含 QUAST 连续性与 BUSCO 完整度评估数据集及对比图表 (`图3.6 - 图3.10`)；
   - 包含用于重构 170 物种系统发育拓扑树的 646 KB 氨基酸多序列比对超级矩阵 (`concatenated_alignment_170.faa`)、涵盖 43 个代表性真核生物的系统发育拓扑骨架树 (`图3.11`) 及 170 物种全景环形进化树 (`species_phylogeny_tree_170_circular.pdf`)；
   - 包含基于 NCBI RefSeq 官方审校 CDS 最长转录本代表序列的代表性物种旁系同源基因 Ka/Ks 选择压力全量数据集 (`表3.1_代表性物种旁系同源重复基因KaKs选择压力统计汇总表.csv` 及 `kaks_species_summary.csv`, `kaks_results_all.csv`)；
   - 包含 181 个共有 BUSCO 单拷贝核心直系同源基因的密码子对齐 (PAL2NAL) 与基于 Nei-Gojobori (NG86) 双重未饱和过滤 (\( K_a < 0.75, 0 < K_s < 0.75 \)) 的各大分类纲探索性试算数据集 (`图3.12`、`181_BUSCO_Orthologs_KaKs_Summary.csv`、`表3.2_各大纲类群BUSCO直系同源基因KaKs选择压力统计汇总表.csv` 及 `real_class_ortholog_kaks_summary.csv`)。

---

## 🔬 两套核心质控与过滤逻辑说明 (Core Filtering Logics)

为确保科研结论的严谨性与真实性，避免传统生信流程中常见的“转录本假性膨胀”与“深部分歧突变饱和假象”，分析管线设立了两套核心质控过滤逻辑：

### 1. 转录本最长异构体提取去冗余 (`step1_longest_cds.py`)
- **科学背景与痛点**：
  在真核生物（尤其是昆虫纲、硬骨鱼类与哺乳纲）基因组注释中，一个蛋白编码基因位点通常对应 3~10 个可变剪接异构体（Alternative Splicing Isoforms）。若直接使用原始 CDS 进行全对全自比对（All-against-all BLASTP），同一基因位点内部的不同转录本之间会产生大量相似度极高但毫无演化分化意义的“假对”，导致候选旁系同源对数呈爆炸式虚假增长（如东方蜜蜂曾虚假输出上万对）。
- **执行过滤标准**：
  1. **位点聚类**：精确解析 RefSeq FASTA 注释中的 `[gene=...]` 结构化标签，按位点聚合；
  2. **最长代表序列提取**：每个基因位点严格仅提取最长的一条 CDS 代表序列；
  3. **读码框生物学质控**：剔除长度非 3 的倍数（`len(seq) % 3 != 0`）序列，剥除末端终止密码子；
  4. **短肽与假基因剔除**：剔除编码长度小于 90 bp（< 30 个氨基酸）或翻译序列内部含有未成熟终止密码子（`*`）的异常序列。

---

### 2. 双重未饱和过滤与年轻旁系对筛选 (`step2_filter_young_paralogs.py`)
- **科学背景与痛点**：
  在利用 Nei-Gojobori (NG86) 算法计算同义替换率 \( K_s \) 时，根据 Jukes-Cantor 公式修正同义距离：
  
  $$
  K_s = -\frac{3}{4} \ln\left(1 - \frac{4}{3} p_S\right)
  $$
  
  当同义位点差异率 \( p_S \ge 0.75 \) 时，对数真数小于等于零，在实数域内无解（突变方差趋向无穷大）。对于真核生物体内古老加倍形成的多基因家族（如蜜蜂体内的气味受体 OR 家族、细胞色素 P450 家族），其分化年代极其久远，**超过 95% 的旁系对已处于同义突变深度饱和区**。若不加甄别地全盘平均，将彻底扭曲选择压力指标。
- **执行过滤标准**：
  1. **Jukes-Cantor 数学可解性约束**：严格要求 \( p_S < 0.75 \)；
  2. **同义未饱和区筛选**：施加 \( 0 < K_s < 0.75 \) 门槛：
     - 剔除 \( K_s = 0 \)：排除序列完全相同、未分化的等位基因组装未坍缩片段或技术扩增假象；
     - 剔除 \( K_s \ge 0.75 \)：排除突变已深度饱和、无法准确估计演化时间的古老重复对；
     - **保留 \( 0 < K_s < 0.75 \)**：精确提取在近期演化历史中活跃产生、能够灵敏反映自然选择净化与适应性创新的**年轻有效旁系对（Young Paralogs）**。

---

## ⚡ 一键复现指南 (One-Click Reproduction Guide)

克隆本仓库后，可通过提供的标准 Python 脚本一键复现报告正文【表 3.1】与【表 3.2】的全量统计指标。

### 1. 环境依赖
仅需标准 Python 3.8+ 及 pandas：
```bash
pip install pandas
```
（注：若需运行 `step1_longest_cds.py` 解析原始 FASTA，需额外安装 `biopython`）。

### 2. 一键复现所有表格 (Tables 3.1 & 3.2)
在仓库根目录下运行：
```bash
python scripts/reproduce_all_tables.py
```

### 3. 分项复现命令
- **复现表 3.1（代表性物种旁系同源基因 Ka/Ks 选择压力统计汇总）**：
  ```bash
  python scripts/step2_filter_young_paralogs.py
  ```
  终端将直接输出 7 个代表物种的候选全集、年轻有效对、Ka、Ks、Ka/Ks 均值/中位数及正选择/纯化选择比例，精确对应正文表 3.1 与第 19 页数值（如西方蜜蜂纯化选择 94.92%，东方蜜蜂纯化选择 93.43%）。

- **复现表 3.2（170 真核生物各大分类纲直系同源基因 Ka/Ks 探索性试算汇总）**：
  ```bash
  python scripts/step3_reproduce_table3_2.py
  ```
  终端将输出被子植物、鸟纲、昆虫纲、辐鳍鱼纲、两栖纲、哺乳纲及真菌类的完整选择压统计与纯化选择比例（96.64% ~ 100%）。

---

## 📂 仓库完整目录结构 (Repository Structure)

```text
NCBI_TGS_and_170_Eukaryotic_Genomes_Analysis/
├── README.md                                             # 仓库综合文档与核心过滤逻辑说明
├── scripts/                                              # 官方复现与数据质控脚本目录
│   ├── step1_longest_cds.py                             # 脚本1：RefSeq CDS 最长转录本提取与去冗余
│   ├── step2_filter_young_paralogs.py                   # 脚本2：双重未饱和过滤 (pS<0.75 & 0<Ks<0.75) 与表3.1复现
│   ├── step3_reproduce_table3_2.py                      # 脚本3：各大分类纲 BUSCO 直系同源选择压与表3.2复现
│   └── reproduce_all_tables.py                          # 一键复现主程序：自动复现并核验表3.1与表3.2
├── 正文图表(Figures)/                                     # 正文全量 12 张出版级高清图表
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
│   ├── 图3.11_170个真核物种最大似然系统发育拓扑进化树.png  # 涵盖43个代表物种骨架树 (Skeleton Tree)
│   └── 图3.12_181个BUSCO直系同源基因KaKs选择压分布与各大类群演化对比热图.png
├── 正文附表与基因清单/                                    # 报告核心结构化附表数据集
│   ├── 181_BUSCO_Orthologs_KaKs_Summary.csv              # 181 个单拷贝直系同源基因 Ka/Ks 明细表
│   ├── NCBI_TGS_Yearly_Submission_Stats_Table.csv        # NCBI 2013-2026 历年提交数量统计表
│   ├── Top25_BUSCO_Genes_Table.csv                       # 建树核心 Top 25 超保守基因注释清单
│   ├── real_class_ortholog_kaks_summary.csv              # 各大纲 BUSCO 直系同源选择压汇总
│   ├── 表3.1_代表性物种旁系同源重复基因KaKs选择压力统计汇总表.csv
│   └── 表3.2_各大纲类群BUSCO直系同源基因KaKs选择压力统计汇总表.csv
├── 系统发育树与比对文件/                                  # 进化树与多序列比对原始超矩阵
│   ├── concatenated_alignment_170.faa                    # 646 KB 氨基酸多序列比对超级矩阵 (25个核心基因)
│   ├── species_tree_170.nwk                              # FastTree 构建的 170 物种 Newick 格式树文件
│   ├── species_phylogeny_tree_170.pdf                    # 170 物种系统发育树矢量图
│   ├── species_phylogeny_tree_170.svg                    # 170 物种系统发育树 SVG 矢量图
│   ├── species_phylogeny_tree_170_circular.pdf           # 170 物种 360° 全景环形最大似然进化树
│   ├── species_phylogeny_tree_backbone.pdf               # 骨架拓扑树矢量 PDF
│   └── species_phylogeny_tree_backbone.svg               # 骨架拓扑树 SVG
└── 统计数据汇总表(stats)/                                 # 汇总统计数据表目录
    ├── busco_summary_stats.csv
    ├── duplicated_genes_functions.csv
    ├── kaks_results_all.csv                              # 7 代表物种全量 122,098 对计算原始数据 (14.2 MB)
    ├── kaks_species_summary.csv                          # 7 代表物种选择压分项汇总表
    ├── ncbi_tgs_assembly_and_domain_stats.csv
    ├── ncbi_tgs_mags_ecological_distribution.csv
    ├── ncbi_tgs_yearly_submission_trends.csv
    ├── quast_summary_stats.csv
    ├── real_class_ortholog_kaks_summary.csv
    ├── repeatmasker_summary_all.csv                      # 170 物种 RepeatMasker 转座子注释全集
    ├── structural_features_all.csv
    ├── taxonomy_summary_stats.csv
    └── tech_summary_stats.csv
```
