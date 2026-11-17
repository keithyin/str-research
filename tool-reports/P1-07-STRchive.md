# STRchive

> **优先级**：P1（重要）｜ **类别**：致病 STR 目录（数据资源）｜ **维护方**：Dashnow lab（Harriet Dashnow，亦为 TRGT 作者之一）｜ **状态**：活跃，2.0.0（73 位点），月度更新
> **一句话**：动态维护的致病 STR 位点数据库（73 个位点，含 motif/致病范围/人群频率/发病年龄），TRGT 已弃用自带 BED、官方改推 STRchive——你的目录层"事实来源"。

## 是什么

STRchive 是社区维护的 STR 疾病位点资源：GitHub 版本控制 + 网站浏览（strchive.org），数据 CC BY 4.0、代码 MIT。输出机器可读 JSON 与 genotyper-ready BED（hg37/hg38/T2T-chm13）。

- 网站：https://strchive.org ｜ 仓库：https://github.com/dashnowlab/STRchive （~30★）
- 论文：*STRchive: a dynamic resource detailing population-level and locus-specific insights at tandem repeat disease loci*, **Genome Medicine 2025**（Hiatt et al.）https://doi.org/10.1186/s13073-025-01454-4

## 解决了什么问题

1. **致病 STR 信息散落且过时**：OMIM/ClinVar/GeneReviews 缺 TR 专属细节或滞后数年；STRipy 只收 65 个、gnomAD 只收 60 个位点，且都缺流行病学字段。
2. **工具自带的目录会"烂掉"**：TRGT 原始打包的 `pathogenic_repeats.hg38.bed` 已被官方标记 deprecated——静态目录跟不上新位点发现（如近年新增的 ABCD3、FGF14、MUC1、TAF1、THAP11、ZFHX3、ZNF713）。
3. **目录定义影响分型质量**：gnomAD 的位点定义是为 ExpansionHunter（短读）校准的，偏窄；长读需要更宽的位点定义。

## 如何解决的

### 1. 数据模型：一个位点 = 一份结构化档案
每位点字段：基因组位置（hg37/hg38/T2T-chm13 三套 BED）、canonical 与 genotyped motif（含临床分类）、正常/中间/致病等位基因大小范围、遗传方式、发病年龄范围（典型+极端）、疾病流行率估计、独立观察计数、结构/中断说明、关键文献引用、到 OMIM/GeneReviews/gnomAD/STRipy/UCSC 的链接。

### 2. 更新机制（"动态"的实现）
- **自动化文献监测**：R/easyPubMed + biomaRt 基因同义词展开，**月度**跑 PubMed 查询（含"新位点"子查询）；
- 专家月度人工评审新文献；既有位点至少季度复查；
- **社区贡献**走 GitHub（discussions + PR），数据文件 `STRchive-loci.json` + `check-loci.py` 校验 + Snakemake workflow（含反向互补 motif 自动填充）——版本控制避免"弃坑软件（abandonware）"。

### 3. 频率与致病性推断
- 用 **gnomAD v3.1.3 的 TR 基因型（~18,500 名非受累个体）**估计致病基因型频率——一个位点的"异常频率"有了人群基线；
- 交叉评审结论示例：作者论证 DMD 的 STR 位点虽在多个目录中列出、但很可能并非致病——目录本身的批判性维护。

### 4. 面向 genotyper 的输出
- `STRchive-disease-loci.hg38.longTR.bed`：**比 gnomAD/EH 校准定义更宽的位点窗口**（最大位点 ~4,000 bp），专为长读分型设计；
- STRdust `--pathogenic` 直接内置消费 STRchive；TRGT 官方文档推荐替换；LongTR README 有配合 STRchive 的示例。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 消费方式 | TRGT/LongTR/STRdust 目录 BED；pipeline 内 JSON 直读 |
| 更新 | `git pull` + 重跑分型目录检查；关注版本号 |
| 浏览 | strchive.org（含诊断蓝图 blueprint：等位基因评估/表型/位点评估流程） |

## 性能与局限

- 73 个位点——**只覆盖疾病位点**，全基因组扫描仍需 Adotto/UCSC/TRF 类目录；
- 人工维护为主，极端前沿位点的收录有滞后窗口；
- 等位基因大小范围是文献汇总值，并非连续更新的人群分布（分布层面参考 pathSTR）。

## 对 CCS-STR 优化的启示

- **把"目录版本"当配置管理**：TR 目录不是常量——升级 STRchive 版本会改变分型结果，应纳入回归测试的变量清单；
- 位点定义宽窄直接影响 spanning read 判定（宽定义 = 更容易捕获跨读）——HiFi 靶向 panel 调优时可对照 STRchive 与 gnomAD 两套定义的差异位点；
- 其"自动化文献监测 + 社区 PR"的模式是任何生物数据库可持续维护的参考架构。

## 参考

- 论文：https://doi.org/10.1186/s13073-025-01454-4 （PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC11938676）
- 网站：https://strchive.org ｜ 仓库：https://github.com/dashnowlab/STRchive
- 配套人群分布：pathSTR（Genome Research 2024）https://pmc.ncbi.nlm.nih.gov/articles/PMC11610575
