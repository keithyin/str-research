# TRGT（Tandem Repeat Genotyping Tool）

> **优先级**：P0（核心必读）｜ **类别**：TR/STR 基因分型 ｜ **维护方**：PacBio ｜ **状态**：活跃，v5.1.0（2026-11）
> **一句话**：PacBio 官方的 HiFi 串联重复分型工具，一次输出**等位基因全长序列 + 长度 + motif 组成/中断 + CpG 甲基化 + 定相证据**，是当前 HiFi STR 分型的事实标准。

## 是什么

TRGT 是 PacBio（Dolzhenko / English / Dashnow，Eberle 实验室）开发的靶向（目录驱动）串联重复分型工具，Rust 编写，随 HiFi-human-WGS-WDL 与 SMRT Link（内部版本号独立）发布。配套可视化工具 **TRVZ**（`trgt trvz`）和队列数据库 **tdb**（原 TRGTdb）。

- 仓库：https://github.com/PacificBiosciences/trgt （~146★，v5.1.0）
- 论文：*Characterization and visualization of tandem repeats at genome scale*, **Nature Biotechnology 2024**（Dolzhenko, English, Dashnow 等）https://www.nature.com/articles/s41587-023-02057-3
- 许可证：**PacBio 专有许可（非开源）**，且 §3 明确「只能用于 PacBio 仪器产生的或 PacBio 提供的数据」——这是 TRGT 不能跑 ONT 的原因，也是 STRkit/LongTR 等开放工具存在的部分理由。

## 解决了什么问题

1. **长度之外的致病信息无法获得**。重复扩增疾病的致病性常取决于 motif 组成与**中断（interruption）模式**（如 FMR1 的 AGG 中断、RFC1 的 AAAAG/AAAGG/AAGGG 混合），而不是纯长度。此前长读工具（Straglr、tandem-genotypes 等）基本只给长度或拷贝数。
2. **短读工具在超读长扩增上失效**：ExpansionHunter 一系虽然用图模型突破读长限制，但短读无法解析重复内部序列，GC-rich 位点（如 FMR1 CGG）覆盖也差。
3. **同一份分析里拿不到甲基化**。HiFi 天然携带 5mC 动力学信息，但没有工具把「重复等位基因 × CpG 甲基化」绑在一起输出（FMR1 全突变与甲基化状态强相关）。
4. **群体规模分析没有合适的表示**：TR VCF 是"替换式"表示（REF/ALT 是整段区域），多样本合并与查询困难。

## 如何解决的

TRGT 的方法是一条「**跨读 → 聚类 → 每等位基因共识 → 结构注释 → 统计分型**」的流水线（以下均出自论文 Methods）：

### 1. 读段选择：靠侧翼定位跨读
在每个目录位点上，TRGT 找出**在重复两侧翼（flanks）都有比对锚**的 HiFi 读段作为 spanning reads。这决定了它的召回上限：等位基因必须被单条读跨越（论文给出的实用上限 ~10 kb；mosaic 致病扩增建议 ≥20× 覆盖度）。

### 2. 等位基因聚类：编辑距离 + Ward 层次聚类
- 对所有 spanning reads 计算两两**编辑距离**；
- 用 **agglomerative clustering（Ward linkage）** 聚类；
- 丢弃支持 <10% 的簇；
- 结合 karyotype（可用 `--karyotype` 提供性别核型，基准中用 Somalier 推断）：二倍体位点取最大两簇为两个等位基因，单倍体（chrX/chrY）取最大一簇。

**设计要点**：HiFi 单读准确率 ~99.9%，同源等位基因间的编辑距离分布是双峰的，这一步本质上是「用聚类替代比对歧义消解」——不需要把读重新比对回参考。

### 3. 每等位基因共识：中位长度骨架 + 多数投票
- 从每簇挑一条**中位长度**的读作比对骨架；
- 其余读比对到骨架；
- 沿骨架扫描，对每个位置的比对操作（match/insertion/deletion）做**多数投票**，得出该等位基因的共识序列。

这个共识序列会作为 VCF 的 ALT 输出——**下游一切（组成注释、可视化、体细胞分析、de novo 检测）都建立在这条序列上**，是整个 TRGT 设计的核心复用点。

### 4. Motif 组成与中断注释
- **简单重复**：用 `(m1)n1 s1 (m2)n2 …` 的 motif 表达式描述。注释算法把「motif 匹配节点 / 中断序列节点」建成 DAG，节点连向下一个拓扑出现位置，**取覆盖最多碱基的最长路径**作为 allele 的 motif 分段。
- **复杂重复**：用 **HMM**（沿用 adVNTR 风格）建模，TRGT 可为任意 motif 集合自动合成 HMM（如 RFC1 的 AAAAG/AAAGG/AAGGG），也支持自定义拓扑。
- 输出每个等位基因的 motif 计数（如 FMR1 的精确 CGG 拷贝数），中断以 mismatch 形式在 TRVZ 图上可见。

### 5. 统计分型与定相：侧翼 SNP 模型
对同源重复、mosaic 等聚类失效的情形，TRGT 回退到**侧翼 SNP 基因分型模型**：
- 每条读编码为「侧翼杂合 SNP 上是否携带 alt」的二元向量；
- 用标准 SNP 分型模型 `P(G|R) ∝ P(R|G)·P(G)`，其中 `P(R|G) = Π_r [P(r|H₁)+P(r|H₂)]/2`，先验 P(G) 可从对照队列估计；
- 选出最可能的基因型并把读分配到单倍型——**等位基因特异性甲基化**即由此按分配后的读集计算。

### 6. CpG 甲基化
对每个等位基因：取支持该等位基因的所有读，对等位基因区间内全部 CpG 的甲基化水平求均值（FORMAT 字段 `AM`）。这是 HiFi kinetics/base-modification 的直接红利。

### 7. 生态位
- 输出 **spanning BAM**（每位点跨读的读段切片，社区称 "BAMlets"）——trgt-instability 与 TRGT-denovo 都消费它；
- `merge` 子命令合并多样本 VCF；队列分析推荐转 **tdb**（Parquet 三表模型，单样本 ~11.4 Mb vs VCF 92.3 Mb，小 87.6%）；
- **TRVZ** 可视化：顶部 track 是「共识 vs 理想纯重复」的对齐（实线=match、灰=mismatch、竖线=insertion），下部是读 pileup + 可选甲基化。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | 比对好的 HiFi BAM（**pbmm2 `--preset HIFI`**，与 WGS 相同即可）+ 参考基因组 + 重复目录 BED |
| 目录来源 | **STRchive**（官方推荐，自带 pathogenic BED 已弃用）、Adotto（全基因组 ~937k 位点）、Broad catalog |
| 主要输出 | VCF（等位基因全长序列 + 长度 + SD/ALLR/AM 等字段）、spanning BAM、`trvz` 图（svg/pdf/png） |
| 常用参数 | `--karyotype`（性别核型）、`--flank-len`（侧翼上下文，基准用 25）、`--max-depth`（位点最大读深，基准 10000） |
| 官方工作流 | HiFi-human-WGS-WDL 内置 TRGT 模块；SMRT Link v25.2 亦有 GUI 入口 |

## 性能与局限

**性能**（论文自报，HG002）：
- vs HG002 T2T 组装：**98.06%** 等位基因完全一致或 ±1 bp；
- 孟德尔一致率（HG002/3/4 trio）：**89.00% 精确 / 98.38% ±1bp**（86.57% 的错误是 off-by-one）；
- 同代对比（长度口径 Mendelian 精确/±1bp）：tandem-genotypes 86.44%/96.08%，Straglr 70.05%/90.34%，短读 GangSTR 63.89%/72.90%；
- 速度：30× 基因组、32 核 35 分钟；
- 临床验证：FMR1（~1 kb 全突变）、DMPK（>5 kb）、STARD7、ATXN10（>4 kb）全部近碱基级正确。

**局限**：
- **必须有跨读**——超过读长（~10 kb 级）或覆盖不足的扩增检不出；STRkit 论文指出「无 spanning read」时 TRGT 直接失败，这正是 LongTR/otter 用非跨读证据补位的场景；
- 聚类在同源重复（homopolymer/双核苷酸）上孟德尔精确率只有 82.95%，off-by-one 常见；
- 专有许可仅限 PacBio 数据；
- STRkit 基准中 TRGT 虽 F1 略低于 STRkit（0.9574 vs 0.9631），但**速度最快**。

## 对 CCS-STR 优化的启示

- **一切下游都挂在 per-allele consensus 上**：若你要优化「体细胞/嵌合检测」，直接复用 spanning BAM + consensus 算 divergence（trgt-instability 就是这么做的）；
- 聚类簇阈值（10%）、`--flank-len`、`--max-depth` 是显式调参入口——改动 ccs/比对层后，用「TRGT 簇大小分布」可以灵敏地检测上游读质量变化；
- off-by-one 错误集中在 homopolymer/双核苷酸 → 上游错误模型（DeepConsensus 对 homopolymer 窗口）是可预期的改进方向；
- 目录质量直接决定分型质量（STRchive 的位点定义比 gnomAD/EH 校准的更宽，专为长读优化）。

## 参考

- 论文：https://www.nature.com/articles/s41587-023-02057-3 （PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC11921810）
- 仓库/文档（含 TRVZ 图解读、PureTarget 指南）：https://github.com/PacificBiosciences/trgt
- 配套：tdb（P1 报告）、trgt-instability（P1 报告）、TRGT-denovo（P2 报告）
- 许可证全文：仓库 `LICENSE.md`（PacBio Software License Agreement）
