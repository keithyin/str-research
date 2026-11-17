# HipSTR

> **优先级**：P2（背景借鉴）｜ **类别**：短读 STR 分型/定相 ｜ **维护方**：Tom Willems / Erlich-Gymrek 一系｜ **状态**：维护模式（最后 release 2020-04）
> **一句话**：短读时代的经典：**EM 学习 per-locus stutter 模型 + HMM 重比对 + SNP 单倍型定相**——LongTR 的直系祖先，其「错误模型可学习」思想仍是 HiFi 工具的参照系。

## 是什么

HipSTR（Haplotype-based STR genotyper）从 Illumina WGS 数据对 STR 做基因分型与定相，C++ 实现。

- 仓库：https://github.com/tfwillems/HipSTR （最后 release 2020-04）
- 论文：*Genome-wide profiling of heritable and de novo STR variations*, **Nature Methods 2017**;14:590–592 https://pmc.ncbi.nlm.nih.gov/articles/PMC5482724
- 性能：vs 毛细管电泳（Marshfield panel）**95.2%**，丢弃最不自信 10% 后 **98.9%**（达到 CE 平台一致性上限）；lobSTR 88.2%、RepeatSeq 57.8%。

## 解决了什么问题

1. **PCR stutter 污染**：PCR 扩增会让读的 STR 长度系统性偏移 ±1..n 个单元——不建模 stutter，所有长度信号都是脏的。
2. **单倍型背景下的分型**：群体数据里等位基因频率先验能显著改善分型。
3. **STR 定相**：STR 变异当时完全无法进入定相/单倍型分析体系。

## 如何解决的

### 1. per-locus stutter 模型 + EM 学习（对后世影响最大）
每位点 θ = {u, d, ρs}：u = stutter 加一个单元的概率、d = 减一个单元的概率、ρs = 几何分布步长参数。观测 stutter δ 单元的概率：δ=0 → 1−u−d；δ>0 → u·ρs(1−ρs)^(δ−1)；δ<0 → d·ρs(1−ρs)^(−δ−1)。
- **E 步**：在群体样本上算基因型后验与每读归属条件概率；
- **M 步**：u/d = 读 STR 尺寸大于/小于底层等位基因的加权比例；ρs = 非零 stutter 读的加权平均步长的倒数；等位基因频率同步更新。

### 2. 双区段比对模型（flank 与 STR 分开建模）
- **flank 区**：Dindel 式 HMM（Match/Ins/Del 三矩阵 + 转移概率）；
- **STR 区**：假设相对单倍型**至多一处 indel、幅度为单元长度整数倍**，大小按 stutter 模型加权；删除按所有可能位置求和，插入碱基按局部重复周期性打分；
- 两区段拼接时要求读在 STR 前后各有一个完全匹配的侧翼碱基——这就是**侧翼锚定**的短读实现；
- 候选等位基因：≥2 条读且 ≥20% 样本读支持；genotyping 迭代中还会把"带 stutter 的比对里反复出现的序列"补回候选集。

### 3. SNP 单倍型定相
- 每个落在已定相杂合 SNP 上的碱基贡献 P = base quality（匹配）或 (1−q)/3（错配）；STR 15 bp 内的 SNP 排除；
- phasing likelihood × alignment likelihood 联合打分基因型；纯合基因型给 2× 先验。

### 4. 过滤语义（后来者的模板）
DFLANKINDEL（侧翼 indel 读占比）、DSTUTTER（stutter 读占比）、MALLREADS（跨读数）——LongTR/TRGT 的输出字段直接继承了这套命名与思想。

## 局限与在 HiFi 上的失效

- 「单次整数倍 indel」假设在长读上崩溃（LongTR 论文的动机）；
- PCR stutter 模型对 HiFi 无 PCR 的场景无意义，但其"可学习错误模型"思想被 LongTR 转化为 homopolymer 错误模型；
- 超读长扩增不可分型（HipSTR 对 >读长等位基因一律报纯合参考）。

## 对 CCS-STR 优化的启示

- **错误模型应该 per-locus/per-context 可学习**（EM 或监督），而不是全局常数——HiFi 的对应物是「per-homopolymer-length / per-motif-context」的错误参数（LongTR 已实践，可进一步细化）；
- 双区段（flank+STR）分开建模、拼接锚定的结构，比单一全局打分更能容忍局部噪声；
- MALLREADS/DFLANKINDEL 式的读级 QC 字段是建立回归测试的好抓手。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC5482724
- 仓库：https://github.com/tfwillems/HipSTR
- WebSTR（hg19 STR 参考区下载）：https://webstr.ucsd.edu/downloads
