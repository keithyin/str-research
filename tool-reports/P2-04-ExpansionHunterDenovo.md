# ExpansionHunter Denovo（EHdn）

> **优先级**：P2（背景借鉴）｜ **类别**：免目录重复扩增**发现**（短读）｜ **维护方**：Illumina（Dolzhenko 等）｜ **状态**：维护模式（Apache 2.0）
> **一句话**：不看任何目录，用「完全在重复内的读（IRR）」的周期性检测 + 病例-对照/离群统计，在短读 WGS 里找已知与全新的重复扩增。

## 是什么

ExpansionHunter Denovo 是 ExpansionHunter 的"发现"姊妹工具：不按位点分型，而是产出每样本的**STR profile**（每 locus / 每 motif 的 IRR 计数），再做群体比较。

- 仓库：https://github.com/Illumina/ExpansionHunterDenovo
- 论文：*ExpansionHunter Denovo: a computational method for locating known and novel repeat expansions in short-read sequencing data*, **Genome Biology 2020**;21:102 https://pmc.ncbi.nlm.nih.gov/articles/PMC7187524

## 解决了什么问题

1. **发现没有目录的新扩增**：新致病位点（如 CANVAS 的非参考 AAGGG motif）在任何现有目录里都不存在，靶向工具永远看不到。
2. **大 cohort 没有长读**：百人–万人级短读队列是当时的现实，需要纯短读的重复离群检测。

## 如何解决的

### 1. IRR（in-repeat read）检测
- 候选：unmapped 或 MAPQ<40 的读中「完全位于重复内部」的读；
- 检测算法：对 motif 长度 k ∈ {2..20}（排除 homopolymer k=1）做周期性检验，要求 S(k) ≥ 0.8 且**加权纯度（WP）≥ 0.9**（低质量错配少罚）；
- 两类 IRR：
  - **Anchored IRR**：mate 以 MAPQ≥50 锚在重复旁 → 可定位到 ~locus（500 bp 内锚合并）；
  - **Paired IRR**：两 mate 都是同 motif 的 IRR → 极长（超片段）重复，只能按 motif 计数，无法定位。

### 2. 统计比较
- profile 按覆盖度归一（anchored IRR 归一到 40×）；
- **Case-control**：每 locus（或每 motif）单侧 Wilcoxon rank-sum；
- **离群分析**（异质队列）：对 95% 分位数分布做 bootstrap，超出者算 z-score 排序——与 exSTRa/STRetch 同思路。

### 3. 免目录的本质
不需要基因组坐标（anchored IRR 的 mate 即坐标）、能发现**组成变化**（motif 替换）、支持最长 20 bp motif——Rafehi 等指出它是**唯一能找到 CANVAS AAGGG 扩增**的计算方法；STRetch 漏掉的 CSTB 12-mer 它能检出。

## 性能

- Case-control：DMPK/FXN/HTT/FMR1 位点显著富集（14–35 病例 vs 150 对照）；
- 离群：单样本数据集中 81% 的扩增进入 top 10；模拟 41/44 致病重复检出（含 8/9 复杂非参考扩增）；
- 77% 的 ≥2 读支持重复在 PacBio 长读中得到验证；每人 ~1574 motifs、~660 个长重复位点；
- **局限：检不出短于读长的扩增**（没有 IRR 就没有信号）。

## 对 CCS-STR 优化的启示

- **IRR 思想直接可移植到 HiFi**：HiFi 读中"完全在重复内的读"同样存在且信息更丰富（可读出 motif 序列）——TRGT 的 in-repeat 读输出已经内置了原料，缺的是 EHdn 式统计层；
- 免目录发现与目录分型是**互补的两条腿**：你的管线若只跑目录分型，等于放弃了 novel 发现能力；
- WP（加权纯度）周期性检验是廉价高效的 motif 探测器，比上 TRF 轻。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC7187524
- 仓库：https://github.com/Illumina/ExpansionHunterDenovo
