# STRling

> **优先级**：P2（背景借鉴）｜ **类别**：k-mer 免目录扩增发现（短读）｜ **维护方**：Quinlan lab（Broad 合作，Harriet Dashnow 一作——后为 TRGT/STRchive 作者）｜ **状态**：开源（MIT），Nextflow/Bpipe/WDL 工作流齐备
> **一句话**：专门回收**比对器丢弃的读**（soft-clipped / unmapped / 错位到"重复黑洞"）做 k-mer 分析，把其他工具看不到的扩增（含全新 motif）找出来。

## 是什么

STRling（STR + ling 从"lingering reads"意）从短读 BAM 中检测已知与新颖 STR 扩增位点。

- 仓库：https://github.com/quinlan-lab/STRling
- 论文：*STRling: a k-mer counting approach that detects short tandem repeat expansions at known and novel loci*, **Genome Biology 2022**;23:257 https://pmc.ncbi.nlm.nih.gov/articles/PMC9753380

## 解决了什么问题

1. **"重复黑洞（sink）"效应**：重复富集读会比对到参考里最长的同源重复上去——错位；扩展读的重复部分常被 soft-clip；全新 motif（参考里根本没有）整读 unmapped。**这些读被所有比对依赖型工具抛弃**。
2. **位置精度**：EHdn 对扩增位点的定位误差平均 ~713 bp——临床解释需要 bp 级定位。

## 如何解决的

### 1. k-mer 特征（k=2–6，非重叠）
- 候选读 = unmapped / MAPQ<40 / 比对到大 STR sink 区 / 带 soft-clip 的读；
- 以步长 k 扫**非重叠 k-mer**（比 assembly 式重叠 k-mer 更适合首尾相接的串联重复）；
- 对每个 k-mer 计算**所有旋转取最小**（容忍插入在重复中打断相位），整数化提速；
- 代表 motif = 覆盖读比例最大的 k-mer（≥80% 信息阈值；homopolymer 记 1 bp 单元）。

### 2. 读分类与定位
- **Anchored pair**：STR 富集读 + 高质量非重复 mate → 按中位 fragment size 从 mate 位置外推；
- **Unplaced pair**：双 mate 都异常 → 无法定位；
- **位点聚类**：≥5 条 informative 读（至少 1 条 anchored）起簇，迭代重算中位中心（丢弃超出 98% fragment +100 bp 的读、吸纳新区间内读）至收敛；soft-clip 边缘精确定界；
- **联合 merge**：多样本先 merge 发现位点（边界更准），再逐样本 call 兜底单样本 novel。

### 3. 大小估计与离群检验
- 亚读长等位基因用 spanning 读；超片段用 anchored/unplaced 计数（计数与大小线性，经模拟 300 个 HTT 等位基因 0–1800 bp 校准）；
- 离群统计：log2((STR 计数+1)/局部深度) 对队列中位数的单侧 z-score + BH 校正。

## 性能

- 134 例 WGS（94 致病携带者）：单样本 87.2%、**联合 calling 98.9%** 检出；离群 83%（除 FXS 外 >96%）；
- **定位精度**：联合调用平均误差 6.14 bp（中位 1 bp）vs EHdn ~713 bp（差 30 倍以上）；
- 新位点案例：CANVAS（AAGGG）5/5、DBQD2 1/1；SCA6 型的小扩展（仅超参考 ~26 bp）漏检；
- 1000 Genomes 联合分析：1549 已知 + 650 新可变位点；~1.43 CPU 时/样本。

## 局限

- 大等位基因系统性低估；~150 bp "灰区"；高 GC 位点（FXS、GGGGGCC）困难；
- 表观 de novo 假阳性偏高（HG002 的 de novo 无 PacBio 支持）；homopolymer 假阳性占 94.6%——过滤后每基因组仅 ~69 个显著非 homopolymer 扩增。

## 对 CCS-STR 优化的启示

- **「比对器丢了什么」是一个独立的信号源**：HiFi 流水里 soft-clipped/unmapped 读同样值得单独检查（虽然 HiFi 上比例小得多）——把它做成监控指标即可；
- 「旋转不变 k-mer」是容忍 motif 内打断的最小实现，可用于读级 motif 分类器的特征工程；
- 联合 merge 精化边界的思路适用于任何 cohort 层位点发现。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC9753380
- 仓库：https://github.com/quinlan-lab/STRling
