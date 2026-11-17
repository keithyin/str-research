# REViewer

> **优先级**：P2（背景借鉴）｜ **类别**：STR 可视化（短读）｜ **维护方**：Illumina（Dolzhenko & Eberle）｜ **状态**：**2026-04 已归档**（read-only）
> **一句话**：把 ExpansionHunter 的图重比对读生成单倍型分辨的 pileup 图——TRVZ 的短读前身；看懂它的图就是看懂 TRVZ 的图。

## 是什么

REViewer 为指定 STR 位点渲染"读-单倍型"pileup 图（svg/png），揭示分型偏差、中断结构、mosaic 比例。C++。

- 仓库：https://github.com/Illumina/REViewer （~92★，2026-04 归档）
- 原理文章：https://www.illumina.com/science/genomics-research/articles/reviewer-alignments-short-reads-long-repeat.html （Dolzhenko & Eberle, 2020-12）
- preprint：*REViewer: Haplotype-resolved visualization of read alignments in and around tandem repeats*, bioRxiv 2021-10-20, DOI: 10.1101/2021.10.20.465046

## 解决了什么问题

1. **VCF 里的 GT 看不出对错**：分型长度对但组成错、pileup 不对称（过/低估）、mosaic 分布——都需要"人眼检查读证据"。
2. **graph-realigned BAM 不可读**：EH 的中间产物需要专用渲染。

## 如何解决的

- 消费 **ExpansionHunter（v3+）生成的 graph-realigned BAMlet** + variant catalog + EH VCF + 参考；
- 渲染：按单倍型分组读 pileup；读与图参考对齐，插入/删除/中断直接可见；配合 EH 的等位基因特异性路径展示单倍型分辨视图；
- 参数化输出（locus、前缀、svg/png）。

## 周边（实用）

- **FlipBook**（Broad，Ben Weisburd）：批量 REViewer 图的图像服务器——大规模人工质检用；
- **Review BAMs**（Andreas Halman）：让 REViewer 能吃普通 BAM（内部替你跑 EH）；
- HiFi 世界的对应物是 **TRVZ**（TRGT 自带）：顶部 track 为「等位基因共识 vs 理想纯重复」对齐，方法上正是 REViewer 图逻辑的长读版。

## 局限

- 依赖 EH 特定输出（README 警告：用 EH 仓库自带示例输出会崩溃）；短读数据为主；
- 归档后不再演进——新项目直接用 TRVZ（HiFi）或 mahi-mahi 生态（ONT）。

## 对 CCS-STR 优化的启示

- **可视化是质检闭环的一半**：论文级/临床级的 STR 结果都应配套 pileup 图审读——建议在你的流程中为高影响位点自动出 TRVZ 图；
- 「共识 vs 理想纯重复」的对齐显示（TRVZ 的顶部 track）是暴露 motif 中断最直观的形式——若你自建可视化，照抄这个设计；
- 归档 ≠ 无用：其单倍型分辨 pileup 的视觉语义已成为领域通用语言（STRipy、TRVZ、FlipBook 都沿用）。

## 参考

- 原理文章：https://www.illumina.com/science/genomics-research/articles/reviewer-alignments-short-reads-long-repeat.html
- 仓库：https://github.com/Illumina/REViewer
- preprint：https://doi.org/10.1101/2021.10.20.465046
