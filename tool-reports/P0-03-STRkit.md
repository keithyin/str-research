# STRkit

> **优先级**：P0（核心必读）｜ **类别**：STR 读级分型 + 定相 ｜ **维护方**：McGill（Lougheed / Pastinen / Bourque）｜ **状态**：活跃，≥v0.25.0（2026-11）
> **一句话**：用**近邻 SNV 把读实时分到单倍型**的读级 STR 分型器，主打「无需 pedigree 的定相」与「读级输出」，GPL-3.0 不限平台（HiFi/ONT 都能跑）。

## 是什么

STRkit 是 McGill 基因组中心开发的 STR 分型与分析工具箱：Python 前端 + Rust 编译扩展（strkit_rust_ext，parasail-SIMD 向量化比对），PyPI/Docker 分发，带 Web 可视化与孟德尔遗传模块。

- 仓库：https://github.com/davidlougheed/strkit （~22★，1476 commits，持续活跃）
- 论文：*Read-level genotyping of short tandem repeats using long reads and single-nucleotide variation with STRkit*, **Genome Research 2026**;36(3):578–588（preprint bioRxiv 2025-03）https://pubmed.ncbi.nlm.nih.gov/41539721
- 许可证：**GPL-3.0，无平台限制**——与 TRGT 的 PacBio 专有许可形成鲜明对照（STRkit 论文明确指出「TRGT 因许可不能在 ONT 上运行」）。

## 解决了什么问题

1. **定相依赖外部信息**：TRGT 聚类分型在同源重复上失效，用侧翼 SNP 模型兜底但需要 karyotype/队列先验；LongTR 依赖外部 haplotag。STRkit 的答案是**现场 call SNV、用 SNV 分组读**——不需要 pedigree、不需要预先 haplotag 的 BAM。
2. **长扩增被读长偏置低估**：覆盖分布天然偏向短等位基因（短等位基因更容易被完整覆盖）。简单取读长均值会系统性低估长扩增。
3. **只给一个 GT 数字不够**：eSTR 关联、体细胞分析需要**读级**拷贝数与分配结果，以及可信的置信区间。
4. **TRGT 许可封闭、ONT 无官方工具**：需要一个开放许可、跨平台的选项。

## 如何解决的

### 1. 候选 tract 生成 + 向量化半全局比对（核心打分引擎）
- 每条读按其 STR 片段大小生成一系列**不同拷贝数的候选重复 tract**；
- 用 **parasail（SIMD）半全局比对**把候选 tract 对到读上：5′/3′ 各留侧翼锚定、允许整段 tract 删除；得分最高的候选即该读的拷贝数；
- Rust 扩展把这一步向量化——这是它比 Straglr/STRdust（数千核-分钟级）快一个量级的原因。

### 2. 读重加权：校正「读长 × 等位基因长度」偏置
每条读的权重 = **该读覆盖此长度 tract 的估计概率的倒数**（由读长与位点覆盖度推出）。直觉：一条 20 kb 的读覆盖 500 bp 等位基因很容易、覆盖 15 kb 等位基因则很罕见，后者应占更高权重。加权后做 **~100 次重采样（bootstrap）**，顺带产出 95%/99% 置信区间。

### 3. 等位基因判定：三档 SNV 策略 + GMM 兜底
| SNV 情形 | 等位基因判定方式 |
|---|---|
| 信息 SNV 充足 | 用 SNV 差异矩阵 + 层次聚类把读分成两组（等位基因） |
| SNV 中等 | 复合距离：SNV 差异权重 5:1 或 10:1 vs 拷贝数差异 |
| 无/极少 SNV | bootstrap **GMM**（默认双峰；一峰权重可忽略时塌缩为纯合模型） |

### 4. 现场 SNV calling 与局部定相（招牌机制）
- 从读本身**现场 call 杂合 SNV**（过滤：读支持、碱基质量、能把读分成两组的能力），不是把 VCF 当真值；
- SNV 直接来自比对读，天然与 STR 位点**同读共现**——HiFi 15–20 kb 读长使「SNV+STR 同读」的概率远高于短读，这是该方法成立的生物学前提；
- 共享相同 SNV 单倍型分组的位点之间做**局部定相**；SNV 还可导出为 VCF；
- 自评 SNV 质量：对 GIAB 小变异基准，错误率 0.008%–0.126%；phase-set 错误率 0.31%–0.51%。

### 5. 输出与周边
- VCF/TSV/JSON（JSON 可在内置 Web 工具里看读数直方图与 motif k-mer 分布）；
- 每读拷贝数、峰分配、motif k-mer 量化、等位基因共识序列；
- `strkit mi`：孟德尔遗传分析 TSV；靶向测序模式（重采样调整）；`--use-hp` 消费已有 haplotag；可选 `--realign` 挽救 soft-clipped 的扩增读。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | BAM/CRAM（minimap2/pbmm2 比对均可）、参考 FASTA、STR 目录 BED（TRF/TRGT/bed4 格式，支持 IUPAC 简并 motif）、可选 dbSNP VCF |
| 目录 | 直接消费 TRGT 格式目录（Adotto/STRchive 可通用） |
| 输出 | VCF + JSON/TSV；读级拷贝数与 SNV；CI；可视化网页 |
| 注意 | **无甲基化输出**（甲基化仍是 TRGT 独占）；不能分型超过单读长度的等位基因；无 de novo 发现 |

## 性能与局限

**性能**（论文基准，HG002，GIAB TR v1.0 的 914,676 位点子集，Truvari/Laytr 评测）：
- HiFi ~32×：**F1 0.9631**（precision 0.9717 / recall 0.9545），其基准中最优；TRGT recall 更高（0.9675）但 F1 0.9574；LongTR 0.9555；STRdust 0.8327；
- ONT R10 simplex ~32×：**F1 0.9544 最优**（去 SNV 模式 0.9440；LongTR 0.9354；STRdust 0.7581）；
- ONT duplex ~12×：0.9051 最优（LongTR 0.8786）；
- 孟德尔一致率：Ashkenazi trio 序列级 **97.56%（HiFi）/ 94.44%（ONT）**——三者中最好；
- 速度：数百核-分钟（TRGT/LongTR 同量级；Straglr/STRdust 需数千）。

**局限**（论文自述）：
- 超读长等位基因不能分型（目录驱动、无 de novo）；
- Δ 等位基因大小 3–500 bp 区间可能遇到非规范 motif/嵌入非重复序列的建模低谷；
- SNV 无法完全解释 ONT 上的性能增益（机理未完全清楚）。

## 对 CCS-STR 优化的启示

- **读重加权是你最值得抄的机制**：任何「按读长聚合估计等位基因」的模块（包括你自己写的任何共识/统计）都存在同样的长度偏置；
- 现场 SNV 定相说明 **HiFi 读的 SNV-STR 联合信息还没有被 TRGT 完全利用**——一个优化方向是给你的分型/嵌合模块加 SNV 锚定；
- 它的 bootstrap-CI 框架可以直接用于量化「上游改动（ccs 参数/比对器）对 STR 分型稳定性」的影响；
- GPL-3.0 + 开放平台：可以放心地魔改源码做实验（对比 TRGT 只能当黑盒调参）。

## 参考

- 论文：https://pubmed.ncbi.nlm.nih.gov/41539721 （PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC12951957；preprint: https://www.biorxiv.org/content/10.1101/2025.03.25.645269）
- 仓库：https://github.com/davidlougheed/strkit
- Zenodo 存档：10.5281/zenodo.12689906
