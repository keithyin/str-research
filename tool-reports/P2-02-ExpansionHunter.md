# ExpansionHunter

> **优先级**：P2（背景借鉴）｜ **类别**：短读 STR/重复扩增分型（序列图）｜ **维护方**：Illumina（Dolzhenko 等，后赴 PacBio）｜ **状态**：活跃（v5 线，~220★）
> **一句话**：用**带环节点的序列图**重比对读段，第一个能稳定分型超读长扩增的短读工具；其 variant catalog 格式成为整个领域的事实接口。

## 是什么

ExpansionHunter 从 BAM/CRAM 中对目录化 STR 位点估计重复大小，C++，Apache 2.0（新版 Illumina 发行采用 PolyForm Strict 包装）。作者 Dolzhenko 后来去了 PacBio 做 TRGT——一脉相承。

- 仓库：https://github.com/Illumina/ExpansionHunter
- 论文：*ExpansionHunter: a sequence-graph-based tool to analyze variation in short tandem repeat regions*, **Bioinformatics 2019**（ISMB）https://pmc.ncbi.nlm.nih.gov/articles/PMC6853681
- 周边：REViewer（可视化）、ExpansionHunter Denovo（发现）、DRAGEN-STR（产品化）。

## 解决了什么问题

1. **短读无法跨越扩增**：读长 100–150 bp，致病扩增动辄数百–数千 bp——此前工具只能报"扩展存在"或依赖 read-pair 统计。
2. **重复区比对歧义**：普通比对器把重复读 clip 或错位；需要"知道重复结构"的定向重比对。
3. **复杂位点结构**：真实致病位点常是复合结构（如 HTT 的 (CAG)\*CAACAG(CCG)\*），单一 motif 建模不够。

## 如何解决的

### 1. 变异目录 → 序列图
- 位点用**受限正则式**描述（支持 IUPAC 简并符号）：`(CAG)*CAACAG(CCG)*`（HTT）、`(A)*(GAA)*`（FXN）、`(GCN)*`（聚丙氨酸）等；
- 每条正则编译成**序列图**：节点存序列、有向边定义可行连接；**重复单元用单节点自环（loop）表示**——一个节点即可代表任意拷贝数（GraphTools C++ 库的关键特性，此前的图方法多为全无环图）。

### 2. 读分类 → 图上重比对 → 路径拼接
- 每位点抽取相关读（spanning / flanking / 完全在重复内的 in-repeat 读）；
- 读重比对到位点图；**对齐路径拼接（stitching）**后进入分型；
- 重复按 2017 方法（EH 前身）分型，SNV/indel 用 Poisson 模型；
- **超读长扩增**：spanning 读消失时，靠 in-repeat 读 + flanking 读的图证据组合分型——首次让短读 WGS 在超读长扩增上达到 precision 0.91 / recall 0.99（GangSTR 0.88/0.83，TREDPARSE 0.84/0.46）。

### 3. 复合位点的联合图
HTT 模拟数据显示：把 CAG 与 CCG 画进**同一张图**联合分析，精度显著高于分开分析——"位点结构要整体建模"的实证。

## 输入输出与使用要点

- 输入：BAM/CRAM + variant catalog（JSON，格式文档完善，被 LongTR/STRkit 等借用作目录格式参考）+ 参考；
- 输出：VCF；配合 REViewer 时额外产出 graph-realigned BAMlet；
- 注意：**短读专用**，不建议直接跑长读 BAM（社区多次踩坑：EH 忽略 unmapped 读、不支持长读 CIGAR 语义）。

## 局限

- 短读平台限制；GC-rich 位点（FMR1）覆盖偏差大；
- 长扩增的 in-repeat 读证据统计精度依赖 fragment size 分布假设；
- 已被 DRAGEN-STR 产品化，开源版功能冻结于核心分型。

## 对 CCS-STR 优化的启示

- **「自环节点表示任意拷贝数」的图压缩**是表示重复的优雅方案——若你做多 motif 位点的重比对，图模型比字符串拼接更干净；
- 正则式目录（含 IUPAC）作为"位点 DSL"被广泛沿用——设计目录格式时可兼容它；
- 「复合位点必须整体建模」教训：拆开分析 CAG/CCG 会丢精度，TRGT 的 motif 表达式/HMM 是其 HiFi 版答案。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC6853681
- 仓库：https://github.com/Illumina/ExpansionHunter
- 产品化：https://www.illumina.com/science/genomics-research/articles/str-expansionhunter.html
