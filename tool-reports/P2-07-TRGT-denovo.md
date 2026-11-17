# TRGT-denovo

> **优先级**：P2（背景借鉴，PacBio 生态新工具）｜ **类别**：家系 de novo TR 突变检测 ｜ **维护方**：PacBio（Mokveld / Dolzhenko / Eberle 等）｜ **状态**：bioRxiv 2024（Rust，repo: PacificBiosciences/trgt-denovo）
> **一句话**：在 TRGT 之上做 trio 分析：把亲子三方的读都对到**子代等位基因共识**上，用 parental-upper-bound 统计把 de novo TR 候选从 ~8 万压到 137 个（95% 验证率）。

## 是什么

TRGT-denovo 消费 TRGT 的输出（VCF + per-locus 读段切片 BAMlets），对每个家系位点独立判定是否存在 de novo 重复突变（长度**或**组成变化）。Rust 实现。

- 仓库：https://github.com/PacificBiosciences/trgt-denovo
- 论文（preprint）：*TRGT-denovo: accurate detection of de novo tandem repeat mutations*, bioRxiv 2024-07-19 https://pmc.ncbi.nlm.nih.gov/articles/PMC11275785

## 解决了什么问题

1. **纯基因型计数无法发现 de novo**：HG002/3/4 trio 在 891k 位点目录上，严格孟德尔不一致位点高达 **81,308 个（9.4%）**——绝大多数是分型噪声（homopolymer/双核苷酸 off-by-one）。等位基因层面 ±1 motif 也还有 2,582 个。
2. **组成级 de novo 完全不可见**：motif 内的 SNV/小 indel 不改变长度，基因型比较漏掉（81,308 vs 2,582 vs 137 三个数量级的差距正来自此）。
3. **候选确认成本高**：罕见病 trio 里 de novo TR 候选动辄数千，实验验证无法穷举。

## 如何解决的

### 1. 六等位基因共识 + 读重分配
- 每个等位基因 consensus + 50 bp 侧翼 → trio 共 6 条共识序列（子代 C0/C1、父 F0/F1、母 M0/M1）；
- 读按三种策略分配到等位基因：**TRGT 原生分配**（默认）／**比对法**（Wavefront Alignment（WFA）端到端比对打分、gap-affine 两段式评分，同分时按定相或 TRGT 分配破平）／**定相法**（HiPhase/WhatsHap phase tags + 侧翼 mismatch）。

### 2. 全员对齐到子代共识 + PU 统计（核心）
- **所有**三方读都对到子代的每条共识序列上（WFA 分数：0=完美，越负越不像）；
- 若子代等位基因是遗传来的：亲代读的分数分布应与子代自身读的自比对分布**重叠**；
- 若是 de novo：子代分布整体右移（远离两亲代分布）；
- **PU（parental upper bound）** = 亲代读分数分布的 **1.0 分位数**；子代读中分数高于 PU 的即"de novo 覆盖"（例：19 条读中 18 条超 PU → de novo ratio ≈0.95）；
- 辅助统计：CL（子代中位数下界，衡量亲代读重叠覆盖 = 遗传证据）、Magnitude（de novo 读与 PU 的平均差）。

### 3. 过滤准则（起始建议值）
- de novo 覆盖 ≥5 条读；等位基因 de novo ratio ≥0.7；
- 子代 de novo ratio ∈ [0.3, 0.7]（防 dropout/mosaic/样本调包）；
- dropout 概率低（每单倍型 ≥1 读，用样本内定相）；
- 官方强调：阈值应随覆盖度与研究设计调整，并对照全基因组分布。

## 性能

- HG002/3/4（30× Revio，891k 位点目录）：81,308 → **4,214 个有证据位点 → 137 个候选**；与纯基因型法重叠仅 0.12%；
- 8 个罕见病 trio 回测：旧目录下 18/18 已测序假阳性全部改判真阴性；937k 目录下每 trio 候选 ~100 个（vs 基因型法的数千–数万）；
- 实验验证：随机 20 个候选长程 PCR + Sequel IIe 测序——**19/20 确认**（唯一失败是双核苷酸收缩，事件幅度小于高深度下的读间波动）；
- 应用案例：AFF3 GCC 重复 trio 分析中提名了全基因组第一的 DNM。

## 局限

- 依赖 trio 与足够覆盖（每成员每等位基因 ≥5 读的门槛）；
- 双核苷酸类小事件受读间波动限制；
- preprint 阶段，阈值语义可能继续演进。

## 对 CCS-STR 优化的启示

- **「对齐到共识序列而不是参考」**在家族内是完全自洽的参照系——这个 trick 让长度/组成突变统一进一个比对分统计，值得在亲子验证、样本一致性检验中复用；
- PU 分位数统计 + 三道 ratio 过滤，是把"分型噪声"和"真突变"分离的可复制模板；
- 反过来看：**81,308 个孟德尔不一致位点里有 95 个是真 de novo**——分型错误率的上限直接决定 de novo 检测的灵敏度，这正是优化 CCS-STR 层的价值所在。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC11275785
- 仓库：https://github.com/PacificBiosciences/trgt-denovo
- 前置：TRGT（P0 报告）
