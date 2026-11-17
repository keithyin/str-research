# PacBio CCS/HiFi STR 分析优化 —— 业界方法调研

> 调研日期：2026-11-17
> 范围：论文（期刊/preprint）、GitHub 仓库、博客与官方文档
> 视角：面向「CCS/HiFi 数据 → STR 基因分型」这条链路上的优化方法
>
> 说明：所有工具与论文均经网络核实；文末汇总全部链接。个别广为流传但不存在的工具（如 ISLPR、STRful）在 §10 中专门说明，避免踩坑。
>
> **配套报告**：文中提到的每个工具的详细调研报告（是什么/解决什么问题/如何解决）见 [tool-reports/](tool-reports/README.md)，按 P0–P3 优先级组织。

---

## TL;DR

**业界现状一句话**：HiFi STR 基因分型的事实标准是 PacBio 官方的 **TRGT**，最强互补者是 Gymrek lab 的 **LongTR**；上游决定性因素是 **CCS/DeepConsensus 的读质量**与 **pbmm2 的比对行为**；下游靠 **STRchive 目录 + tdb 队列库 + GIAB TR 真值集**闭环。

核心结论（详细依据见正文）：

1. **算法套路已经收敛**：主流工具都在做同一组事——侧翼锚定（flanking anchoring）→ 按等位基因聚类读段 → 构建 per-allele consensus → 用 HMM/似然模型（带 homopolymer 感知错误模型）重比对基因分型 → 定相。差异在于各自的取舍（速度 vs. 长扩增召回 vs. 序列解析度）。
2. **TRGT 与 LongTR 互补而非替代**：TRGT 强在 motif 组成、中断（interruption）、CpG 甲基化、可视化和群体数据库；LongTR 强在 >500 bp 的长扩增、null allele（整段缺失）检测和序列解析单倍型。
3. **上游读生成层仍有红利**：DeepConsensus 把 HiFi 残余错误减少约 42%（Nature Biotechnology 2023），2026 年 AlphaEvolve 优化训练后 Q30 读产出从 47.9%→53.2%；HiFi 残余错误集中在 homopolymer/STR 区，这正是 STR 分型的误差来源。
4. **别只信一个 caller**：Nature Protocols 2025 的官方建议是对关联研究产出**多工具 consensus 基因型**；基准测试结论也高度依赖真值集和指标，横向比较要小心。

---

## 1. 背景：HiFi 上做 STR，难在哪

先定义问题。PacBio CCS（Circular Consensus Sequencing）产出 HiFi 读：单分子多次聚合酶 passes → 高精度共识（中位 QV≈Q30+，~99.9%）。对 STR/VNTR 而言，HiFi 有天然优势（读长足够跨越大多数重复、无 PCR stutter），但仍有四类难点：

| 难点 | 原因 | 业界应对（详见后文） |
|---|---|---|
| **残余错误集中在重复区** | HiFi 的剩余错误偏好 homopolymer indel 与低复杂度区（综述实证） | homopolymer 感知 HMM 错误模型（LongTR/Dindel 式）；DeepConsensus 上游纠错 |
| **比对歧义** | 读段起点/终点落在重复内部时，比对器 clip 或错位；pbmm2 特定参数可能让 spanning read 不比对 | 侧翼锚定；切掉重复的人工参考重比对（STRdust sensitive）；otter 局部组装 |
| **超读长扩增** | 扩增超过 read length（如 FMR1 全突变）时无 spanning read | in-repeat/flanking 读类别建模（GangSTR 一脉）；graph realign（ExpansionHunter）；完全跨读长的 ONT 侧走 assembly |
| **复杂 motif 与嵌合** | 纯长度无法刻画致病性（如 RFC1、ATXN10 的 interruptions）；体细胞 mosaic 需要读级证据 | motif 表达式/HMM 定义重复（TRGT）；读级不稳定统计（trgt-instability） |

体量参考： tandem repeats 约占 hg38 的 **8.1%**，携带个体 **~24.9%** 的变异（GIAB TR 基准论文，English et al. 2024）。

另外一个实用分类（来自 [PacBlog: Sequencing 101 – Tandem repeats](https://www.pacb.com/blog/sequencing-101-tandem-repeats)）：**STR（motif 1–6 bp）** 对单读错误最敏感，HiFi 也未必稳；**VNTR（motif ≥7 bp）** HiFi 优势最大。选工具前先看你的目标位点属于哪类。

---

## 2. 全景图：一条 HiFi STR 管线的分层

把业界方案放到管线上看，每一层都有独立的优化抓手：

```
Subreads (ZMW)
   │
   ├─[1] 读生成        ccs / pbccs (QV 阈值、min-passes、--by-strand)
   │                   DeepConsensus (gap-aware transformer, +42% 错误削减)
   ▼
HiFi BAM
   │
   ├─[2] 比对          pbmm2 (HIFI/CCS preset, 默认 -Y 软剪切、无 secondary)
   │                   minimap2 对照 (spanning read 丢失时的排查手段)
   ▼
HiFi 比对 BAM
   │
   ├─[3] STR 基因分型   TRGT ★标准 │ LongTR │ STRkit │ STRdust │ Straglr │ tandem-genotypes
   │                   (目录驱动: STRchive / Adotto / Broad catalog)
   ▼
VCF + spanning BAM (+ 甲基化/组成)
   │
   ├─[4] 队列与群体     tdb (TRGTdb) │ gnomAD 式 variation clusters │ Multi-caller consensus
   │
   ├─[5] 体细胞/嵌合    trgt-instability │ Trio de novo (TRGT-denovo) │ MonSTR
   │
   ├─[6] 可视化/报告    TRVZ (TRGT 自带) │ trviz │ REViewer (短读侧) │ VarSeq
   │
   └─[7] 基准回归       GIAB HG002 TR truth set │ truvari │ 组装对照 (T2T / otter)
```

后面各节按这个分层展开。

---

## 3. 工具盘点

### 3.1 PacBio 官方栈（你优化工作的"主场"）

| 工具 | 作用 | 版本（2026-11） | 备注 |
|---|---|---|---|
| [ccs](https://github.com/PacificBiosciences/ccs) | Subreads → HiFi 共识 | v8.x（文档在 [ccs.how](https://ccs.how)） | BSD-3；Revio/Sequel IIe 已上机运行；关键参数 `--min-passes`、`--min-rq`、`--by-strand` |
| [DeepConsensus](https://github.com/google/deepconsensus) | CCS 后 ML 纠错 | 1.2.0 | Google 出品，Python/TF；注意：**ccs 默认丢弃 Q20 以下读，DeepConsensus 救不回来**，需按其 quick-start 调 ccs 参数 |
| [pbmm2](https://github.com/PacificBiosciences/pbmm2) | 官方 minimap2 前端 | v26.x | `HIFI/CCS` preset 为默认；强制 `-Y` 软剪切、/=X CIGAR、默认无 secondary；≥1.13.0 内置 minimap2 v2.26 |
| [TRGT](https://github.com/PacificBiosciences/trgt) | TR 基因分型（长度+组成+甲基化+定相） | v5.1.0，~146★，Rust | ★事实标准；VCF + spanning BAM + `trvz` 可视化；论文为 Nature Biotechnology 2024 |
| [trgt-instability](https://github.com/PacificBiosciences/trgt-instability) | 体细胞不稳定量化（基于 TRGT 输出） | 2025 预发布 v0.5.0 | read-vs-allele-consensus divergence + Dirichlet-multinomial 基线模型 + 参数 bootstrap 检验 |
| [HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) | 官方 WGS 流程（含 TRGT 模块） | v4.0.0 | 单样本/家系 workflow；TRGT v5.0.0 已入 CI；旧 snakemake 版已归档 |
| SMRT Link Tandem Repeat 分析 | GUI 入口 | v25.2 | 内部跑 TRGT（SMRT Link 版本号独立，v25.2 内为 TRGT v1.5.1） |

**TRGT 输出要点**（对优化最有用）：`--in-repeat-reads` 类输出可拆 spanning/flanking/in-repeat；per-allele consensus 用于组成与中断注释；复杂重复可用 **motif 表达式或 HMM** 定义（如 RFC1）；`merge` 子命令合并多样本 VCF。

### 3.2 长读长学术工具（论文 + GitHub）

| 工具 | 机构/实验室 | 核心算法 | 亮点/性能（自报口径） | 论文 | 仓库 |
|---|---|---|---|---|---|
| **tandem-genotypes (tdg)** | Frith lab (东京大学) | LAST/last-split 概率比对 + mismap 过滤 + 侧翼双锚定 + 按重复单元整数倍读长度变化 | 免目录 de novo 排序（coding>UTR>intron）；SAMD12 5kb 扩增案例 | [Genome Biology 2019](https://pubmed.ncbi.nlm.nih.gov/30890163) | [mcfrith/tandem-genotypes](https://github.com/mcfrith/tandem-genotypes) (~52★, Python) |
| **TRiCoLOR** | Magi lab (佛罗伦萨) | CIGAR indel 模式检测 TR + 读级聚合定长 | 早期 error-prone 长读（CLR/ONT）方案 | [GigaScience 2020](https://pubmed.ncbi.nlm.nih.gov/33034633) | （见论文） |
| **Straglr** | BCGSC (Birol lab) | TRF+blastn 定 motif → GMM 按等位基因长度聚类 | 兼做 de novo 发现与靶向分型；ONT EPI2ME 采用 | [Genome Biology 2021](https://pubmed.ncbi.nlm.nih.gov/34389037) | [bcgsc/straglr](https://github.com/bcgsc/straglr) (~174★, Python, v1.3.0) |
| **TRGT** | PacBio (Dolzhenko/English/Dashnow, Eberle) | 编辑距离聚类 → per-allele majority consensus → motif 表达式/HMM 拆解组成 → CpG 甲基化 | HG002: Mendelian 98.38% (±1bp)；2,000+ 临床位点案例全检出 | [Nature Biotechnology 2024](https://www.nature.com/articles/s41587-023-02057-3) | [PacificBiosciences/trgt](https://github.com/PacificBiosciences/trgt) |
| **TRGT-denovo** | PacBio + 合作 | Trio：子代读段重比对到等位基因 consensus + parental-upper-bound 统计 | de novo 候选 ~81k→137 个，95% 实验验证率 | [bioRxiv 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11275785) | （TRGT 生态） |
| **LongTR** | Gymrek lab (UCSD) + NIST | 改造 HipSTR：编辑距离聚类 → POA(spoa) 单倍型 consensus → Dindel 式 HMM（几何分布 homopolymer 错误模型） | HG002 组装对照：序列一致率 84.7%；Mendelian 86% vs TRGT 79%（Ashkenazi trio）；>500bp 重复优势最大；可检 514 个 null allele | [Genome Biology 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11229021) | [gymrek-lab/longtr](https://github.com/gymrek-lab/longtr) (~41★, C++, v1.2) |
| **STRdust** | De Coster/Höijer/Ameur (Uppsala) | fast=CIGAR walk 直接取重复序列；sensitive=比对到"切除重复"的人工参考；Ward/DBSCAN 单倍型聚类 | 30x ONT 全基因组 2,000 位点 7 秒；`--pathogenic` 直连 STRchive | [Bioinformatics 2024](https://doi.org/10.1093/bioinformatics/btae268) | [wdecoster/STRdust](https://github.com/wdecoster/STRdust) (~28★, Rust) |
| **STRkit** | McGill (Lougheed/Pastinen/Bourque) | GMM 估重复数 + parasail-SIMD 向量化比对 + **近邻 SNV 锚定定相**（无需 pedigree）+ bootstrap 置信区间 | F1: 0.9631 (HiFi) / 0.9544 (ONT)，其基准中最优；孟德尔一致性高于对手 | [Genome Research 2026](https://pubmed.ncbi.nlm.nih.gov/41539721)（preprint 2025） | [davidlougheed/strkit](https://github.com/davidlougheed/strkit) (~22★, Python+Rust 扩展, GPL-3.0) |
| **TREAT + otter** | Holstege lab (Radboud) | TREAT=端到端分型/motif/病例-对照离群分析；**otter=靶向局部组装器**（自适应错误率与覆盖度） | 5–30x 多平台基准：<500bp PacBio 最佳，≥1000bp ONT 更好；0.06% TR 出现覆盖度丢失（富集 GC-rich/G-四联体），加非 HiFi 读可救（ABCA7 相关性 0.58→0.99） | [Genome Research 2024](https://genome.cshlp.org/content/34/11/1942) | [holstegelab/treat](https://github.com/holstegelab/treat) |
| **ATaRVa** | CCMB + Dashnow lab | 平台无关（HiFi/ONT）、主打速度 | Genome Biology 2026 基准中 91.7%、约快一个数量级 | [bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.05.13.653434v2)（未经同行评审） | [SowpatiLab/ATaRVa](https://github.com/SowpatiLab/ATaRVa) |

> ⚠️ 各工具性能数字来自**各自论文的基准**，真值集与指标不同，**不可直接横向比较**；中性对比见 §5。

### 3.3 短读时代经典（值得借鉴的设计）

| 工具 | 论文 | 借鉴价值 |
|---|---|---|
| **lobSTR** | [eLife 2012](https://elifesciences.org/articles/00313) | 开创"识别重复读段→抽重复区→ML 分型"范式 |
| **HipSTR** | [Nature Methods 2017](https://pubmed.ncbi.nlm.nih.gov/28436466) | **EM 学习 per-locus stutter 模型 + HMM 重比对 + 侧翼 SNP 定相**——LongTR 的直系祖先；95.2% vs 毛细管电泳 |
| **ExpansionHunter** | [Bioinformatics 2019](https://pubmed.ncbi.nlm.nih.gov/31134279) | **序列图（loop node）realign**：首个稳定分型超读长扩增的短读方法；其 variant catalog 格式成为事实接口 |
| **GangSTR** | [NAR 2019](https://pubmed.ncbi.nlm.nih.gov/31194863) | **enclosing/spanning/FRR(完全重复)/flanking 四类读的联合最大似然**——后来所有长读工具分类读段的源头 |
| **ExpansionHunter Denovo** | [Genome Biology 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02017-z) | 免目录 IRR（in-repeat read）周期检测 + 病例-对照/离群统计 |
| **STRling** | [Genome Biology 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753380) | k-mer 计数回收 soft-clipped/未比对读段中的扩增证据——对"比对器丢弃了什么"这一问题的通用提醒 |

> 仓库：[tfwillems/HipSTR](https://github.com/tfwillems/HipSTR)、[Illumina/ExpansionHunter](https://github.com/Illumina/ExpansionHunter) (~220★, C++, PolyForm 许可)、[gymreklab/GangSTR](https://github.com/gymreklab/gangstr)、[Illumina/REViewer](https://github.com/Illumina/REViewer)（可视化，2026-04 已归档）。Gymrek lab 还有整套下游生态：[TRTools](https://trtools.readthedocs.io)、[EnsembleTR](https://github.com/gymrek-lab/EnsembleTR)（多 caller 合并！）、[strdenovotools/MonSTR](https://github.com/gymreklab/strdenovotools)。

### 3.4 ONT 侧对照（一览）

| 工具 | 说明 | 出处 |
|---|---|---|
| **NanoSatellite** | 直接用 raw squiggle 信号 DTW 分型（GC-rich 受益） | [Genome Biology 2019](https://pubmed.ncbi.nlm.nih.gov/31727106) |
| **HMMSTR** | 靶向 panel 的 profile-HMM 拷贝数估计（Boyle lab） | [NAR 2025](https://pubmed.ncbi.nlm.nih.gov/38746091) |
| **Mahi-mahi** | ONT 官方：长度+组成+中断+甲基化一次性输出（TRGT 的 ONT 对位） | [medRxiv 2026](https://pubmed.ncbi.nlm.nih.gov/42996439)、[nanoporetech/mahi-mahi](https://github.com/nanoporetech/mahi-mahi) |
| **Medaka Tandem** | ONT 官方轻量分型 | Genome Biology 2026 基准中 ONT 序列级第一（97.4%） |

---

## 4. 论文精读：可借鉴的算法优化思路

按主题归纳（每条注明出处），这些是"如果我要优化 CCS-STR，业界已经试过什么"。

### 4.1 读段处理与局部重比对

- **侧翼锚定（flanking anchoring）**：所有可靠 caller 的共同前提——要求读段在重复两侧都有比对锚，避免重复内部比对歧义（HipSTR 2017、tdg 2019、TRGT 2024 一脉相承）。tdg 额外用 **mismap probability < 1e-6** 过滤 + 屏蔽比对中的低复杂度小写区（[Genome Biology 2019](https://pubmed.ncbi.nlm.nih.gov/30890163)）。
- **按等位基因聚类 → per-allele consensus**：TRGT 用成对编辑距离 + Ward 聚类把读段分到两个等位基因，再对每簇取多数共识；STRdust 用长度加权 Levenshtein + DBSCAN(k-mer 组成) 双通道并互为 QC（[TRGT NB 2024](https://www.nature.com/articles/s41587-023-02057-3)；[STRdust Btae268](https://doi.org/10.1093/bioinformatics/btae268)）。**这是 HiFi 高一致读数下的红利：共识把单读误差抹平后再做下游决策。**
- **"切除重复"的人工参考重比对**：STRdust sensitive 模式把参考中的重复序列挖掉后重比对，专治"读段 clip 穿过扩增"的情形（[STRdust 2024](https://doi.org/10.1093/bioinformatics/btae268)）。
- **POA（partial order alignment）精修候选单倍型**：LongTR 用 spoa 对读簇取共识，INEXACT_ALLELE 标记哪些候选是 POA 近似出来的（[LongTR GB 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11229021)）。
- **局部组装兜底**：otter 靶向局部组装，绕开参考比对伪影；TREAT 证明它能救回 GC-rich/G-四联体位点的覆盖度丢失（[TRESI GR 2024](https://genome.cshlp.org/content/34/11/1942)）。

### 4.2 基因分型模型

- **读类别分解**：GangSTR 首创 enclosing/spanning/FRR(完全重复)/flanking 四类读的联合似然（[NAR 2019](https://pubmed.ncbi.nlm.nih.gov/31194863)）；ExpansionHunter 用序列图把各类读统一 realign 后拼接路径分型（[Bioinformatics 2019](https://pubmed.ncbi.nlm.nih.gov/31134279)）。**优化启示：in-repeat 读是超读长扩增的唯一证据源，丢弃它=丢召回。**
- **homopolymer 感知错误模型**：LongTR 的 Dindel 式 HMM 显式建模几何分布的 indel 错误，参数为 HiFi 调过默认（`--stutter-align-len` 可选开）；HiFi 残余错误确实集中于 homopolymer/STR 区的实证见 [GPB 2025 综述](https://pmc.ncbi.nlm.nih.gov/articles/PMC12257948)。
- **stutter 模型作为可学习参数**：HipSTR 用 EM 逐位点学 stutter（PCR 伪影模型）——HiFi 无 PCR stutter，这套模型的对应物是"聚合酶 indel 偏好"，是 CCS 错误模型可以借鉴的方向。
- **免目录发现 vs 目录分型**：EH Denovo/STRling 走 k-mer 离群路线（[GB 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02017-z)、[GB 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753380)）；tdg 走全基因组重复表扫描+重要性排序。两条路线与目录驱动分型（TRGT）互补。

### 4.3 超读长与复杂 motif

- **motif 表达式/HMM 定义复杂重复**：TRGT 允许用正则式 motif 表达或 HMM 定义重复结构（如 RFC1），并对每等位基因 consensus 拆解出 **motif 组成与中断（interruptions）**——致病性往往由"中断模式"而非纯长度决定（[TRGT NB 2024](https://www.nature.com/articles/s41587-023-02057-3)）。
- **长度≠一切的实证**：TRVZ/REViewer 一类的 pileup 可视化专门用于发现"分型长度对但组成错"的案例（[REViewer 原理文](https://www.illumina.com/science/genomics-research/articles/reviewer-alignments-short-reads-long-repeat.html)）。

### 4.4 定相与家系

- **近邻 SNV 锚定定相**：STRkit 的招牌——现场 call 局部 SNV、用它们把读分配到单倍型，无需 pedigree 或外部 phasing（[GR 2026](https://pubmed.ncbi.nlm.nih.gov/41539721)）。HiFi 读长让"SNV+STR 同读"成为可能，这是短读做不到的优化点。
- **Trio de novo**：TRGT-denovo 把子代读重比对到等位基因 consensus，用 parental-upper-bound 比对分统计筛 de novo，81k→137 个候选（[bioRxiv 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11275785)）；Gymrek 侧有 MonSTR（基于 GangSTR）。

### 4.5 体细胞不稳定与嵌合

- **读级 divergence 统计**：trgt-instability 用"读 vs 等位基因 consensus 的长度归一化编辑距离"构造 15-bin 不稳定谱，Dirichlet-multinomial 拟合对照基线后做检验，支持 `--adaptive-fdr-q` 提前终止 bootstrap（[PacificBiosciences/trgt-instability](https://github.com/PacificBiosciences/trgt-instability)）。**这是"CCS 读质量好 → 单读误差小 → 体细胞信号可分辨"的直接应用。**

### 4.6 CCS/读生成层（你的上游）

- **DeepConsensus**：gap-aware transformer（GATE）+ alignment loss，pbccs 之上再削 42% 残余错误；固定下机数据量时 Q20/Q30/Q40 读产出 +9%/+27%/+90%（[Nature Biotechnology 2023](https://www.nature.com/articles/s41587-022-01435-7)；⚠️ 社区资料常误写为 Nature Methods）。
- **AlphaEvolve 优化训练**（2026-04 PacBlog）：banded AlignmentLoss 经演化搜索调优后，实测达 Q30 的读比例 47.9%→53.2%（[PacBlog](https://www.pacb.com/blog/improving-hifi-sequencing-accuracy-with-google-deepconsensus-and-alphaevolve)）。
- **ccs 参数与 DeepConsensus 的交互**：ccs 默认丢弃预测 QV<20 的分子，DeepConsensus 无法救回被丢弃者——官方建议按其 quick-start 调 ccs 参数（[google/deepconsensus](https://github.com/google/deepconsensus)）。
- **`--by-strand` 与 QV 分箱 FAQ**：[ccs.how](https://ccs.how) 文档讨论了 passes 数与实测准确率的关系——对"STR 区单读错误率"的建模是直接输入。

### 4.7 多 caller 一致（ meta 层）

- Nature Protocols 2025 的 TR-关联研究指南明确推荐**跨多个 caller 取 consensus 基因型**再进 GWAS（[Nature Protocols 2025](https://experiments.springernature.com/articles/10.1038/s41596-025-01231-y)）；短读时代已有 EnsembleTR 做同样的事。

---

## 5. 基准与真值集

### 5.1 真值集（回归测试的"地基"）

| 资源 | 内容 | 出处 |
|---|---|---|
| **GIAB HG002 TR benchmark** ★标准 | 124,728 个小变异 + 17,988 个大变异，附表示归一化比较算法（处理不同工具的等位基因编码差异） | [English et al., Nature Biotechnology 2024](https://www.nature.com/articles/s41587-024-02225-z) |
| Broad 全基因组 TR truth set | 146,640 个变异（68,304 扩增 / 68,993 收缩），Mendelian 过滤 | [Weisburd et al., bioRxiv 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10197592)（未正式发表） |
| Variation clusters | 把 per-locus 等位基因折叠成可比单元，供群体数据库/GWAS 用 | [AJHG 2026](https://www.sciencedirect.com/science/article/pii/S0002929726001503) |
| HG002 T2T 组装 | TRGT/LongTR/STRkit 论文共同的对照锚 | — |

配套比较工具：[truvari](https://github.com/ACEnglish/truvari)（ACEnglish，即 tdb 作者）。

### 5.2 基准测试结论（汇总）

| 基准 | 平台 | 结论要点 |
|---|---|---|
| TRGT 论文（自报） | HiFi | Mendelian 98.38%(±1bp) > tdg 96.08% > Straglr 90.34% > GangSTR(短读) 72.90% |
| LongTR 论文（自报） | HiFi | 序列一致率 84.7%；>500bp 长重复与 null allele 是 TRGT 短板 |
| STRkit 论文（自报） | HiFi+ONT | F1 0.9631 (HiFi)；**TRGT 在其基准中速度最快** |
| TREAT/otter 论文 | 5–30x, SeqII/Revio/ONT | **<500bp PacBio 最佳，≥1000bp ONT 更好**；0.06% TR 覆盖度丢失（GC-rich/G4 富集）可用非 HiFi 读救回 |
| [Genome Biology 2026 大横评](https://link.springer.com/article/10.1186/s13059-026-04210-y)（Aliyev et al.） | ONT | 盘点 25 个工具、实测 7 个：Medaka Tandem 97.4% > LongTR 96.6% > ATaRVa 91.7%；含选型指南 |

**选型速查**：全基因组标准流程 → **TRGT**；追求长扩增/序列解析单倍型 → **+LongTR**；关联研究/需定相 → **+STRkit**；临床靶向 → PureTarget(=TRGT 引擎) + STRchive 目录；速度敏感的大 cohort 扫描 → STRdust / ATaRVa / tdg。

---

## 6. 重复目录与注释资源

- **[STRchive](https://github.com/dashnowlab/STRchive)**（[strchive.org](https://strchive.org)，[Genome Medicine 2025](https://doi.org/10.1186/s13073-025-01454-4)）：动态维护的致病 STR 目录，**TRGT 官方已弃用自带 BED 改推 STRchive**；含 genotyper-ready BED（最大位点 ~4kb）、UCSC track、等位基因频率。
- **Adotto catalog** 与 **Broad TR catalog**：全基因组位点目录（TRGT README 推荐），配套 Weisburd et al. 的 variation clusters。
- **全基因组重复图谱**：Chiu et al. 用 272 人长读数据建了 **>5M 个 TR 位点**的目录（[Nature Communications 2026](https://www.nature.com/articles/s41467-025-66153-5)）——给 TRGT/LongTR 喂位点的好素材。
- **[pathSTR](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610575)**（Genome Research 2024）：医学相关 TR 的正常人群等位基因分布（判"异常"的基线）。
- 经典起点：[Tandem Repeats Finder](https://github.com/Benson-Genomics-Lab/TRF)（TRF）与 UCSC simpleRepeat/microsat 表（tdg 直接吃这两种格式）。
- 队列层：**[tdb](https://github.com/ACEnglish/tdb)**（TRGTdb）：把多样本 TRGT VCF 变成紧凑数据库，`create/merge/bigmerge/query`，合并时建议每线程 5–10GB 内存。

---

## 7. 博客、官方文档与社区资源

### PacBlog（按实用度排序）

1. [Sequencing 101: Tandem repeats](https://www.pacb.com/blog/sequencing-101-tandem-repeats)（2023-11）——最实用的一篇：STR vs VNTR 的选择逻辑 + TRGT/TRVZ 工作流（长度→组成→甲基化→可视化）。
2. [Improving HiFi accuracy with DeepConsensus and AlphaEvolve](https://www.pacb.com/blog/improving-hifi-sequencing-accuracy-with-google-deepconsensus-and-alphaevolve)（2026-04）——上游读质量最新进展。
3. [A roadmap to expand PureTarget repeat expansion panels](https://www.pacb.com/blog/a-roadmap-to-expand-puretarget-repeat-expansion-and-carrier-panels-enabling-high-throughput-clinical-research-solutions)（2025-05）——靶向 panel 20→38 基因、96-plex、~707x 覆盖度。
4. [Beyond the Bench: repeat expansion diseases with HiFi](https://www.pacb.com/blog/beyond-the-bench-exploring-the-terrain-of-repeat-expansion-diseases-with-hifi-sequencing)（2026-03）。
5. TRGT 讲座视频（Dolzhenko 主讲，含 FMR1/mosaic/甲基化演示）：[PacBio 视频页](https://programs.pacb.com/l/1652/2025-08-27/456yrx)。

### 官方文档/教程

- [TRGT GitHub 文档](https://github.com/pacificbiosciences/trgt)（DevNet 已并入）：教程+MyBinder notebook、[TRVZ 图解读](https://github.com/pacificbiosciences/trgt/blob/main/docs/trvz-plots.md)、PureTarget 分析指南、spanning BAM 说明。
- [ccs.how](https://ccs.how)：ccs 官方用户指南（参数、QV 语义、准确率-vs-passes FAQ）。
- [google/deepconsensus](https://github.com/google/deepconsensus) quick start：ccs 参数与 DeepConsensus 的配合。
- SMRT Link v25.1/25.2 用户指南 PDF（pacb.com）：Tandem Repeat 分析 App 说明。
- [HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL)：官方 WGS 流程（含 TRGT 模块），生产化参考。
- SIB 瑞士 NGS 培训：[Repeat expansion analysis of PacBio data](https://sib-swiss.github.io/NGS-longreads-training/2022.3/course_material/group_work/project2)（2022 版，TRGT 出现前的老流程，可作历史参考）。
- 可公开测试数据：`downloads.pacbcloud.com/public/dataset/RepeatExpansionDisorders_NoAmp/` 与 `.../PureTargetRE/RFC1/`——调参现成素材。

### 第三方博客/社区

- [ASeq: Repeat Expansions – Long and Short Read Sequencing](https://aseqer.com/p/repeat-expansions-long-and-short)（2023-12）：HG002 FMR1 实测——HiFi pileup 肉眼可数 ~30 个 CCG，ExpansionHunter 给 31，ONT 比对明显跑偏。
- [Golden Helix: Long Read Part 3 – Repeat Expansions](https://www.goldenhelix.com/blog/long-read-part-3-repeat-expansions)（2026-05）：TRGT 输出 → VarSeq 注释（STRchive）→ 报告的下游消费视角。
- [AWS HealthOmics 基准 PacBio WGS 流程](https://aws.amazon.com/blogs/publicsector/benchmarking-pacbio-whole-genome-sequencing-variant-pipeline-analysis-with-aws-healthomics-workflows)（2025-07）：含 TRGT 模块的端到端耗时/成本（HG002 全流程 8.67h/$21，g5.2xlarge）。
- [Illumina: REViewer 原理](https://www.illumina.com/science/genomics-research/articles/reviewer-alignments-short-reads-long-repeat.html)（2020）：graph realign 可视化如何暴露分型偏差。
- GitHub issue 实战：[TRGT #93](https://github.com/PacificBiosciences/trgt/issues/93)（多样本/多 run 合并后跑 TRGT 的标准做法）；pbmm2 相关 issue（社区报告过 `min-gap-comp-id-perc` 偏高导致重复 spanning read 不比对——排查时用原生 minimap2 对照）。
- 目录站：[long-read-tools.org](https://long-read-tools.org)（1000+ 工具）、[B-UMMI/long-read-catalog](https://github.com/B-UMMI/long-read-catalog)、[Awesome-Bioinformatics](https://github.com/danielecook/Awesome-Bioinformatics)。

### 中文资源（检索确认，未逐篇深读）

- [三代测序时代的临床科研（知乎）](https://zhuanlan.zhihu.com/p/553831995)：30+ 重复扩增疾病与长读长必要性。
- [Subreads, HiFi 与 CCS 深度解析（知乎）](https://zhuanlan.zhihu.com/p/1910281983575164389)：CCS 读生成机制。
- **注意**：中文世界尚无成体系的 HiFi-STR 管线教程（TRGT 实战类），这本身是个内容空白。

---

## 8. GitHub 仓库总表

| 仓库 | 语言 | 版本/状态（2026-11） | Star≈ | 一句话 |
|---|---|---|---|---|
| [PacificBiosciences/trgt](https://github.com/PacificBiosciences/trgt) | Rust | v5.1.0，活跃 | 146 | ★HiFi TR 分型事实标准（长度+组成+甲基化+TRVZ） |
| [PacificBiosciences/ccs](https://github.com/PacificBiosciences/ccs) | C++ | v8.x（develop 分支活跃） | 127 | Subreads→HiFi 共识 |
| [google/deepconsensus](https://github.com/google/deepconsensus) | Python/TF | 1.2.0 | 265 | CCS 后 transformer 纠错 |
| [PacificBiosciences/pbmm2](https://github.com/PacificBiosciences/pbmm2) | C++ | v26.x，活跃 | 224 | 官方 minimap2 前端（HIFI preset） |
| [PacificBiosciences/trgt-instability](https://github.com/PacificBiosciences/trgt-instability) | — | 预发布 | 7 | TRGT 输出→体细胞不稳定检验 |
| [PacificBiosciences/HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) | WDL | v4.0.0 | — | 官方生产流程（含 TRGT） |
| [ACEnglish/tdb](https://github.com/ACEnglish/tdb) | Python | 活跃 | 9 | TRGT VCF→队列数据库 |
| [dashnowlab/STRchive](https://github.com/dashnowlab/STRchive) | 数据+Snakemake | 活跃 | 30 | 致病 STR 动态目录（TRGT 官方推荐） |
| [gymrek-lab/longtr](https://github.com/gymrek-lab/longtr) | C++ | v1.2 | 41 | HipSTR 改造的 HiFi/ONT 分型器 |
| [davidlougheed/strkit](https://github.com/davidlougheed/strkit) | Python+Rust | ≥v0.25.0，活跃 | 22 | SNV 锚定读级分型+CI；GPL-3.0 不限平台 |
| [wdecoster/STRdust](https://github.com/wdecoster/STRdust) | Rust | v0.16.x | 28 | 极速分型（CIGAR walk / 切除重比对双模式） |
| [bcgsc/straglr](https://github.com/bcgsc/straglr) | Python | v1.3.0 | 174 | GMM 聚类分型+de novo 发现 |
| [mcfrith/tandem-genotypes](https://github.com/mcfrith/tandem-genotypes) | Python | 1.5.x | 52 | LAST 比对驱动的免目录扫描 |
| [holstegelab/treat](https://github.com/holstegelab/treat) | — | 活跃 | — | 端到端 TR 分析（含 otter 局部组装） |
| [tfwillems/HipSTR](https://github.com/tfwillems/HipSTR) | C++ | 2020-04 最后发布 | — | 短读经典（LongTR 前身） |
| [gymreklab/GangSTR](https://github.com/gymreklab/gangstr) | C++ | — | — | 短读四类读似然分型 |
| [Illumina/ExpansionHunter](https://github.com/Illumina/ExpansionHunter) | C++ | 活跃（v5） | 220 | 短读序列图分型（catalog 格式事实标准） |
| [Illumina/REViewer](https://github.com/Illumina/REViewer) | C++ | 2026-04 已归档 | 92 | 短读 pileup 可视化（TRVZ 的前身概念） |
| [gymrek-lab/EnsembleTR](https://github.com/gymrek-lab/EnsembleTR) | Python | — | — | 多 caller VCF 合并 |
| [ACEnglish/truvari](https://github.com/ACEnglish/truvari) | Python | 活跃 | — | VCF 比较/基准（SV+TR） |
| [Jong-hun-Park/trviz](https://github.com/Jong-hun-Park/trviz) | Python | 活跃 | — | 重复序列分解/编码/可视化库（区别于 TRGT 内置的 TRVZ） |
| [Benson-Genomics-Lab/TRF](https://github.com/Benson-Genomics-Lab/TRF) | C | 经典 | — | 重复发现元老（Straglr 依赖它） |

> Star 数为抓取时近似值，仅用于判断社区体量，不代表质量。

---

## 9. 对「CCS STR 部分」的可落地优化清单

结合上述调研，按投入产出排序：

1. **读生成层**
   - 若还没开 DeepConsensus：先开（残余错误 −42%，Q40 级读 +90%），并按 [quick start](https://github.com/google/deepconsensus) 调 ccs 参数（避免"低 QV 读先被 ccs 丢弃、DeepConsensus 无从救回"）。
   - 关注 DeepConsensus 版本升级（AlphaEvolve 一代对重复/端粒区有专门优化）。
2. **比对层**
   - 用 pbmm2 `HIFI` preset；**排查"spanning read 数量异常"时用原生 minimap2 对照**（社区报告过 pbmm2 序列一致性参数导致重复读不比对）；确认 `--sort`、多样本多 run 合并后再跑 TRGT（[TRGT #93](https://github.com/PacificBiosciences/trgt/issues/93)）。
3. **分型层**
   - 主力 TRGT（v5.x）+ **STRchive 目录**（自带 BED 已弃用）。
   - 对 >500bp 扩增与 null allele 疑例用 LongTR 复核；对需要定相/关联分析的场景加 STRkit。
   - 重要位点（临床/发表）用 **otter 局部组装或 T2T 组装**交叉验证；GC-rich/G4 位点是覆盖度丢失高发区。
4. **队列与体细胞层**
   - 多样本用 tdb 建 cohort 库；体细胞不稳定用 trgt-instability（留意其预发布状态）。
   - 关联研究遵循 Nature Protocols 2025：多 caller consensus 基因型。
5. **基准回归**
   - 建立回归测试：GIAB HG002 TR truth set + truvari；任何 ccs/比对参数改动后跑一遍，量化对 STR 位点 Mendelian 一致率/序列一致率的影响——这是"CCS 层优化是否伤害 STR"的最直接证据。
6. **数据素材**
   - 用 PacBio 公开数据集（RepeatExpansionDisorders_NoAmp、PureTargetRE/RFC1）做调参的固定考卷。

---

## 10. 调研中的"防坑"记录

以下名字在社区/文档中偶有流传，但**未能核实为真实存在的工具**（检索于 2026-11）：

- **ISLPR**、**STRful**、**TReCAT**：未找到任何论文/仓库。
- **NanoSaturn**：不存在；真实工具是 **NanoSatellite**（ONT 信号 DTW，Genome Biology 2019）。
- **STRtyper**：法医 PCR 试剂盒品牌，不是生信工具。
- **tandem-genotypes 的归属**：常被误传为 Eichler lab——实为 **Martin Frith lab**（Genome Biology 2019）。
- **TRGT 论文期刊**：是 **Nature Biotechnology 2024**（不是 Genome Research）；LongTR 是 **Genome Biology 2024**（不是 Bioinformatics）。

---

## 11. 参考链接汇总

**方法论文（长读）**：TRGT ([NB 2024](https://www.nature.com/articles/s41587-023-02057-3)) · LongTR ([GB 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11229021)) · STRkit ([GR 2026](https://pubmed.ncbi.nlm.nih.gov/41539721)) · STRdust ([B 2024](https://doi.org/10.1093/bioinformatics/btae268)) · Straglr ([GB 2021](https://pubmed.ncbi.nlm.nih.gov/34389037)) · tandem-genotypes ([GB 2019](https://pubmed.ncbi.nlm.nih.gov/30890163)) · TRiCoLOR ([GigaSci 2020](https://pubmed.ncbi.nlm.nih.gov/33034633)) · TREAT/otter ([GR 2024](https://genome.cshlp.org/content/34/11/1942)) · TRGT-denovo ([bioRxiv 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11275785)) · ATaRVa ([bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.05.13.653434v2)) · Mahi-mahi ([medRxiv 2026](https://pubmed.ncbi.nlm.nih.gov/42996439)) · HMMSTR ([NAR 2025](https://pubmed.ncbi.nlm.nih.gov/38746091))

**方法论文（短读背景）**：lobSTR ([eLife 2012](https://elifesciences.org/articles/00313)) · HipSTR ([NM 2017](https://pubmed.ncbi.nlm.nih.gov/28436466)) · ExpansionHunter ([B 2019](https://pubmed.ncbi.nlm.nih.gov/31134279)) · GangSTR ([NAR 2019](https://pubmed.ncbi.nlm.nih.gov/31194863)) · EH Denovo ([GB 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02017-z)) · STRling ([GB 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9753380))

**CCS/读质量**：CCS/HiFi 原始论文 ([NB 2019](https://www.nature.com/articles/s41587-019-0217-9)) · DeepConsensus ([NB 2023](https://www.nature.com/articles/s41587-022-01435-7)) · HiFi 综述 ([GPB 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12257948)) · [DeepConsensus/AlphaEvolve PacBlog (2026)](https://www.pacb.com/blog/improving-hifi-sequencing-accuracy-with-google-deepconsensus-and-alphaevolve)

**基准/真值/目录**：GIAB TR benchmark ([NB 2024](https://www.nature.com/articles/s41587-024-02225-z)) · Broad truth set ([bioRxiv 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10197592)) · Variation clusters ([AJHG 2026](https://www.sciencedirect.com/science/article/pii/S0002929726001503)) · 全基因组 TR 目录 ([NC 2026](https://www.nature.com/articles/s41467-025-66153-5)) · pathSTR ([GR 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610575)) · STRchive ([GM 2025](https://doi.org/10.1186/s13073-025-01454-4)) · 大横评 ([GB 2026](https://link.springer.com/article/10.1186/s13059-026-04210-y)) · 多 caller 指南 ([Nat Protocols 2025](https://experiments.springernature.com/articles/10.1038/s41596-025-01231-y))

---

*本文档由网络调研汇总而成，链接于 2026-11-17 核实可访问；工具版本号以当日 GitHub 页面为准。*
