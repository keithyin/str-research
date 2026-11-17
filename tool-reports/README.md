# PacBio CCS-STR 工具调研报告集

> 生成日期：2026-11-18 ｜ 上游文档：[pacbio-ccs-str-optimization-survey.md](../pacbio-ccs-str-optimization-survey.md)
> 每份报告统一按「**是什么 → 解决了什么问题 → 如何解决的**」展开，辅以输入输出要点、性能与局限、对 CCS-STR 优化的启示。

## 优先级划分

按与你「PacBio CCS 的 STR 部分优化」工作的相关度划分：

- **P0 核心必读**：直接构成 HiFi-STR 主链路的 6 个工具（读生成 2 + 比对 1 + 分型 3）。
- **P1 重要**：决定分型上限的资源与外延能力（目录、队列库、体细胞、组装式替代、极速扫描）。
- **P2 背景借鉴**：短读经典与衍生工具——它们的方法是长读工具的设计来源，借鉴思想而非直接使用。
- **P3 外围速览**：特定场景（发现、基准、可视化、对照平台、流程运维）才有用的工具，单文件合并速览。

## 报告索引

### P0 — 核心必读

| # | 报告 | 工具 | 一句话定位 |
|---|---|---|---|
| 1 | [P0-01-TRGT.md](P0-01-TRGT.md) | TRGT | HiFi TR 分型事实标准（长度+组成+甲基化+定相）；注意专有许可仅限 PacBio 数据 |
| 2 | [P0-02-LongTR.md](P0-02-LongTR.md) | LongTR | 改造 HipSTR 的 HiFi 分型器；序列级单倍型、>500bp 长扩增、null allele |
| 3 | [P0-03-STRkit.md](P0-03-STRkit.md) | STRkit | 近邻 SNV 锚定的读级分型与定相；读长偏置重加权；GPL-3.0 跨平台 |
| 4 | [P0-04-ccs.md](P0-04-ccs.md) | ccs / pbccs | Subreads→HiFi 共识；残余错误 92% 是 homopolymer indel（STR 的头号误差源） |
| 5 | [P0-05-DeepConsensus.md](P0-05-DeepConsensus.md) | DeepConsensus | gap-aware transformer 后纠错：残余错误 −42%、Q40 读 +90% |
| 6 | [P0-06-pbmm2.md](P0-06-pbmm2.md) | pbmm2 | 官方 minimap2 前端；HIFI preset 与强制参数决定了 spanning read 的判定 |

### P1 — 重要

| # | 报告 | 工具 | 一句话定位 |
|---|---|---|---|
| 7 | [P1-01-STRdust.md](P1-01-STRdust.md) | STRdust | 极速分型（CIGAR walk / 切除重比对双模式）；QC 标签进 VCF 的典范 |
| 8 | [P1-02-Straglr.md](P1-02-Straglr.md) | Straglr | GMM 分型 + de novo 扫描双模式；AIC 选等位基因数 |
| 9 | [P1-03-tandem-genotypes.md](P1-03-tandem-genotypes.md) | tandem-genotypes | LAST mismap 过滤 + 整数倍取整的免目录扫描（Frith lab，非 Eichler） |
| 10 | [P1-04-TREAT-otter.md](P1-04-TREAT-otter.md) | TREAT + otter | 靶向局部组装分型；G4/GC-rich 覆盖度丢失的实证与补救 |
| 11 | [P1-05-trgt-instability.md](P1-05-trgt-instability.md) | trgt-instability | 读-vs-consensus divergence + Dirichlet-multinomial 的体细胞不稳定检验 |
| 12 | [P1-06-tdb.md](P1-06-tdb.md) | tdb | TRGT VCF → Parquet 三表队列库（体积 −87.6%） |
| 13 | [P1-07-STRchive.md](P1-07-STRchive.md) | STRchive | 动态维护的致病 STR 目录（TRGT 官方推荐，替代自带 BED） |

### P2 — 背景借鉴

| # | 报告 | 工具 | 借鉴价值 |
|---|---|---|---|
| 14 | [P2-01-HipSTR.md](P2-01-HipSTR.md) | HipSTR | EM 学习 stutter 模型 + HMM + SNP 定相——LongTR 直系祖先 |
| 15 | [P2-02-ExpansionHunter.md](P2-02-ExpansionHunter.md) | ExpansionHunter | 自环节点序列图 + 复合位点整体建模；catalog 格式事实接口 |
| 16 | [P2-03-GangSTR.md](P2-03-GangSTR.md) | GangSTR | 四类读（E/S/FRR/F）联合似然——「读类别」框架源头 |
| 17 | [P2-04-ExpansionHunterDenovo.md](P2-04-ExpansionHunterDenovo.md) | EH Denovo | IRR 周期检测 + 离群统计的免目录发现 |
| 18 | [P2-05-STRling.md](P2-05-STRling.md) | STRling | 回收比对器丢弃读的 k-mer 分析（旋转不变 k-mer） |
| 19 | [P2-06-REViewer.md](P2-06-REViewer.md) | REViewer | 单倍型分辨 pileup 可视化（TRVZ 的短读前身） |
| 20 | [P2-07-TRGT-denovo.md](P2-07-TRGT-denovo.md) | TRGT-denovo | 全员对齐子代共识 + PU 统计的 trio de novo 检测（8.1万→137） |

### P3 — 外围速览

[**P3-外围工具速览.md**](P3-外围工具速览.md)（单文件 14 条）：lobSTR、TRiCoLOR、ATaRVa、HMMSTR、Mahi-mahi、NanoSatellite、EnsembleTR、MonSTR、truvari、trviz、TRF、pathSTR、HiFi-human-WGS-WDL。

## 建议阅读顺序

1. 先读 P0-01/02/03（三个分型器的方法对比——聚类共识 vs POA+HMM vs SNV 锚定 GMM 三条路线）；
2. 再读 P0-04/05（读生成层：ccs 错误谱 → DeepConsensus 如何压错误——这是你优化工作的上游杠杆）；
3. P0-06 + P1-04（比对层的坑与组装式替代）；
4. P1-05/06/07（体细胞、队列、目录——产品化所需的外围）；
5. P2 按兴趣抽读（每篇末尾都有「启示」小结）；
6. P3 需要时再查。

## 关键事实速记（调研中的核对结果）

- **TRGT 许可证**：PacBio 专有许可（`LICENSE.md` §3），仅限 PacBio 仪器产生的数据 → 解释了 STRkit（GPL-3.0）论文中「TRGT 不能跑 ONT」的说法。
- **DeepConsensus 论文期刊**：**Nature Biotechnology** 41(2):232–238（DOI 10.1038/s41587-022-01435-7）——多处资料误写为 Nature Methods。
- **tandem-genotypes 归属**：Frith lab（Genome Biology 2019），非 Eichler lab。
- **HiFi 残余错误谱**（CCS 论文）：92% homopolymer indel / 4.6% 其他 indel / 3.4% mismatch——一切 STR 错误模型设计的前提。
- **工具性能数字**均为各自论文自报口径，真值集与指标不同，不可直接横向比较；中性基准见 Genome Biology 2026 大横评与 STRkit 论文。
