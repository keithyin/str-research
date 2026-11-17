# 相关论文扩展调研

> 日期：2026-11-18
> 方法：以本目录 4 篇论文（Jeanjean 2025 / Li 2017 / HMMSTR 2025 / Ma 2011）为种子，用 Europe PMC 引文 API 挖掘其**参考文献与被引文献**，再辅以定向检索核实。共整理 ~55 篇新增相关论文。
> 标注说明：**[已收录]** = 已在 [pacbio-ccs-str-optimization-survey.md](pacbio-ccs-str-optimization-survey.md) 中覆盖；★ = 对你的 CCS-STR 优化最相关，见文末 Top 10。

---

## A. 沿 Jeanjean 2025（平台比较 / 错误谱 / stutter）出发

Jeanjean 的核心关切是「**测出 STR 的真实长度分布**」——其 78 篇参考文献里埋着一条完整的「错误建模与校正」文献线：

### A1. stutter/错误机制的直接前身与同类工作

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Zavodna et al., 2014** — The accuracy, feasibility and challenges of sequencing STRs using NGS platforms | PLoS One, [PMID 25436869](https://pubmed.ncbi.nlm.nih.gov/25436869) | Jeanjean 的直接前身：早期跨平台 STR 测序准确性系统实验 |
| **Raz et al., 2019** ★ — Short tandem repeat **stutter model inferred from direct measurement** of in vitro stutter noise | NAR, [PMID 30698816](https://pubmed.ncbi.nlm.nih.gov/30698816) | 用体外实验直接实测 stutter 并推断统计模型——HipSTR/GangSTR stutter 模型的实证基础；给「如何为噪声建模」提供了范式 |
| **Shinde et al., 2003** — Taq polymerase **slippage mutation rates** ((CA/GT)n, (A/T)n) | NAR, [PMID 12560493](https://pubmed.ncbi.nlm.nih.gov/12560493) | 经典定量：PCR 滑移率的原始测量（Li 2017 的滑移概率模型与之呼应） |
| **Woerner et al., 2020** ★ — Reducing noise and **stutter in STR loci with UMI** | Forensic Sci Int Genet, [PMID 33429137](https://pubmed.ncbi.nlm.nih.gov/33429137) | Jeanjean 的 UMI 思路在法医 STR 上的先行验证 |
| **Crysup et al., 2022** — UMI to improve allele calling in low-template mixtures | Forensic Sci Int Genet, [PMID 36462297](https://pubmed.ncbi.nlm.nih.gov/36462297) | UMI + 混合物分解（法医 MPS 场景） |
| **Daunay et al., 2019** — **低温等温扩增**大幅降低 stutter | NAR, [PMID 31584085](https://pubmed.ncbi.nlm.nih.gov/31584085) | 湿实验层面"消灭 stutter"的另一个流派（低温扩增） |
| **Stoler & Nekrutenko, 2021** — Sequencing **error profiles of Illumina** instruments | NAR Genom Bioinform, [PMID 33817639](https://pubmed.ncbi.nlm.nih.gov/33817639) | 短读错误谱参考数据（Jeanjean 错误模型的对照） |
| **Laehnemann et al., 2015** — Denoising DNA deep sequencing data（综述） | Brief Bioinform, [PMID 26026159](https://pubmed.ncbi.nlm.nih.gov/26026159) | 高通量测序去噪的系统综述 |

### A2. 分子条码/共识的起源经典（UMI/duplex 的根）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| Kivioja et al., 2011 — **Counting absolute molecules using UMI** | Nat Methods, [PMID 22101854](https://pubmed.ncbi.nlm.nih.gov/22101854) | UMI 概念的原始出处 |
| Kinde et al., 2011 — **Detection of rare mutations** with massively parallel sequencing | PNAS, [PMID 21586637](https://pubmed.ncbi.nlm.nih.gov/21586637) | Safely-Seq：分子条码纠错开端 |
| Schmitt et al., 2012 — **Duplex sequencing** 检测超稀有突变 | PNAS, [PMID 22853953](https://pubmed.ncbi.nlm.nih.gov/22853953) | Duplex 共识的鼻祖（Jeanjean Type III 文库的源头） |
| Hiatt et al., 2013 — **Single molecule MIP** 靶向高精度低频变异 | Genome Res, [PMID 23382567](https://pubmed.ncbi.nlm.nih.gov/23382567) | MIP 靶向 + 单分子共识 |
| **Bae et al., 2023** ★ — **CODEC**: single **duplex DNA sequencing** | Nat Genet, [PMID 37106072](https://pubmed.ncbi.nlm.nih.gov/37106072) | 现代 duplex 全基因组测序方法（Jeanjean duplex 结论的方法学近亲） |
| Wang TT et al., 2019 — **High efficiency error suppression** for low-frequency variants | NAR, [PMID 31127310](https://pubmed.ncbi.nlm.nih.gov/31127310) | 错误抑制的通用框架 |
| Fang et al., 2014 — **Reducing INDEL calling errors** (GIAB 系) | Genome Med, [PMID 25426171](https://pubmed.ncbi.nlm.nih.gov/25426171) | indel 错误削减的经典工程实践 |

### A3. 长读层错误（与你 CCS 工作直接相关）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Zhang et al., 2020** ★ — **A comprehensive evaluation of long read error correction methods** | BMC Genomics, [PMID 33349243](https://pubmed.ncbi.nlm.nih.gov/33349243) | 长读纠错方法横评——CCS/DeepConsensus 之外的方法族地图 |
| **Sacristán-Horcajada et al., 2021** — **ARAMIS**: from systematic errors of NGS long reads to accurate assemblies | Brief Bioinform, [PMID 34013348](https://pubmed.ncbi.nlm.nih.gov/34013348) | 长读**系统性错误**的刻画与利用（对「为什么 consensus 偏向插入」这类问题有直接参考价值） |
| [已收录] GPB 2025 HiFi 综述、DeepConsensus NB 2023、ccs NB 2019 | — | 见调研总文档 §4.6 |

### A4. 微卫星突变生物学（背景，理解"真值分布"长什么样）

Sainudiin 2004 Genetics（微卫星突变模型，[PMID 15454552](https://pubmed.ncbi.nlm.nih.gov/15454552)）· Ellegren 2004 Nat Rev Genet（[PMID 15153996](https://pubmed.ncbi.nlm.nih.gov/15153996)）· Mirkin 2007 Nature（可扩增重复与疾病）· Eckert 2009（每位点突变差异）· Kelkar 2010 GBE（"什么是微卫星"的定义性工作）· Sun 2012 Nat Genet（微卫星直接测人类突变率）· 酵母经典 Strand 1993 / Wierdl 1997 / Sia 1997 · **McGinty et al., 2025 Nat Commun** ★（[PMID 41365873](https://pubmed.ncbi.nlm.nih.gov/41365873)）——「简单重复的内在不稳定性塑造进化稳定的长度分布」，**引用了 Jeanjean**，把"真实等位基因分布"上升为群体遗传学问题。

---

## B. 沿 Li 2017 AmpSeq-SSR（大规模 SSR 分型）出发

33 篇被引文献基本就是「**靶向测序 SSR 分型**」这个细分领域的地图，按价值筛选：

### B1. 计算方法/工具（对管线最有用）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Liu et al., 2024** ★★ — **Seq2Sat and SatAnalyzer** toolkit: comprehensive microsatellite genotyping | Mol Ecol Resour, [PMID 38289068](https://pubmed.ncbi.nlm.nih.gov/38289068) | 引用 Li 的最新一站式微卫星分型工具箱（SatAnalyzer），生态位与 Li 的计算流程重叠且更新 |
| **Huo et al., 2021** ★ — An integrated strategy for **target SSR genotyping tolerating nucleotide variations** | BMC Bioinformatics, [PMID 34496768](https://pubmed.ncbi.nlm.nih.gov/34496768) | 解决 AmpSeq 类方法的痛点：SSR 内碱基变异导致等位基因分型失败的容错算法 |
| **Guo et al., 2020** — **MultiplexSSR**: pipeline for developing multiplex SSR-PCR assays from resequencing data | Ecol Evol, [PMID 32211176](https://pubmed.ncbi.nlm.nih.gov/32211176) | 多重 PCR 引物设计的自动化（Li 的 3105 对引物人工成本 → 流程化） |
| **Lepais & Weir, 2020** ★ — **Fast sequence-based microsatellite genotyping** development workflow | PeerJ, [PMID 32411534](https://pubmed.ncbi.nlm.nih.gov/32411534) | 37 引，与 AmpSeq-SSR 并列的主流 SSR-测序分型工作流（含 stutter 处理细节） |

### B2. 平行技术路线（扩音器 vs 捕获 vs 降样本复杂度）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Yang et al., 2019** ★★ — **Target SSR-Seq**: SSR genotyping associated with perfect SSRs | Front Plant Sci, [PMID 31105728](https://pubmed.ncbi.nlm.nih.gov/31105728) | 39 引，Li 的最直接竞品/姊妹技术（同样多重 PCR+测序思路） |
| **Li T et al., 2019** — Application of high-throughput amplicon sequencing-based SSR genotyping（辣/甜椒） | BMC Genomics, [PMID 31159719](https://pubmed.ncbi.nlm.nih.gov/31159719) | AmpSeq-SSR 的应用复制 |
| **Baggett et al., 2021** — **De novo SSR identification → targeted sequencing** fingerprints Sorghum | PLoS One, [PMID 33684158](https://pubmed.ncbi.nlm.nih.gov/33684158) | 「先挖 SSR 再靶向」的完整闭环，衔接 de novo 目录问题 |
| **Shen et al., 2021** / **Guo et al., 2019** — **GBTS**（目标测序基因分型面板）+ KASP | Front Plant Sci, [PMID 34149754](https://pubmed.ncbi.nlm.nih.gov/34149754) / Mol Breed | 中国育种界的主流替代路线：固相捕获面板（GBTS）vs 扩增子面板（AmpSeq）——两条路线的成本/精度对比值得读 |
| **Deng et al., 2024** — GenoBaits®WheatplusEE targeted capture panel | Theor Appl Genet, [PMID 38291310](https://pubmed.ncbi.nlm.nih.gov/38291310) | 捕获面板在品种鉴定上的最新实践 |

### B3. 质量/污染/数据库（呼应 Ma 2011）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Gao et al., 2023** ★ — **Carryover contamination-controlled amplicon sequencing** workflow | Microbiol Spectr, [PMID 37098913](https://pubmed.ncbi.nlm.nih.gov/37098913) | Ma 2011「污染控制」问题的现代测序版（carryover 残留对 amplicon 测序的影响与控制流程） |
| Jiang et al., 2020 — **PIDS** 植物 DNA 指纹数据库管理系统 | Genes, [PMID 32235513](https://pubmed.ncbi.nlm.nih.gov/32235513) | Li 指纹图谱思路的数据库化 |
| Azizi et al., 2021 — Plant variety/cultivar identification 技术综述 | J Biosci, [PMID 34544910](https://pubmed.ncbi.nlm.nih.gov/34544910) | 全景定位 AmpSeq-SSR 在品种鉴定谱系中的位置 |
| Putman & Carstens, 2014 — **Challenges in analysis and interpretation of microsatellite data** | Ecol Evol, [PMID 25540699](https://pubmed.ncbi.nlm.nih.gov/25540699) | 微卫星数据分析的通用陷阱清单（等位基因丢失、null allele 等） |
| Miah et al., 2013 —水稻 SSR 与抗病育种综述 | Int J Mol Sci, [PMID 24240810](https://pubmed.ncbi.nlm.nih.gov/24240810) | Li 论文的水稻语境 |

### B4. 微卫星功能生物学（背景）

Bagshaw & Pitt 2017 GBE「Functional mechanisms of microsatellite DNA in eukaryotic genomes」（110 引，[PMID 28957459](https://pubmed.ncbi.nlm.nih.gov/28957459)）——Li/Jeanjean 都引它，理解 SSR 为何重要的一站式综述。

---

## C. 沿 HMMSTR 2025（靶向长读 STR 分型）出发

被引 9 篇几乎全是 2025–2026 的最新工作——这是当前最活跃的方向：

### C1. 方法论文（最高优先）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Su et al., 2026** ★★★ — **MASTR-seq**: multiplexed analysis of STRs with sequencing | Cell Rep Methods, [S2667-2375(26)00041-X](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(26)00041-X)（preprint [PMC11092654](https://pmc.ncbi.nlm.nih.gov/articles/PMC11092654)） | **Cas9 介导、尺寸选择、PCR-free** 的多重 STR 测序——把 Li 的"多重靶向"、Jeanjean 的"PCR-free 抗 stutter"、HMMSTR 的"靶向长读"三个关键词合到了一个方法里；引用了 HMMSTR，是用户四篇论文的技术交汇点 |
| **Benarroch et al., 2026** ★★★ — Toward the clinical application of long-read sequencing in repeat-expansion disorders | **Nat Genet** 58:2087–2103 | 2026 年最新的长读-重复扩增临床应用权威综述（CRISPR 靶向 + 长读的完整图景），引用 HMMSTR |
| **Lee et al., 2026** ★★ — Optimized **Cas9-enriched nanopore** sequencing workflow for clinical diagnosis of REDs | Adv Sci, [PMID 42667167](https://pubmed.ncbi.nlm.nih.gov/42667167) | nCATS 类靶向富集的流程优化（HMMSTR 实验策略的下一代） |
| **Fang et al., 2022** ★★ — **DeepRepeat**: direct quantification of STRs on nanopore **signal data** | Genome Biol 23:108, [DOI](https://link.springer.com/article/10.1186/s13059-022-02670-6)，代码 [WGLab/DeepRepeat](https://github.com/WGLab/DeepRepeat) | HMMSTR 论文里对比的信号法 caller：离子电流→RGB 图像→CNN 判重复；绕过碱基调用，专治 GC-rich/超长 |
| **Liu et al., 2020** ★ — Genome-wide detection of STR expansions by long-read sequencing（**RepeatHMM 长读版 + RepeatHMM-scan**） | [PMID 33371889](https://pubmed.ncbi.nlm.nih.gov/33371889) | RepeatHMM 的长读扩展：用 21 个长读数据建立 432,604 个 STR 的正常范围——「正常范围」资源思路与 pathSTR 相通 |
| **Wang et al., 2023** — **TRcaller**: precise, ultrafast TR variant detection | [PMC10390829](https://pmc.ncbi.nlm.nih.gov/articles/PMC10390829) | 靶向/扩增子场景的快速 TR 检测（与 USAT 同生态位） |
| **Ge et al., 2022** — **USAT**: universal STR allele toolkit（长度↔序列等位基因转换 + 比较可视化） | BMC Bioinformatics 23:497 | 法医/靶向扩增子 STR 的下游解析工具——你做等位基因序列比较时可能用得上 |

### C2. 短读时代被反复引用的基线（补齐工具谱系）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Liu et al., 2017** ★ — **RepeatHMM**（HMM 重复分型，WGLab） | Bioinformatics，代码 [WGLab/RepeatHMM](https://github.com/WGLab/RepeatHMM) | HMMSTR/tdg/Jeanjean 论文共同的对比基线；RepeatHMM→DeepRepeat 是同一实验室的两条路线 |
| **Bakhtiari et al., 2018** ★ — **adVNTR**: targeted genotyping of VNTRs（HMM，短读+PacBio CLR） | NAR Genom Bioinform，代码 [mehrdadbakhtiari/adVNTR](https://github.com/mehrdadbakhtiari/adVNTR) | TRGT 的 HMM 建模方式直接继承 adVNTR；后续 **adVNTR-NN**（2021，浅神经网络 18 s/VNTR）与 **Code-adVNTR**（2022 iScience，编码区 VNTR） |
| **Dashnow et al., 2018** ★ — **STRetch**: detecting and discovering pathogenic STR expansions | GigaScience, [PMC6102892](https://pmc.ncbi.nlm.nih.gov/articles/PMC6102892)，代码 [Oshlack/STRetch](https://github.com/Oshlack/STRetch) | PCR-free 短读 WGS 的全基因组扩展筛查（EHdn/STRling 的同族）——注意作者 Harriet Dashnow 就是后来 TRGT/STRchive 的作者 |
| Tankard et al., 2018 — **exSTRa**: detecting expansions of tandem repeats in cohorts | Bioinformatics | ECDF 离群图做队列扩展筛查（EHdn 的平行方案） |
| **Tang et al., 2017** ★ — **TREDPARSE**: probabilistic approach to profile STR disease alleles in WGS | PLOS Comput Biol | 概率模型扩展读长限制的早期实现（read-depth + 背景分布模型）；同组 **Tang et al. 2017 AJHG**「12,632 全基因组 STR 致病等位基因谱」[PMID 29100084](https://pubmed.ncbi.nlm.nih.gov/29100084) 是大规模 WGS STR 应用的标杆 |
| **Anvar et al., 2014** — **TSSV**: characterization of complex allelic variants in pure and mixed genomes | Bioinformatics 30(12):1651, [DOI](https://doi.org/10.1093/bioinformatics/btu068) | 复杂 STR 等位基因（含混合样本）的序列级解析，法医 MPS 常用 |
| Fungtammasan et al., 2015 ★ — Accurate typing of STRs from genome-wide sequencing data | Genome Res, [PMID 25823460](https://pubmed.ncbi.nlm.nih.gov/25823460) | 早期全基因组 STR typing 系统方法（含错误特征刻画） |
| **Bahlo et al., 2018** — 短读重复扩增检测综述（EH/exSTRa/STRetch/TREDPARSE 四工具对比） | [PMC6008857](https://pmc.ncbi.nlm.nih.gov/articles/PMC6008857) | 上面四篇的一页式导读 |

### C3. 综述与跨领域基准（2025–2026）

| 论文 | 出处 | 为什么相关 |
|---|---|---|
| **Maestri et al., 2025** ★★ — **Navigating triplet repeats sequencing**: concepts, challenges, perspective for HD | NAR, [PMID 39676657](https://pubmed.ncbi.nlm.nih.gov/39676657) | 11 引；以亨廷顿病为例的三联重复测序方法学综述（湿实验+生信+验证） |
| **Tran et al., 2026** ★ — Detection of STRs in the **cattle genome**: comparison of bioinformatic tools | BMC Genomics, [PMID 41882523](https://pubmed.ncbi.nlm.nih.gov/41882523) | **非人物种**的 STR 工具横评——如果你要处理水稻/作物 SSR（Li 的语境），这是"非人基因组上谁好用"的现成参考 |
| **Khojakulov et al., 2026** — Computational STR genotyping in a large Turkish neurodegeneration cohort | Int J Mol Sci, [PMID 42196324](https://pubmed.ncbi.nlm.nih.gov/42196324) | 大 cohort 计算 STR 分型的最新实战 |
| Rajan-Babu et al., 2021 ★ — Genome-wide sequencing as **first-tier screening** for STR expansions | Genome Med, [PMID 34372915](https://pubmed.ncbi.nlm.nih.gov/34372915) | 「WGS 一线筛查 + 靶向确认」的策略框架（Straglr 一作） |
| **Shin et al., 2017** ★ — **CRISPR-Cas9 targeted fragmentation** + selective sequencing for massively parallel microsatellite analysis | Nat Commun, [PMID 28169275](https://pubmed.ncbi.nlm.nih.gov/28169275) | CRISPR 靶向微卫星富集的鼻祖（MASTR-seq/nCATS 的前身） |
| Waalkes et al., 2018 — **Single-molecule MIP** capture for pan-cancer MSI | Clin Chem, [PMID 29632165](https://pubmed.ncbi.nlm.nih.gov/29632165) | 靶向捕获 + 单分子共识在 MSI 上的临床级实现 |
| Gallon et al., 2019（×2）— 测序 MSI 检测（6 标记高通量 / 体质性 MMR 缺乏症） | Hum Mutat | MSI 临床测序检测的验证框架 |
| [已收录] Maestri 同类综述 Chaushevska 2025 Front Genet；ONT 侧 Mahi-mahi/NanoSatellite/WarpSTR/ATaRVa | — | 见调研总文档 |

---

## D. 沿 Ma 2011（PCR 污染控制）出发

Ma 2011 的主题（amplicon 假阳性/污染）在测序时代主要由「B3」节承接：

- **Gao 2023**（carryover contamination-controlled amplicon workflow，见表 B3）★ —— 最直接的现代对应
- Li 2017 自身的「双平台一致性 + 二项检验剔除 error-prone 位点」也是一种污染/错误控制设计
- 若需要更深的湿实验线，Ma 2011 的引用网络（GMR 2011，引文多为大豆 MAS 应用）技术含量有限，建议跳过

---

## E. Top 10 推荐阅读（按你 CCS-STR 优化的相关性排序）

1. **MASTR-seq**（Su 2026, Cell Rep Methods）——Cas9+尺寸选择+PCR-free 的多重 STR 测序：你四篇论文关键词的合流点
2. **Benarroch 2026**（Nat Genet 58:2087）——长读重复扩增临床应用最新权威综述，快速补齐全景
3. **Raz 2019**（NAR）——stutter 噪声的实测与建模：错误模型设计的范式来源
4. **DeepRepeat**（Fang 2022, Genome Biol）——信号层 STR 检测：与 ccs/碱基层互补的另一维度
5. **Woerner 2020**（Forensic Sci Int Genet）——UMI 降 stutter 的实证：如果考虑在 amplicon 层做分子条码
6. **Liu 2020**（RepeatHMM 长读版）——432k STR 正常范围资源 + 长读扩展检测
7. **Seq2Sat/SatAnalyzer**（Liu 2024, Mol Ecol Resour）——植物/生态微卫星分型最新工具箱（若 SSR 分型是目标场景）
8. **McGinty 2025**（Nat Commun）——重复长度分布的群体遗传学：为"真实等位基因分布"提供理论背景
9. **Zhang 2020**（BMC Genomics）——长读纠错方法横评：定位 ccs/DeepConsensus 在纠错方法族中的位置
10. **Tran 2026**（BMC Genomics）——非人物种 STR 工具横评：非人基因组应用选型参考

---

## F. 检索方法与覆盖说明

- 种子：本目录 4 篇 PDF 的实际 PMID（Jeanjean: 40036507；Li: 28184437；HMMSTR: 39676678；Ma: GMR 2011，未做引文挖掘，理由见 §D）
- 引文数据：Europe PMC REST API（references + citations，全量拉取后人工筛选）
- 补充定向检索：DeepRepeat、MASTR-seq、USAT、RepeatHMM、adVNTR、STRetch、TSSV/exSTRa、TREDPARSE、Benarroch 2026 的出处与链接均经二次核实
- 未覆盖：纯临床病例报道、纯群体遗传学应用（如具体物种多样性调查）仅择要收录；Jeanjean 引用的癌症 MSI 临床文献（Kim 2013、Hause 2016、Boland 2010 等）未单独展开——若你做 MSI 方向可从 §A1/§C3 的 MSI 条目顺藤摸瓜
