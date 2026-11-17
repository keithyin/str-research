# LongTR

> **优先级**：P0（核心必读）｜ **类别**：TR/STR 基因分型 ｜ **维护方**：Gymrek lab（UCSD）+ NIST ｜ **状态**：活跃，v1.2（conda）
> **一句话**：把短读经典 HipSTR 的「单倍型 + HMM」体系改造成 HiFi/ONT-Duplex 版，主攻**序列级单倍型解析**与**长重复/null allele**，与 TRGT 互补。

## 是什么

LongTR（*Long-read TR genotyping*）是 Gymrek lab 对自家 HipSTR（短读 STR 分型器）的改造版，C++/HTSlib 实现，GPL-2.0，输入 HiFi 或 ONT Duplex 的比对 BAM + 重复目录，输出序列级单倍型的 VCF。

- 仓库：https://github.com/gymrek-lab/longtr （~41★，v1.2，README 自述仍在活跃开发）
- 论文：*LongTR: genome-wide profiling of genetic variation at tandem repeats from long reads*, **Genome Biology 2024**（Ziaei Jam 等）https://pmc.ncbi.nlm.nih.gov/articles/PMC11229021
- 生态：Gymrek lab 还有 TRTools、EnsembleTR（多 caller 合并）、MonSTR（de novo）。

## 解决了什么问题

1. **HipSTR 的比对假设在长读上失效**。HipSTR 的 STR 区模型假设「读相对单倍型最多只有一处 indel，且幅度是重复单元长度的整数倍」——这对 PCR stutter 主导的短读成立；但 LongTR 作者实测发现**长读的错误以单碱基 indel 为主、且出现在任意位置**，整数倍假设直接崩掉。
2. **长度之外要序列**：TRGT 之外的另一条路线是「直接给出等位基因的完整序列」，用于等位基因序列关联、null allele（整段重复缺失）检测。
3. **>500 bp 长重复与跨读不足的情形**：TRGT 依赖 spanning reads，某些截断读、低支持大插入会漏。
4. **大 cohort 的工程问题**：单核、内存小、直接读云端 URL 的 BAM。

## 如何解决的

### 1. 候选单倍型生成：两轮聚类 + POA 精修
- **Trimming**：把比对到 TR 位点的读裁剪为「重复 + 用户定义侧翼上下文」；
- **第一轮**：任何被 ≥2 条读、且 >20%（单样本）或 >5%（跨样本）支持的序列成为候选单倍型；
- **第二轮（补救漏掉的单倍型）**：当某样本 >25% 的比对读没有归属时，用排除读再跑一遍贪心聚类——按重复区比对长度排序，以**编辑距离阈值 T**（初始 10）归簇或开新簇；编辑距离用 Needleman-Wunsch（gap=1，match=0，mismatch=1），带逐行下界剪枝加速；
- **POA 精修循环**：对每簇用**部分顺序比对（partial order alignment, Lee 2003）**取共识更新质心 → 重排序 → 合并质心编辑距离 < T 的簇 → 迭代至收敛；保留读数 > min(10, 0.1×ns) 的簇；若剩余读 < 0.8×ns，则**逐步放宽 T（20→50→80→…→700）**重来。

### 2. HMM 重比对：Dindel 模型 + 几何分布 homopolymer 错误模型
这是 LongTR 最核心的方法创新：
- 用 **Dindel（Albers et al. 2011）式 HMM** 把每条读比对到每个候选单倍型上——HMM 天然允许任意位置、任意大小的 indel，摆脱了 HipSTR「整数倍、单次事件」的限制；
- **错误参数从数据中估计**：ρ = #InDels/Σ|s|（s 为读长与组装真值基因型长度的差）、i = 插入概率、d = 缺失概率，**按 homopolymer 长度区间分别建模**（错误大小服从几何分布）；HiFi 模型用 HipSTR 参考集的 840,248 个 homopolymer、对照 HG002 组装训练；
- 效果：错误建模使 trio 一致率 83% vs 81%，组装一致率 +13%。

### 3. 基因分型：全基因型空间打分 + 归一化质量
- 对每个单倍型对 (Hᵢ, Hⱼ)（含纯合组合）计算得分 `S = Π_r [P(h_r=1)·S_r(Hᵢ) + P(h_r=2)·S_r(Hⱼ)]`，S_r 为 HMM 比对分；
- 对所有可能基因型的得分总和做归一化，输出质量分 Q；
- **定相支持**：沿用 HipSTR 的 phased-read 机制（haplotag 把 P(h_r=k) 钉成 0/1），条件是每个单倍型至少有一条读、未定相读 ≤20%。

### 4. 工程设计
- C++ + HTSlib，**直接从 URL/云端流式读 BAM**，无需本地下载——为 gnomAD 级别 cohort 设计；
- 输出 VCF 含推断的**完整等位基因序列**与 ALLREADS 字段（读级长度观察，也是其错误模型的自标定原料）。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | 排序索引的 BAM/CRAM（带 RG/SM 标签）、参考 FASTA、TR 目录 BED（4-5 列：chrom/start/end/motif[/name]） |
| 关键参数 | `--min-reads 4`、`--indel-flank-len 25`、`--max-tr-len 10000`、`--skip-assembly`、`--phased-bam`、`--haploid-chrs` |
| 输出 | bgzip+tabix VCF：GT、Q、PQ、DP、GL、PL、ALLREADS、START/END/MOTIF、DFLANKINDEL 等 |
| 平台 | HiFi 与 ONT Duplex（README 明确支持） |

## 性能与局限

**性能**（论文自报，30× HG002 HiFi，Adotto 937k 位点）：
- vs HG002 组装：**98.5% 长度一致（±1bp）/ 84.7% 序列一致**；trio Mendelian **86%**（TRGT 同基准 79%）；homopolymer 位点 78.3%；
- **>500 bp 重复上优势最大**；独特检出 **514 个 null allele**（整段重复缺失）；TRGT 漏掉的是截断读/低支持大插入的情形；
- 速度：428 min ≈ TRGT 436 min（单核 4GB）；比 adVNTR 快 15× 且一致率更高（96%）；
- 跨平台：HiFi vs ONT Duplex 基因型一致率 90%；组装一致率 HiFi 99% vs ONT 88%；ONT 能解决一些 HiFi 难比对区的大扩增（Fisher p=1.75e−292）；
- 临床：HTT（n=4）、FMR1（n=3）拷贝数全对。

**局限**：
- 序列一致率（84.7%）绝对值不高——motif 复杂位点仍是难点；
- 目录驱动，无 de novo 发现；
- STRkit 基准中 HiFi F1 0.9555，略低于 STRkit/TRGT；
- 单核实现，caching 大 cohort 时并行度受限。

## 对 CCS-STR 优化的启示

- **错误模型可以自标定**：LongTR 用「读 vs 组装真值」的长度差直接拟合 homopolymer 错误参数——你在 CCS 侧优化错误模型时，这是现成的评估框架（ALLREADS 字段就是原料）；
- 「POA 共识 → HMM 重比对」两段式与 TRGT 的「聚类 → 多数投票共识」是同一思想的不同实现，可对照测试哪种共识对复杂 motif 更稳；
- 免费拿到的 **null allele 检测**是 TRGT 没有的能力，值得在流程中作为互补通道；
- 单倍型聚类阈值 T 的自适应放宽（20→700）策略，可用于你自己的任何读聚类场景。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC11229021
- 仓库：https://github.com/gymrek-lab/longtr
- 前身：HipSTR（见 P2 报告）
- 配套分析：TRTools https://trtools.readthedocs.io 、EnsembleTR https://github.com/gymrek-lab/EnsembleTR
