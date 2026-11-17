# ccs / pbccs（Circular Consensus Sequencing）

> **优先级**：P0（核心必读）｜ **类别**：读生成（Subreads → HiFi）｜ **维护方**：PacBio ｜ **状态**：活跃，v8.x（默认分支 develop，文档站 ccs.how）
> **一句话**：把同一 SMRTbell 分子的多次聚合酶 pass 合成单分子高精度共识读（HiFi）的工具——你的整条 STR 流水的读都从这里出来。

## 是什么

ccs（内部名 pbccs/unanimity 系）是 PacBio 官方的 HiFi 读生成器：输入一个 movie 的 `.subreads.bam`，输出高精度共识读（unaligned BAM / bgzip FASTQ），带预测质量值（QV）。BSD-3-Clause-Clear 许可。

- 仓库：https://github.com/PacificBiosciences/ccs （~127★）｜ 官方文档：**https://ccs.how**
- 论文：*Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome*, **Nature Biotechnology 2019**（Wenger & Peluso 等）https://www.nature.com/articles/s41587-019-0217-9
- 运行位置：Sequel IIe / Revio **上机自动执行**（Revio 上中间 CLR BAM 不再落盘，因为单 SMRT Cell 25M 超 1 TB）。

## 解决了什么问题

1. **单 pass 读精度不够**：PacBio 单 pass 读精度 ~85–90%，直接用于 STR 分型不可行（STR 区错误会直接变成假拷贝数变化）。
2. **错误谱偏科**：即使做成了 HiFi，**残余不一致里 92.0% 是 homopolymer indel**（每 477 bp 一个），而 mismatch 只有每 13,048 bp 一个——错误恰好集中在 STR 最敏感的错误类型上。这就是为什么 STR 分型器都要 homopolymer-aware 错误模型。
3. **QV 可信度**：下游全部按预测 QV 过滤，需要 QV 与实测一致率严格对齐。

## 如何解决的

### 1. 共识生成（统计模型 + Arrow）
- 一个 ZMW 的多条 full-pass subreads 被用于构建共识：`ccs` 早期版本走「选种子读 → 其余 subreads 比对到种子 → **Arrow（HMM 共识模型）**逐窗口产出共识碱基与 QV」的路线；当前版本在此基础上叠加 **DeepConsensus**（见下）；
- **`np` 标签**：计入共识的 full-length subread 数（须通过长度过滤与 polishing）；v3.3.0 起按窗口取最小 np，v4.0.0 起改为取 mode；`ec` 标签是含部分 subread 的有效覆盖（约 np+1）。

### 2. DeepConsensus 集成（Revio，SMRT Link v12.0 起）
官方文档（ccs.how/faq/revio）描述的三段式流程，**这就是当前 Revio HiFi 读的真实生成路径**：
1. **Arrow** 生成初始共识与 QV；
2. 初始共识中**低质量窗口**交给 DeepConsensus（gap-aware transformer）精修——典型人类 WGS run 中 **30–70% 的窗口**被 DC 润色；
3. 再用 Arrow 对润色后序列**重算全读 QV**（下游工具期望单一模型的 QV；Arrow 在 DC 窗口会略微低估，故最终值不计入 rq）。
- 实现：C++ 前端 + ONNXruntime，跑在数据中心 GPU 的 TensorRT core 上；Arrow 润色同样 GPU 化（比双 64 核 EPYC 快 10×）。
- QV 是「DC 窗口 QV + Arrow 窗口 QV 的平均」，碱基质量分**7 个 bin** 压缩文件体积。

### 3. 质量控制语义
- `--min-passes`：最少 full-pass 数；文档记录了 v4.2.0 的一个关键修复——低 min-passes（如 1）时若只有一条 subread 比对上草稿就可能放过差草稿，v4.2.0 起额外要求 **>50% subreads 比对到草稿**；
- **`--top-passes` 默认 60**：按中位长度排序最多取 60 条 full-pass，超过无精度收益；
- **Q20 阈值语义**：预测 QV<20 的读默认丢弃（进 `fail_reads.bam`）——注意 DeepConsensus 救不回被丢弃的读，所以官方建议配合 DC 时调整 ccs 参数；
- `--by-strand`：按正反链分别要求证据（对特殊分子结构有用）。

### 4. QV 校准
论文实测：预测 QV 与实测一致率平均偏差 **−1.2（Q 单位）**，即 QV 基本诚实；**约 10 个 full-pass 才能达到 Q30 中位数**——这是 min-passes/模板长度权衡的依据。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | 单 movie 的 subreads BAM |
| 输出 | HiFi BAM（unaligned）/ FASTQ / SMRT Link XML |
| 关键参数 | `--min-passes`、`--min-rq`（QV 阈值）、`--by-strand`、`--top-passes`、`--max-drops` |
| 性能注意 | `--chunk` 并行 + 环境变量调优（文档站 FAQ） |

## 性能与局限

**性能**（2019 论文，ccs v3.0.0，`--minPasses 3 --minPredictedAccuracy 0.99`）：
- HG002 ~28×：读长 13.5±1.2 kb，中位 Q30 / 均值 Q27（99.8%）；
- 残余不一致：92.0% homopolymer indel、4.6% 非 homopolymer indel、3.4% mismatch；vs NovaSeq：mismatch 低 17×、**indel 高 181×**；
- 读级伪影 ~0.6%（0.5% 嵌合、0.1% 低质量段）；
- 下游：DeepVariant SNV 99.91%/99.96%；组装一致率 Q44–48；99.64% 变异定相（block N50 206 kb）；
- 成本：平均 3,035 CPU 核时/SMRT Cell。

**局限**：
- homopolymer/低复杂度区仍是主要残余错误源（对 STR 是核心痛点）；
- QV 是「预测」——极端区（G4、超长 homopolymer）可能失准（TREAT 论文的覆盖度丢失证据）；
- 上机版中间产物（CLR）不可回放，调参实验需要用离线 ccs 重跑 subreads。

## 对 CCS-STR 优化的启示

- **你的优化空间直接映射到这份错误谱**：STR 分型的 off-by-one 错误 ≈ homopolymer indel；ccs/DC 的任何 homopolymer 改进都应直接体现为 STR 分型孟德尔一致率提升——用 GIAB TR 真值集量化；
- **np/ec/top-passes 是现成的读级特征**：STR 位点上的 np 分布可以成为「该位点共识可信度」的协变量（类似 trgt-instability 用读 divergence）；
- `--by-strand` 与 >50% subread 支持规则说明 ccs 内部有「证据一致性」概念——STR 优化若要下探到读生成层，这些语义是接口；
- QV 7-bin 压缩意味着**单读碱基质量分辨率有限**，自己做 per-base 错误建模时要注意。

## 参考

- 论文：https://www.nature.com/articles/s41587-019-0217-9 （PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6776680）
- 官方文档：https://ccs.how （FAQ：accuracy-vs-passes、revio、fail reads、QV binning、low complexity）
- DeepConsensus 集成细节：https://ccs.how/faq/revio
