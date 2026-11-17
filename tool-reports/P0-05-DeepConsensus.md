# DeepConsensus

> **优先级**：P0（核心必读）｜ **类别**：读生成（ML 后纠错）｜ **维护方**：Google（与 PacBio 合作）｜ **状态**：r1.2 / 1.2.0（pip），Revio 上已集成进 ccs
> **一句话**：在 ccs 之上叠加一层 **gap-aware transformer**，专门修正 HiFi 共识的残余错误（-42%），把高 Q 读的产出提高一个量级——STR/homopolymer 区正是它的主战场。

## 是什么

DeepConsensus 是 Google 出品的 HiFi 读后纠错模型：输入 subreads（与 ccs 相同的原料），输出**修正后的 FASTQ 读**，替代 ccs 原始输出进入下游。Python/TensorFlow 实现，BSD-3-Clause。

- 仓库：https://github.com/google/deepconsensus （~265★，v1.2.0）
- 论文：Baid, Cook, Shafin et al., *DeepConsensus improves the accuracy of sequences with a gap-aware sequence transformer*, **Nature Biotechnology 2023**;41(2):232–238. DOI: 10.1038/s41587-022-01435-7
  - ⚠️ 注意：社区/部分资料常误写为 Nature Methods——已核实为 **Nature Biotechnology**（s41587 DOI 前缀即 NB）。
- 运行形态：独立流程（`pip install deepconsensus[gpu]==1.2.0`）或 **Revio 上机集成**（ccs 内部 C++/ONNXruntime/TensorRT，见 ccs 报告）。

## 解决了什么问题

1. **Arrow（HMM）的残余错误**：HiFi 平均 Q27–30 意味着每条 15 kb 读仍有 ~10–50 个错误；对 STR 分型，单个 homopolymer indel 就足以造成 off-by-one。
2. **产出损失**：为了让共识更准而过滤低质量分子，会损失数据量；错误修正比丢弃更划算。
3. **柱状编码的盲区**：传统共识/纠错模型把序列编码成对齐后的列（column encoding），**gap（indel）位置没有自然表示**——而 HiFi 的主要残余错误恰恰是 indel（92% homopolymer indel）。

## 如何解决的

### 1. Gap-aware 输入编码（GATE 的"gap"从这来）
- 把窗口内的 subreads **比对到 ccs 草稿序列**，逐列编码；
- 特征包括：每列的碱基、碱基质量、正反向、以及 **gap 的显式编码**——插入用"插入伙伴碱基（insertion partner）"表示，删除用占位符表示；
- 这使 transformer 能直接"看到" indel 证据，而不是把它们塞进 mismatch 通道。

### 2. 模型与损失（论文核心）
- 架构：**GATE（Gap-Aware Transformer Encoder）**，对每个输出位置做多分类头预测「match / mismatch→哪一种碱基 / insertion→哪种碱基 / deletion」；
- 训练用 **alignment-based loss**：以 HG002 真值组装为参考，把模型预测与真值的编辑操作对齐后计损失，而不是逐位置交叉熵——gap-aware 的另一半；
- 2026 年更新（PacBlog）：**banded AlignmentLoss** 经 **AlphaEvolve**（演化搜索）调优后，实测达 Q30 的读比例从 **47.9% → 53.2%**，吞吐 ~18 Gb/hour；此前版本已针对端粒/重复区优化。

### 3. 集成方式（与 ccs 的分工）
- **独立跑**：ccs 全量输出后，DeepConsensus 用 subreads 重算并输出修正 FASTQ；
- **Revio 上机**：见 ccs 报告——只有 **30–70% 的低质量窗口**被 DC 润色，再统一用 Arrow 重算 QV；
- **关键交互**：ccs 默认丢掉预测 QV<20 的读，**DC 无法救回已被丢弃的分子**——官方 quick start 给出了配套的 ccs 参数（放宽过滤，让 DC 接手）。

### 4. 版本配伍
- 下游配 DeepVariant 时：官方建议 **DeepVariant v1.5 + 正常 pacbio 模型**（不要用 DC v0.1 训练的旧模型）；DC v1.2 输出将用于训练下一代 DV 模型；
- 组装：v1.2 优于 v0.3（contiguity/accuracy 双升）。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | subreads BAM（与 ccs 同源）+ ccs 输出 |
| 输出 | 修正后 FASTQ（替代 HiFi 读进入下游） |
| 算力 | 每 SMRT Cell ≥500 shards，每 shard ≤16 CPU（更多无增益）；GPU 版走 `[gpu]` extra |
| 配套 | `calibration.md`（分数校准）、`yield_metrics.md`（3 种化学/读长分布的实测产量） |

## 性能与局限

**性能**（论文 + 官方 metrics）：
- 相比 pbccs 单独：**残余错误减少 42%**；
- 固定测序量下的读产出：**Q20 +9%、Q30 +27%、Q40 +90%**——越高质量的读增益越大；
- 下游：DeepVariant 小变异更准；组装 contiguity/accuracy 提升；SNP/indel 谱更干净。

**局限**：
- 模型是「预测器」，不是统计共识——极端 motif/超长 homopolymer 的行为依赖训练分布（CANVAS 类 G4 位点的覆盖度丢失仍存在，TREAT 论文实证）；
- QV 语义在 Revio 集成模式下的 rq 与碱基质量均值可能轻微不一致（ccs 文档自述）；
- 独立跑需要显著算力（500 shards × 16 CPU / SMRT Cell）；
- "非 Google 官方产品、不用于临床"（README 自述）。

## 对 CCS-STR 优化的启示

- **你的 STR 优化如果上探到读生成层，DC 是最大的已知杠杆**：-42% 错误中相当部分是 homopolymer indel，直接对应 STR 分型的 off-by-one；
- **评估闭环现成**：用「ccs-only vs ccs+DC」两版读分别跑 TRGT/LongTR，在 GIAB TR 真值集上对比孟德尔一致率——这是最干净的 A/B 实验设计；
- DC 只润色低质量窗口（30–70%）的策略说明「哪里错修哪里」——若你发现 STR 位点系统性失败，值得检查这些窗口是否落入 DC 的处理范围；
- AlphaEvolve/alignment loss 的思路（把纠错目标定义为对齐后编辑操作而非逐列交叉熵）对任何序列纠错任务（包括你自训模型）都适用。

## 参考

- 论文：https://www.nature.com/articles/s41587-022-01435-7 （NB 41:232–238）
- 仓库：https://github.com/google/deepconsensus （docs: quick_start / train_model / calibration / yield_metrics）
- Revio 集成：https://ccs.how/faq/revio
- 最新进展博客：https://www.pacb.com/blog/improving-hifi-sequencing-accuracy-with-google-deepconsensus-and-alphaevolve （2026-04）
