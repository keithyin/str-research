# STRdust

> **优先级**：P1（重要）｜ **类别**：TR 分型（极速/临床位点点火）｜ **维护方**：Wouter De Coster / Iñes Höijer / Adam Ameur（Uppsala 大学）｜ **状态**：活跃，v0.16.x（2026-11）
> **一句话**：Rust 写的极速 TR 分型器——默认从现有比对的 CIGAR 直接"抠"出等位基因，sensitive 模式把参考中的重复**挖掉再重比对**，30× 基因组 2,000 位点只要 7 秒。

## 是什么

STRdust 是 ONT 圈著名工具作者 De Coster（NanoPlot/Mosdepth 作者）开发的 TR 分型器，Rust 实现，MIT 许可，BAM/CRAM 进、VCF 出。起源于 GOLGA8A 重复扩增研究，后泛化。

- 仓库：https://github.com/wdecoster/STRdust （~28★）
- 论文：*STRdust: tandem repeat genotyping from long reads*, **Bioinformatics 2024**（ISMB/ECCB）https://doi.org/10.1093/bioinformatics/btae268
- 姊妹资源：pathSTR（医学相关 TR 的正常人群等位基因分布，Genome Research 2024）。

## 解决了什么问题

1. **速度**：全基因组扫描级应用（大 cohort、pipeline 中间步骤）需要分钟级而非小时级的 TR 分型。
2. **比对器 clip 穿过扩增**：读比对在重复边界被剪断时，CIGAR 里根本看不到扩增——需要一条"重新比对"的补救通道。
3. **临床位点的即开即用**：不想维护目录时，希望一行命令对 STRchive 全部致病位点直接分型。

## 如何解决的

### 1. fast 模式（默认）：CIGAR walk
直接沿现有比对的 CIGAR 走读段，把落在重复区间内的部分抽出来作为等位基因序列。前提是比对本身把重复区连续比对上了——快，但对比对质量敏感。

### 2. sensitive 模式：切除重复的人工参考重比对
把参考中该重复序列**挖掉**生成人工参考，将读用 minimap2 重新比对到这个"无重复参考"。重复内部的比对歧义消失，clip 穿过扩增的读也能连续比对——然后重复部分自然暴露为大的比对差异。慢一个量级，但对扩张读更鲁棒。

### 3. 单倍型聚类（与 TRGT/LongTR 同思路的不同实现）
- 有外部定相时直接用；否则 **Ward 层次聚类**（长度加权 Levenshtein 距离）与 **DBSCAN**（k-mer 组成空间）聚类，两者可互为 QC 模式；
- 输出 ≤2 个单倍型 + 每簇共识序列（写进 VCF 的 ALT）。

### 4. 工程加速
- **QUICKREF**：先用浅检查跳过明显纯合参考位点，省掉大多数位点的完整计算；
- `--priority` 模式（sensitive/balanced/precise/expanded）显式在召回与精度间取舍。

### 5. QC 语义（值得借鉴）
VCF 里带显式 QC 标签：`UNSPANNED_EXPANSION`（扩增无跨读支持）、`EXPANSION_OUTLIER`、`IMPRECISE_LENGTH` 等——把"这个位点我不确定"编码进数据而非日志。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | 参考 FASTA + BAM/CRAM（本地或 URL）；位点：`-r` 单区间 / `-R` BED / `--pathogenic`（直接对 STRchive GRCh38 致病 STR） |
| 输出 | VCF：共识序列 ALT、RB/FRB/MRL 长度字段、MOT（motif）、LPS/LPSM（最长纯重复段）、QC 标签 |
| 注意 | 坐标是目录区间而非最小变异表示——多样本合并时要做归一化 |

## 性能与局限

- **速度**：30× ONT 基因组、2,000 个目录位点 **7 CPU 秒**（同场景 TRGT/LongTR 是分钟–小时级）；
- STRkit 论文基准：HiFi F1 0.8327 / ONT 0.7581——**精度明显低于 TRGT/LongTR/STRkit**，速度换精度；
- 论文作者报告其在致病等位基因检出上较灵敏（de Coster 用它处理 GOLGA8A 等）；
- 无甲基化、无 de novo 发现、超读长等位基因不能分型。

## 对 CCS-STR 优化的启示

- **"挖掉重复再比对"是一个便宜可复用的消歧技巧**——任何重复区间内的比对歧义问题都可以用这招做个对照通道；
- QC 标签进 VCF 的做法值得你在自己的输出格式里采用；
- 用作**大 cohort 初筛 + 异常位点上报**的角色最合适：先 STRdust 扫全基因组，异常位点再用 TRGT/LongTR/otter 精析——一个分层分型架构。

## 参考

- 论文：https://doi.org/10.1093/bioinformatics/btae268
- 仓库：https://github.com/wdecoster/STRdust
- pathSTR：https://pmc.ncbi.nlm.nih.gov/articles/PMC11610575
