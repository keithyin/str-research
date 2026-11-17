# tandem-genotypes（tdg）

> **优先级**：P1（重要）｜ **类别**：TR 扩增扫描（免目录 de novo）｜ **维护方**：Martin C. Frith lab（东京大学；论文一作 Mitsuhashi）｜ **状态**：活跃，~v1.5+（master，bioconda 分发）
> **一句话**：用 LAST 概率比对的 mismap 过滤 + 「按重复单元整数倍取整」的读长度变化估计，把全基因组 ~70 万个重复按"重要性"排序找扩增——免目录、轻量、快。

## 是什么

tandem-genotypes（命令名 `tandem-genotypes`，俗称 tdg）从**已比对的长读**（MAF 格式，通常 LAST/last-split 产生）中找串联重复的长度变化，Python 实现（画图需 R），GPL-3.0。常被误传为 Eichler lab 作品——实为 Frith lab。

- 仓库：https://github.com/mcfrith/tandem-genotypes （~52★）
- 论文：*Tandem-genotypes: robust detection of tandem repeat expansions from long DNA reads*, **Genome Biology 2019**（Mitsuhashi … Frith）https://pmc.ncbi.nlm.nih.gov/articles/PMC6425644
- 配套：`tandem-genotypes-join`（病例-对照合并重排）、`-plot`（直方图）、`-merge`（等位基因共识，需 lamassemble）。

## 解决了什么问题

1. **免目录发现**：不想依赖致病位点目录或全基因组注释目录，希望直接对 UCSC simpleRepeat/RepeatMasker 输出的 ~70 万重复做扫描。
2. **error-prone 读的比对歧义**：重复区读比对错位会把"比对错误"误判为"重复变化"——需要显式的歧义概率控制。
3. **批量病例的优先级排序**：单患者扫描出上百个候选，需要能区分"真扩增"与"比对伪影"的排序机制。

## 如何解决的

### 1. 比对层：LAST + last-split + mismap 过滤
- `lastdb -uNEAR` 建索引、**last-train** 从数据样本学习插入/缺失/替换错误率、`lastal | last-split` 产出**概率分裂比对**：每段比对带 **mismap probability**（错位概率）；
- 丢弃 mismap > 10⁻⁶ 的比对；
- **tantan 软屏蔽**：比对中"简单序列"（如 atatat…）由 tantan 标成小写，主要小写的比对直接丢弃——防止把重复内部的比对漂移当成重复变化；可选项 WindowMasker 加速；
- 同染色体同链的共线比对间隔 ≤10⁶ bp 时连接。

### 2. 长度变化估计：整数倍取整规则（核心 ad hoc 逻辑）
定义两个距离：far `f = max(100, u)`、near `n = max(60, u)`（u = 重复单元长度）：
- **单比对情形**：比对必须**两侧都超出重复 ≥ f**（侧翼锚定），否则该读不用于该重复；检查比对 gap：
  - 不与重复重叠且 i ≤ u/2、或整体距离 ≥ f、或整体 > n：忽略；部分 ≥ f：放弃该读；
  - r = 重复内部未比对的参考碱基数；净缺失 D = min(d−i, r)；gap 贡献 = **D 四舍五入到最近的 u 的整数倍**（平局归零）；
- **多比对情形**：左右比对必须各自越过重复 ≥ f；插入大小 = 左右比对之间的 read 碱基数 − reference 碱基数，同样取整到 u 的整数倍。

**直觉**：重复区的 indel 错误大多数不是单元长度的整数倍；取整规则天然滤掉非整数倍的比对噪声，只保留"真的少了/多了若干个单元"的信号。

### 3. 优先级排序
- 每个重复丢弃最极端的预测（若 ≥3 个读）；
- 分数 = 长度变化碱基数 ÷（参考重复长度 + 30）× 基因注释权重（coding 50 / UTR 20 / 启动子 15 / ncRNA 外显子 15 / 内含子 5），聚谷氨酰胺/聚丙氨酸编码区 ×2；
- **病例-对照模式（-join）**：病例分數 d = 患者分数立方均值，健康 h = 对照立方均值，联合分 = max(d − max(h,0), 0)——对照可用 NA12878 PacBio / NA19240 ONT / 黑猩猩 panTro5。

### 4. 双等位基因粗估（-o2）与链分离
- `-o2` 给出粗糙的双等位基因预测；
- 输出按正反链分离——防止把链偏倚误判成两个等位基因（作者在同行评议中感谢 Onoguchi 指出这一点）。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 比对输入 | MAF 格式（LAST）；FASTA/FASTQ（画图/合并用）；gz 支持 |
| 重复输入 | UCSC simpleRepeat.txt / microsat.txt / rmsk.txt、RepeatMasker .out、tantan 输出、4 列 BED；自带 hg19/hg38 疾病重复文件 |
| 输出 | 文本表（BED3 + 单元 + 基因 + 基因部分 + 每读每链拷贝数变化 + 可选双等位预测）；直方图 PDF |

## 性能与局限

**性能**（论文）：
- 质粒混合读：CAG/CAA/GGGGCC/CCTG 拷贝数变化全对（GGGGGCC 有链偏倚双峰）；
- **BAFME 患者 ~6× PacBio**：3 条读支撑的 ~5 kb SAMD12 扩增在全基因组 0.7M 重复中排第 4（无对照第 13）；RepeatHMM/PacmonSTR/NGMLR-sniffles/nanoSV/PBSV 全部失败；
- 排序基准：7 个编码区扩增经对照合并后**全部进入全基因组 top 10**；
- 速度：LAST 13 h（vs NGMLR 37 h）+ tdg ~1 h；
- 同行评议期间顺手发现了 NIID 的 NOTCH2NLC GGC 扩增。

**局限**：
- 长度变化口径，无序列/motif 组成解析；STRkit 论文基准中它在 HiFi 上落后于主流工具；
- 输出是"每读变化表"而非标准 VCF 基因型（-o2 的双等位预测较粗）；
- 需要 LAST 生态（MAF），与主流 BAM 流程衔接需要转换。

## 对 CCS-STR 优化的启示

- **「整数倍取整」是把错误模型变成过滤器的最便宜实现**——HiFi 错误率低，这一招在 HiFi 上应该更干净；你的任何读级拷贝数估计模块都可以加这层;
- mismap-probability 过滤与 pbmm2（无 secondary 输出）体系不兼容——若想引入 tdg 思路，需要 LAST 独立通道或用 minimap2 的 mapq 代理；
- 「病例-对照立方均值差」的排序统计，是 cohort 层面找离群位点的轻量方案（对比 EHdn 的 Wilcoxon/z-score）。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC6425644
- 仓库：https://github.com/mcfrith/tandem-genotypes （README 极详尽，含疾病位点文件）
