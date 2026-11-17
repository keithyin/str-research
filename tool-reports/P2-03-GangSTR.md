# GangSTR

> **优先级**：P2（背景借鉴）｜ **类别**：短读全基因组 TR 分型（读类别似然）｜ **维护方**：Gymrek lab（UCSD，Mousavi 等）｜ **状态**：维护模式
> **一句话**：把读分成**四类**（enclosing/spanning/FRR/flanking）各建统计模型、联合最大似然分型——「读类别」框架被后来所有长读工具继承。

## 是什么

GangSTR 从短读 WGS/WES 对全基因组 TR（含扩增）做分型，C++。

- 仓库：https://github.com/gymreklab/gangstr
- 论文：*Profiling the genome-wide landscape of tandem repeat expansions*, **Nucleic Acids Research 2019**;47(15):e90 https://pmc.ncbi.nlm.nih.gov/articles/PMC6735967
- 配套：TRTools 下游生态、EnsembleTR、strdenovotools/MonSTR。

## 解决了什么问题

1. **超过读长与片段长的扩增无法用跨读分型**：lobSTR/HipSTR 只用 enclosing/spanning 读，对大扩增直接失效。
2. **不同证据要不同统计模型**：read-pair 的片段长度收缩、in-repeat 读的锚距离、off-target 救回——混在一个模型里比单通道强。

## 如何解决的

### 四类读 × 四个统计量（框架核心）
| 读类 | 统计量 | 模型 |
|---|---|---|
| **Enclosing（E）**：整段重复 + 两侧 ≥20bp 侧翼在一条读内 | 拷贝数 n | 直接计数 + PCR stutter（几何分布，P=0.9，沿用 HipSTR） |
| **Spanning（S）**：片段整体跨过重复，两 mate 各锚一侧 | 表观片段长 Δ | 正态分布（均值随重复长度变化）——扩增使 Δ 收缩 |
| **FRR（完全重复读）**：至少一条 mate 全由重复 motif 构成 | 非重复 mate 的锚距 Ω + FRR 总数 | Ω 越小重复越长；FRR 计数 ~ Poisson（λ 线性于等位基因总长 × 覆盖 × 读长 × motif 长度） |
| **Flanking（F）**：读一端伸入重复 | 可见拷贝数 k | k ~ Uniform(1, A)，给出"某等位基因至少长 k"的下界证据 |

### 联合似然
- 基因型 ⟨A,B⟩ 的对数似然 = Σ 读对 [P(C(r);A) × P(r|C(r);A)] + Poisson(FRR 计数项)；
- 每个读对等可能来自两条等位基因；跨类读对同时贡献各类；
- 类概率 P(C(r);A) 有闭式解（Supplementary）；优化默认网格搜索、大范围用 COBYLA 多步法；bootstrap 给 CI 与扩展概率 P0/P1/P2。

### 超读长扩增
- 无 enclosing 读时靠 **Δ 收缩 + Ω 锚距 + FRR 计数**三路证据；还支持扫描预定义 off-target 区域救回"锚在别处"的 FRR 读；
- 宣称 "not limited by fragment or read length"（统计意义上有精度上限）。

## 性能

- 模拟 14 个致病 TR（40×）：各位点 RMSE 最低；ExpansionHunter 在双等位基因都长时不稳，Tredparse 超片段长即低估；
- 实测 HTT RMSE 7.9（Tredparse 8.3 / EH 10.1）；WES 200 HD 患者 RMSE 5.4（EH 9.1 / Tredparse 96.4）；
- 基因组级：CEU trio ~673k TR，Mendelian 98.9%（ML 口径）；发现 43,571 个 HipSTR 漏掉的 TR；NA12878 的 56 个候选扩增 53 个有 PacBio 支持、9/11 CE 确认；
- 速度：全基因组 ~16.6 CPU 时（vs STRetch ~157）。

## 局限

- 短读平台统计模型的分辨率上限（GC-rich FMR1 所有工具都差）；后续基准中长读工具全面超越（TRGT 论文里 GangSTR Mendelian 仅 72.90%）；
- 扩增很长的位点用 FRR 计数反推，误差大、且受覆盖建模影响。

## 对 CCS-STR 优化的启示

- **「按证据来源分通道建模」的框架永不过时**：HiFi 上你的 in-repeat 读（超长扩增）同样需要 Poisson/锚距类模型补位——TRGT 只用 spanning 读，GangSTR 思路是天然扩展方向；
- 四类读的判定规则（尤其 FRR 的锚距分布）可直接移植为 HiFi 读分类器；
- 扩展概率 P0/P1/P2 的输出设计对临床报告友好（"多大置信度认为扩展了"）。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC6735967
- 仓库：https://github.com/gymreklab/gangstr
- 下游：https://trtools.readthedocs.io
