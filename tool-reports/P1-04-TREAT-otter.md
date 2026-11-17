# TREAT + otter

> **优先级**：P1（重要）｜ **类别**：组装式 TR 分型 + 端到端注释/分析套件 ｜ **维护方**：Holstege lab（Radboud UMC / TU Delft，Tesi 等）｜ **状态**：开源（otter: C++；TREAT: Python/R）
> **一句话**：otter 用**靶向局部组装**替代"读-参考比对"来做 TR 分型（对大 TR 比 TRGT 更准）；TREAT 在其上搭了 motif 注释、离群分析、病例-对照的全套分析层。

## 是什么

- **otter**（One Targeted-locus Read-assembler）：C++ 局部组装器，htslib 读 BAM、WFA2-lib 重比对，输入 BAM+FASTA+BED。
- **TREAT**（Tandem REpeat Annotation Toolkit）：以 otter 为核心的端到端 workflow（Python/R），只吃 BAM+参考+BED。
- 仓库：https://github.com/holstegelab/otter 、https://github.com/holstegelab/treat
- 论文：*Characterizing tandem repeat complexities across long-read sequencing platforms with TREAT and otter*, **Genome Research 2024**;34(11):1942–1953 https://pmc.ncbi.nlm.nih.gov/articles/PMC11610572

## 解决了什么问题

1. **参考比对在重复区是伪影之源**：所有"比对→聚类"路线（TRGT/LongTR）都先验地信任比对器；长重复内部比对断裂/歧义直接传导到分型。otter 干脆**在位点内组装**，绕开参考比对。
2. **平台间没有公平基准**：HiFi vs ONT 在不同重复长度上的相对优劣此前只有零散说法——这篇论文给出了 5–30× × 四平台 × 161k/864k 位点的系统对比。
3. **"覆盖度丢失"没人量化过**：为什么某些临床位点 HiFi 也会没数据？——GC-rich/G-四联体重复阻断聚合酶 multiple passes。

## 如何解决的

### otter：靶向局部组装
1. **读选择**：htslib 抽取跨读与非跨读；非跨读若出现 clip 信号，用 **WFA2-lib** 把 ~100 bp 参考侧翼对上去，相似度 ≥90% 则升级为跨读（默认参数）；
2. **局部基线错误率估计**：对跨读做两两比对，用**高斯核密度估计**估计成对距离的局部分布——**单峰在 0 = 纯合；多峰 = 多等位基因**。这是"错误率自适应"的关键：不同位点/平台的噪声水平不同，先估背景再分信号；
3. **层次聚类**：按成对距离层次聚类，停在密度最高峰处；一簇保留 curate 步骤（合并欠支撑簇、或按最大等位基因数 α 截断），**随局部覆盖度自适应**；
4. **共识**：每簇用 pseudo-POA 构建共识序列。

### TREAT：分析层
- **Assembly analysis**：otter 分型 → pytrf 做 motif/拷贝数表征 → 多样本 gVCF（等位基因序列+大小+重复内容）；
- **Reads analysis**：按 TR 大小迭代聚类 + 每读 motif 表征（适合体细胞不稳定、MSA）；
- **分析模块**：离群检测（归一化到中位等位基因大小的距离 + χ² 检验）、病例-对照 logistic 回归（短/长/联合等位基因大小）、MotifScope 可视化。

### 平台基准结论（论文最有引用价值的部分）
- **PacBio vs ONT 按重复长度**：<500 bp **PacBio 更准**；500–1000 bp 打平；**≥1000 bp ONT 更好**；覆盖度永远越高越好；Duplex > Simplex；
- otter vs TRGT（PacBio）：总体持平，但 **>500 bp 大 TR otter 更准**（错误率 0.2–2.5% vs TRGT 0.6–3.8%）；
- motif 一致性：TREAT vs TRGT ~96%；
- **覆盖度丢失**：864k TR 中 486 个（0.06%）HiFi 覆盖比 <0.25——更长（平均 +493 bp）、GC 更富（+17.4%）、AG/CC/CG/CT/GG 双核苷酸富集、**G-四联体富集（+3.76%）**；437/486 为 PacBio 特有（OR=9.4）；
- **ABCA7 补救案例**：HiFi-only（局部覆盖 1–7×）与 Southern blot 相关性 0.58；**加入 non-HiFi 读（覆盖升至 ~22×）后相关性 0.99**；RFC1 同样被救回。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| otter 输入 | BAM + 参考 FASTA + BED |
| TREAT 输入 | 同上 + 可选注释 |
| 输出 | gVCF（等位基因序列）、离群/回归结果表、可视化 |

## 性能与局限

- 5×–30× 剂量曲线完整给出；高覆盖永远更好，Simplex ONT 受错误率限制；
- 组装路线天然比比对路线慢（有针对性位点数的场景最合适）；
- 局部组装对覆盖度的要求更高——超低覆盖位点更适合比对式工具。

## 对 CCS-STR 优化的启示

- **「先估局部噪声、再分信号」**是组装路线的精髓——你做任何聚类/分型模块时，先把"这个位点的真实错误水平"估出来，比全局统一阈值更稳；
- **G4/GC-rich 位点是 HiFi 的系统性盲区**（0.06% 但集中在临床相关位点）：管线里应有一层「覆盖度异常检测」，对这类位点自动降置信度或提示补 non-HiFi 数据；
- 「非跨读 + 侧翼 WFA 重比对升级为跨读」与 STRkit 的 `--realign`、Straglr 的 clip 补救同源——这是社区已收敛的廉价技巧；
- 论文的平台×长度结论可直接写进你选择工具组合的依据（如：HiFi 为主 + ONT 补长重复）。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC11610572 （Genome Research 34(11):1942）
- 仓库：https://github.com/holstegelab/otter 、https://github.com/holstegelab/treat
