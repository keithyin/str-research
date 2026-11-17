# trgt-instability

> **优先级**：P1（重要）｜ **类别**：体细胞不稳定量化 ｜ **维护方**：PacBio ｜ **状态**：预发布 v0.5.0（2025，活跃开发中，官方警告仅供实验）
> **一句话**：站在 TRGT 输出之上，用「读 vs 等位基因共识的编辑距离谱 + Dirichlet-multinomial 对照基线」把**体细胞 TR 不稳定**变成一个有 p 值的统计检验。

## 是什么

PacBio 2025 年新发布的工具（二进制分发，2 commits，极早期），消费 TRGT 的 spanning BAM + VCF，对每个等位基因输出不稳定显著性。研究用途，非诊断。

- 仓库：https://github.com/PacificBiosciences/trgt-instability （7★，v0.5.0 pre-release）
- 输入：重复目录 BED + 对照样本的 TRGT 输出 + 病例样本的 TRGT 输出。

## 解决了什么问题

1. **mosaic/体细胞突变没有现成的量化框架**：TRGT 给出等位基因序列与读支撑，但「这个样本的重复比对照不稳定多少、显著吗」没有答案——体细胞重复扩增疾病（如 CANVAS、部分 C9orf72）和癌症里的 TR 不稳定都需要。
2. **没有对照基线的离群检测不可信**：不同位点天然 divergence 差异巨大（长重复 vs 短重复、纯 vs 复杂 motif），直接比编辑距离没有意义。

## 如何解决的

### 1. 读级 divergence 提取
- 对每条读与其所属等位基因的 **TRGT consensus**（+50 bp 侧翼上下文）计算编辑距离；
- 定义**读 divergence 率** = 长度归一化编辑距离；
- 对每个等位基因，把所有读的 divergence 率离散到 **15 个 bin**，得到「不稳定谱」（instability profile）——一个 15 维计数向量。

### 2. 位点特异性基线模型
- 用**对照样本**（正常、非不稳定）的同一等位基因拟合 **Dirichlet-multinomial 模型**——建模"正常情况下这些 bin 的计数如何波动"；
- Dirichlet-multinomial 比 multinomial 多一层"过度离散"参数，能容纳读间真实的生物学差异，避免把正常波动当成显著。

### 3. 检验与输出
- 病例等位基因的 profile 对基线做参数检验：p 值小 = 该 profile 在基线下 unlikely；
- `--report-effect-size`：输出深度感知的效应量 **d_median（读分布与基线的 Wasserstein 距离中位数）+ 95% credible interval**——p 值负责显著性，d 值负责排序；
- **自适应提前终止** `--adaptive-fdr-q`：参数 bootstrap（默认 `--n-sim 100000`）中若某等位基因已不可能达到 BH 校正显著，提前停算——大规模位点扫描的算力控制。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | 重复目录 BED；对照组 TRGT 输出（spanning BAM + VCF）；病例组 TRGT 输出 |
| 中间输出 | divergence 记录（TSV，可 gzip）、每重复的不稳定模型 |
| 最终输出 | 每等位基因：p 值（+可选 d_median、d_ci） |
| 依赖 | 完整 TRGT 前置流程（含 spanning BAM） |

## 性能与局限

- 预发布工具（v0.5.0），接口与统计口径可能变动；
- 需要**成规模的对照队列**来拟合基线（每重复一个模型）；
- 只能检测**相对基线的不稳定**，不解决 de novo 突变检测（那是 TRGT-denovo 的任务）；
- divergence 谱对上游读质量敏感——ccs/比对层改动会直接改变 bin 分布（这既是风险也是可用的回归信号）。

## 对 CCS-STR 优化的启示

- **「读 vs per-allele consensus 的 divergence」是 HiFi STR 上体细胞信号的标准读法**——HiFi 单读错误低到足以让"读间真实差异"（体细胞）从"错误噪声"里分离出来；这正是 CCS 质量红利兑现的地方；
- Dirichlet-multinomial 基线 + Wasserstein 效应量的组合是可移植的统计模板（任何"病例 vs 对照的读级特征"检验都能用）；
- 如果你优化的是 ccs/DC 层，可以用 divergence 谱的**对照样本稳定性**作为无真值集的质量回归指标。

## 参考

- 仓库（README 即文档）：https://github.com/PacificBiosciences/trgt-instability
- 前置：TRGT（P0 报告）
