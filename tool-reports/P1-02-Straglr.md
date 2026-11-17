# Straglr

> **优先级**：P1（重要）｜ **类别**：TR 扩增发现 + 分型（双模式）｜ **维护方**：BCGSC / Birol lab（Readman Chiu 等）｜ **状态**：活跃，v1.3.0（pip）
> **一句话**：同时支持「全基因组扫描发现新扩增」与「靶向目录分型」的长读工具，GMM 聚类等位基因，被 ONT EPI2ME 官方流程采用。

## 是什么

Straglr 是加拿大 BC Cancer 基因组科学中心（BCGSC，Birol lab）的工具，Python 实现（依赖 TRF 与 blastn 在 PATH 中），输入长读比对 BAM（推荐 minimap2 `-Y`），输出 BED/TSV/VCF。

- 仓库：https://github.com/bcgsc/straglr （~174★，300 commits）
- 论文：*Straglr: discovering and genotyping tandem repeat expansions using whole genome long-read sequences*, **Genome Biology 2021**（Chiu, Rajan-Babu, Friedman, Birol）https://pmc.ncbi.nlm.nih.gov/articles/PMC8361843

## 解决了什么问题

1. **多数工具是"目录驱动"的，发现不了目录外的新扩增**——疾病研究中经常是"不知道哪里扩了"。
2. **噪声读的 motif 识别**：error-prone 长读上，复杂/不完美重复需要宽松的 TR 检测参数（TRF 默认参数不适用）。
3. **比对器在长重复边界丢 split alignment**：minimap2 偶尔不报 ≥1 kb 重复两端的 split 比对，需要补救。

## 如何解决的

### 1. 基因组扫描模式（de novo 发现）
- 从比对 BAM 中找**相对参考的大插入**：单比对内直接读 CIGAR 插入段；或从 **split 比对**重建（断点位置接近、read 位置相距远）；
- 跨读合并邻近插入事件；可提供 bypass 列表跳过着丝粒/片段重复等难点区；
- 每个插入序列用 **TRF（宽松参数）**筛查：要求「≥70% 由单个 TR 构成」才保留；
- 用插入的 motif 与基因组邻区的 motif **松匹配**（VNTR 用降种子长的 BLASTN）确认位点。

### 2. 分型模式（靶向）
- 每位点抽取「夹在两侧邻核苷酸坐标之间」的读序列；
- **clip 读补救**：把参考 80 bp 侧翼与 clip 部分做 BLASTN，找到唯一命中即可重建缺失边界，提取完整重复子序列；
- 提取序列再过 TRF，motif 与目标匹配（先正则，失败降级到允许错配的 BLASTN）；不匹配则不出分型（宁可缺失不硬给）。

### 3. GMM 等位基因判定
- 每读的重复大小（TRF 外推的拷贝数/长度）做样本；
- 用 scikit-learn 拟合 **1..k 个高斯分量的 GMM**（k 默认最大 2），**按最小 AIC 选分量数**——把"几倍体等位基因"变成模型选择问题；
- 每簇取**中位数**作等位基因大小；输出每等位基因的支撑读名清单（mosaic 比例可直接看读数分布）。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | BAM（minimap2 `-Y` 推荐）+ 参考 FASTA；可选 `--loci` BED、`--exclude`、染色体过滤 |
| 关键参数 | `--min-supports 2`、`--min-ins-size 100`、motif 长度范围默认 2–50 bp、最大簇数 |
| 输出 | TSV（每支撑读明细）、BED（位点级等位基因大小/拷贝数/支撑数）、VCF |
| 周边 | `straglr_compare.py`（样本间比较）、`extract_repeats.py`（提取重复序列） |

## 性能与局限

**性能**（论文自报）：
- 模拟（NanoSim 30×，17 个疾病位点）：灵敏度 99%（805/817），特异性完美；tandem-genotypes 98%；RepeatHMM 大等位基因系统性偏小且有假阳性；
- 分辨力：2000 bp 基底上可分辨 100 bp 差异（9/17 位点）；可处理三等位基因位点；
- FXS mosaic 模拟：FM（500 拷贝）vs PM（150）正确分开且按比例分配读；
- 实测 WGS（HG00733，vs 二倍体组装）：418 个杂合位点灵敏 87%、特异 93%（RepeatHMM 40%/84%）；扫描模式找到 146 个未注释扩增位点；
- 速度：minimap2 ~3h + Straglr ~1.5h（32 线程，~15 GB）@30×。

**局限**：
- 在 TRGT 论文基准上表现垫底（Mendelian 精确 70.05%/±1bp 90.34%）——GMM 对 HiFi 读长分布的拟合不如编辑距离聚类稳；
- 依赖外部 TRF/blastn；Python 实现偏慢；长读噪声下 motif 匹配失败即无分型。

## 对 CCS-STR 优化的启示

- **AIC 选簇数**是把「等位基因数」从硬编码 2 变成数据驱动的方法——处理三等位基因/嵌合位点时可借用；
- 「clip 读 + 参考侧翼 BLASTN 补边界」是廉价可靠的跨读补救，比完整重比对轻量；
- HiFi 上它不是首选分型器（精度落后），但**扫描模式仍是"目录之外找扩增"的第一梯队**——与 TRGT 的靶向模式互补。

## 参考

- 论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC8361843
- 仓库：https://github.com/bcgsc/straglr
