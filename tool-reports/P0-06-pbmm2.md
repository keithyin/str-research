# pbmm2

> **优先级**：P0（核心必读）｜ **类别**：比对 ｜ **维护方**：PacBio ｜ **状态**：活跃，v26.x 线（2026-11）
> **一句话**：PacBio 官方的 minimap2 前端（C API 包装而非 fork），`HIFI/CCS` preset 默认开启一组**模仿 BLASR 的强制参数**——理解这些参数是排查「STR spanning read 丢失」的前提。

## 是什么

pbmm2 是 PacBio 对 minimap2 的 C++ 包装：原生 PacBio 数据进出（ BAM 内标签兼容 SMRT 生态）、内置 on-the-fly 排序、preset 化参数。内置 minimap2 v2.26（≥1.13.0；更早版本 v2.15）。BSD-3-Clause-Clear。

- 仓库：https://github.com/PacificBiosciences/pbmm2 （~224★）
- 定位：官方 **BLASR 替代品**（序列一致性、mappable bases、速度均优于 BLASR）。

## 解决了什么问题

1. **裸 minimap2 的输出不合 PacBio 生态**：SMRT Link/下游工具期望 `=/X` CIGAR、BLASR 式行为与 PacBio 标签体系。
2. **参数选择门槛**：minimap2 的参数空间大，preset 让「HiFi 比对怎么配」有唯一官方答案（`--preset HIFI`，v1.12.0 起即默认 CCS/HIFI）。
3. **重复序列的 query 重叠问题**：比对器允许 primary 与 supplementary 共享 query 区间（最多 50%），违反「一个碱基只比对一次」——在 LINE/转座/重复插入场景污染下游分析。

## 如何解决的

### 1. Preset 体系
| Preset | 用途 | 备注 |
|---|---|---|
| `SUBREAD` | CLR/subread 比对 | 开启 homopolymer-compressed k-mer |
| `CCS` / `HIFI`（默认） | HiFi 比对 | v1.14.99 移除了 HiFi preset 的最小 DP 长度阈值；**HIFI 只是 CCS 的别名**，同一参数集 |
| `ISOSEQ` | 转录本/可变剪接 | 允许覆盖 `-G`（最大内含子）、`-C` 等剪接参数 |
| `UNROLLED` | 滚环 DNA 等 | 同开 homopolymer-compressed k-mer |

HiFi preset 可覆盖项：`-k/-w`（k-mer/窗口）、`-u`（关 homopolymer 压缩）、`-A/-B`（match/mismatch 分）、`-z/-Z`（Z-drop）、`-r`（bandwidth）、`-g`（chain 终止）、`-o/-O/-e/-E/-L`（分段线性 gap 罚分）。注意：**用索引文件（.mmi）时 `-k/-w/-u` 不可覆盖**。

### 2. 强制（隐式）参数——模仿 BLASR
- **`-Y` 软剪切**：保留 soft-clipped 序列在 CIGAR 中（对 STR 至关重要：扩增读的 clip 部分不被丢弃）；
- **`--eqx`**：CIGAR 用 `=/X` 而非 `M`（下游可精确区分 match/mismatch 列）；
- **`-L`**：长 CIGAR 支持（`CG` 标签）；
- **repeated-matches trimming**：primary 与 supplementary 的 query 重叠区按 minimap2 顺序修剪非 primary 比对，现场重写 CIGAR 与参考坐标（`rm` 标签）——提升 mappable bases 与 SV 召回；
- **默认关闭 secondary 比对**（新版可用 `--secondary --max-secondary-alns` 打开，默认上限 5）。

### 3. 一致性过滤（STR 排查的关键嫌疑）
三条过滤 AND 叠加：`--min-concordance-perc`（遗留，隐藏）、`--min-id-perc`（BLAST 一致性，隐藏）、**`--min-gap-comp-id-perc`（默认 70%，把 indel 算作单一事件的 gap-compressed 一致性）**。
- 社区报告过：重复区 spanning read 因 gap-compressed 一致性不达标而被过滤成 unmapped 的现象——**排查 TRGT「spanning read 偏少」时，先用裸 minimap2 对照**（TRGT issue 区的常见处方）；
- 相关标签：`mc`/`mg`/`mi`（一致性度量）。

### 4. 工程细节
- 内置 samtools 排序（推荐 on-the-fly，排序线程最多 ~8 有收益）；
- 默认不输出 unmapped（`--unmapped` 保留）；
- v2.26 内置 minimap2 的 2.26 引入了若干比对改进（HiFi 场景的 gap 代价模型更好）。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | subreads/CCS/HiFi BAM/FASTQ + 参考 FASTA |
| 输出 | 排序好的 PacBio BAM（兼容索引并行生成 `--parallel-bai`？——26.x 加入并行 .bai 生成） |
| STR 流水推荐 | `pbmm2 align --preset HIFI ref.fa hifi.bam out.sorted.bam`（与 TRGT 官方指南一致：**WGS 与靶向用同一 preset**） |

## 性能与局限

- 官方基准：一致性、mappable bases、runtime 均优于 BLASR；
- 局限：强制参数组合使行为与裸 minimap2 有微妙差异（尤其重复区）——**做方法对比/基准时要固定比对器与版本**；`min-gap-comp-id-perc` 是隐性过滤源；索引文件会锁定部分参数。

## 对 CCS-STR 优化的启示

- **比对是 STR 分型的隐形瓶颈**：TRGT 的 spanning read 全部来自「两侧翼有锚」的比对结果；pbmm2 的 `-Y`/`--eqx`/trimming/gap 一致性过滤共同决定了谁被算作 spanning read；
- 建议在你的回归测试里加一个「裸 minimap2 v2.26 对照」通道：两者 spanning read 数差异大的位点，就是比对层优化（或个性化参数）的目标位点；
- 版本敏感：v1.14.99 移除 HiFi 最小 DP 长度阈值这类变更会静默改变结果——升级比对器必须重跑 STR 基准；
- 若要在重复区做「多比对策略投票」（如 STRdust 的切除重比对），注意 pbmm2 默认丢 secondary，需要 `--secondary` 才能拿到备选比对。

## 参考

- 仓库/README（preset、强制参数、FAQ）：https://github.com/PacificBiosciences/pbmm2
- TRGT 官方分析指南（推荐 HIFI preset）：https://github.com/pacificbiosciences/trgt
- 社区排查案例：pbmm2/pbbioconda issues（spanning read 不比对 → 对照裸 minimap2）
