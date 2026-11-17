# tdb（原 TRGTdb）

> **优先级**：P1（重要）｜ **类别**：队列数据层（VCF → 可查询数据库）｜ **维护方**：Adam English（ACEnglish，PacBio；truvari 作者）｜ **状态**：活跃（pip 安装，BSD-3）
> **一句话**：把多样本 TRGT VCF 归一化成 Parquet 三表数据库（locus/allele/sample），让万级样本的 TR 队列可以像查表一样做分析——TRGT 官方文档推荐的队列方案。

## 是什么

tdb 是一个 Python 工具 + 文件格式：把 tdb 兼容 VCF（TRGT 输出）转换为紧凑数据库（`project.tdb/` 目录，内含 parquet 文件），提供 `create/append/merge/bigmerge/query` 子命令与 Python API（pandas DataFrame）。

- 仓库：https://github.com/ACEnglish/tdb （~9★，BSD-3，notebooks 有入门教程）
- 定位：TRGT 生态的"下游 cohort 层"，对标 gnomAD 式群体数据库的本地版。

## 解决了什么问题

1. **VCF 不是规范化的表**（wiki 原话）：INFO 字段违反第一范式的原子性——一列里塞了多值；单条 VCF 行混合了三种信息（Locus / Allele / Sample）。
2. **TR 的 VCF 用"替换式"表示**：CHROM/POS/REF 是整个目标区域、ALT 是完整等位基因序列——等位基因长度虽然一眼可得（ALT 长度），但跨样本比较/汇总要在非规范结构上反复解析。
3. **万级样本的体量**：TRGT VCF 单样本 92.3 MB（bgzip），万样本 = ~1 TB 未压缩 VCF；tdb 归一化后 **~11.4 MB/样本（小 87.6%）**，且各样本独立文件（可去标识化、可增量合并）。

## 如何解决的

### 1. 三表归一化（数据模型）
| 表 | 文件 | 字段（来自 TRGT FORMAT） |
|---|---|---|
| Locus | `locus.pq` | LocusID, chrom, start, end（0-based） |
| Allele | `allele.pq` | LocusID, allele_number（0=参考）, allele_length, sequence |
| Sample | `sample.<id>.pq` | LocusID, allele_number, **spanning_reads**（SD）, length_range_lower/upper（ALLR）, **average_methylation**（AM） |

### 2. 命令体系
- `tdb create -o out.tdb in.vcf.gz`（可多 VCF 一并建库）；`tdb append` 增量加样本；`tdb merge` / `bigmerge`（≥10 个 tdb，建议**每线程 5–10 GB 内存**，如 `--mem 20 --threads 4`）；
- `tdb query allele_cnts out.tdb` 等命名查询；`gtmerge` 导出基因型矩阵（TSV/CSV）。

### 3. Python API + 谓词下推
```python
import tdb
data = tdb.load_tdb("family.tdb",
    samples=["proband","mother"],
    lfilters=[("chrom","=","chr1")],
    afilters=[("allele_length",">=",1000)])
ac = tdb.allele_count(data)
```
pyarrow 过滤器下推——只加载需要的样本/染色体/等位基因长度，万样本库也能交互式查询（Jupyter notebook 演示了「每个位点 ≥3 等位基因」「每位点 ≥20× 跨读覆盖率」等自定义查询）。

## 输入输出与使用要点

| 项 | 内容 |
|---|---|
| 输入 | tdb 兼容 VCF（TRGT 原生输出；`SD/ALLR/AM` 缺失也能转换） |
| 输出 | `.tdb` 目录（parquet）；查询导出 TSV/CSV |
| 运维 | 每线程 5–10 GB 合并内存；样本文件独立 → 可并行构建、可去标识化 |

## 性能与局限

- 单样本体积 −87.6%；交互式查询靠列式存储 + 谓词下推；
- 局限：绑定 TRGT 的 VCF 字段语义（LongTR/Straglr VCF 不直接兼容）；不是分布式系统——单机内存为界的合并规模；9★ 的小社区意味着踩坑要自己啃源码。

## 对 CCS-STR 优化的启示

- **队列层要提前设计**：TR 分析的"数据爆炸"不在分型而在汇总——如果你要跑 cohort（尤其加上 AM/ALLR 这类读级衍生统计），tdb 是现成答案，自研格式前先评估；
- 三表模型（locus/allele/sample 分离）本身就是对 TR 数据本质的抽象，你的内部 schema 可以直接借用；
- `spanning_reads/length_range/average_methylation` 这三个字段是 TRGT 官方认可的"最小充分统计量"——回归测试只需要比这三个字段。

## 参考

- 仓库 + wiki（数据模型、查询教程）：https://github.com/ACEnglish/tdb
- TRGT 文档中的推荐：https://github.com/pacificbiosciences/trgt
