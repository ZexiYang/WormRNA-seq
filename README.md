# WormRNA-seq

This is a playground for practising RNAseq data analysis for C. elegans, starting on 20260801.

单线虫（C. elegans）Smart-seq3 RNA-seq 定量流程与四方法对比。
本仓库收录四种定量流程的**运行脚本、环境说明与关键补丁**，以及四方法对比分析脚本。
大型 conda 环境与参考基因组索引不包含在仓库内（见各流程 README 的版本清单重建）。

## 目录结构

```
pipelines/
├── bulk/        公司标准 bulk 流程: fastp → STAR → featureCounts（read 计数, 不去重）
├── umi_tools/   自建 Smart-seq3 流程: umi_tools + HTSeq（含 E. coli / rRNA 去除）
├── zumis/       zUMIs 2.9.7e（Smart-seq3 原生, UMI hamming 校正去重）
│                scripts/ 下含通用版 run_zumis.sbatch 与 785 样本批量脚本
└── umite/       umite 0.1.1（umiextract + umicount, UMI directional 校正去重）
docs/            四种定量流程整理 + zUMIs 作业提交指南
comparison/      四方法对比分析脚本（计数合并、检出分析、详细报告生成）
patches/         上游源码的关键补丁（zUMIs / umite，含应用说明）
```

## 快速开始

每个流程一个通用脚本，参数风格一致（样本名、R1/R2、STAR 索引、GTF、输出目录、线程数）：

```bash
# bulk（fastp → STAR → featureCounts）
bash pipelines/bulk/scripts/run_bulk.sh <SAMPLE> <R1.fq.gz> <R2.fq.gz> <STAR_INDEX> <GTF> <OUTDIR> [THREADS]

# umi_tools + HTSeq（含污染去除，需额外两个 bowtie2 索引）
bash pipelines/umi_tools/scripts/run_umi_tools.sh <SAMPLE> <R1> <R2> <STAR_IDX> <GTF> \
  <BT2_ECOLI_PREFIX> <BT2_RRNA_PREFIX> <OUTDIR> [THREADS]

# zUMIs（Slurm 提交；索引须用 zUMIs 自带 STAR 2.7.3a 构建）
sbatch pipelines/zumis/scripts/run_zumis.sbatch <SAMPLE> <R1> <R2> <STAR_INDEX> <GTF> <OUTDIR>

# umite（GTF 缺 exon_id 时自动转换）
bash pipelines/umite/scripts/run_umite.sh <SAMPLE> <R1> <R2> <STAR_INDEX> <GTF> <OUTDIR> [THREADS]
```

## 关键补丁（重上环境时必打，见 `patches/`）

- **zUMIs** `fqfilter_v2.pl`：bulk 配置（无 `BC()`）时 BC tag 为空导致 Counting 崩溃 →
  空 BC 赋常量 `CELL1`；`zUMIs-BCdetection.R`：`BarcodeBinning: 0` 时也写 binned 文件。
- **umite** `umiextract.py`：剥除 DNBSEQ 读名 `/1` `/2` 后缀，否则 R1/R2 读名校验报错。

## 文档

- [四种定量流程 pipeline 整理](docs/四种定量流程_pipeline整理.md) — 四流程详细说明、
  特性对比表、ID 桥接（CELE_ ↔ WBGene）、方法选择推荐（zUMIs 主流程 + umite 正交验证）
- [zUMIs 作业提交指南](docs/README_zumis_submit.md) — Slurm 提交新手向教程（785 样本批量）

## 环境版本速查

| 流程 | 核心工具版本 |
|---|---|
| bulk | fastp 0.24.0 / STAR 2.7.11b / featureCounts 2.0.8 / samtools 1.21 |
| umi_tools | umi_tools 1.1.6 / HTSeq 2.0.9 / cutadapt 5.2 / bowtie2 2.5.4 / STAR 2.7.11b |
| zUMIs | zUMIs 2.9.7e 自带环境（STAR 2.7.3a / samtools 1.9 / R 3.6.3） |
| umite | umite 0.1.1 / HTSeq 2.0.9 / RapidFuzz 3.13.0 / STAR 2.7.11b / samtools 1.21 |

参考基因组：WormBase WS298（推荐）或 NCBI GCF_000002985.6（公司流程原版，
Geneid 为 CELE_ 序列名，跨流程对比需 ID 桥接，见 docs）。
