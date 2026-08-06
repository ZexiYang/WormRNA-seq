# Smart-seq3 单线虫 RNA-seq 四种定量流程（Pipeline）整理

> 整理日期：2026-08-04（2026-08-06 增加流程四 umite）
> 数据：N2_Adult-1 / N2_Adult-2（150bp PE，Smart-seq3，R1 5' 端 pattern `ATTGCGCAATG` + UMI(12-19)）
> 集群：zhllab（ircbc），所有路径均为集群路径
> 四流程独立环境与通用脚本：`/share/home/zhllab_student/singleworm_smartseq3/`

---

## 流程一：公司标准 bulk 流程（不去重、不去污染）

```
raw fq → fastp（接头/质量修剪）
       → STAR（NCBI GCF_000002985.6 / WBcel235，sorted BAM）
       → featureCounts（read 水平计数，不做 UMI 去重）
```

- **位置**：`/share/home/zhllab_student/yzx/swRNAseq_test_202606/F26A040005318_CAEofpsyT_0615/`
  （`01raw_data/` → `02clean_data/` → `03align_data/` → `04counts/`）
- **参考基因组**：NCBI RefSeq `GCF_000002985.6`（`Reference_genome/GCF_000002985.6/genomic.gtf`），
  **Geneid 为 CELE_ 序列名**（如 `CELE_Y74C9A.6`），与 WormBase WBGene ID 不同名空间
- **特点**：
  - 流程最简单，无 UMI 概念、无 rRNA / 细菌污染去除
  - 多重比对 reads 被 featureCounts 丢弃（约占全部比对的 60-70%）
  - 产出为 read 水平计数，PCR 重复未去除，绝对定量虚高
- **产出**：`04counts/<sample>_Aligned.sortedByCoord.out_read.count`（+ `.summary`）

## 流程二：自建 Smart-seq3 流程（umi_tools + HTSeq）

```
raw fq → FastQC
       → umi_tools extract（R1 5' 端 8nt UMI 移入 read 头）
       → cutadapt（接头 + 低质量修剪）
       → bowtie2 vs E. coli HT115(DE3)（去细菌，保留 unmapped）
       → bowtie2 vs 线虫 rRNA（去 rRNA，保留 unmapped）
       → STAR 2.7.11b（WormBase WS298）
       → HTSeq-count（raw read 计数）
       → umi_tools dedup → HTSeq-count（UMI 分子计数）
```

- **脚本**：`/share/home/zhllab_student/rnaseq_project/scripts/01_run_sample.sh`
  （批量提交 `02_submit_array.sbatch`，合并 `03_merge_counts.py`）
- **产出**：`/share/home/zhllab_student/rnaseq_project/results/<sample>/`
  + 合并矩阵 `/share/home/zhllab_student/rnaseq_project/counts/{raw,umi}_count_matrix.csv`
- **参考基因组**：WormBase WS298 `canonical_geneset.gtf`（Geneid 为 WBGene ID）
- **特点**：
  - 唯一做污染去除的流程（E. coli + rRNA 两步 bowtie2）
  - UMI 去重按 "UMI + 比对位置"，对 8bp UMI 的碰撞不做序列错误校正
  - cutadapt 后大量短读长被 STAR 丢弃（"too short" 约 40%）
- **工具版本**：umi_tools 1.1.6 / HTSeq 2.0.9（smartseq3 env）、cutadapt 5.2（rnaseq env）、
  bowtie2 2.5.4 / STAR 2.7.11b（aligning env）

## 流程三：zUMIs（本项目的 Smart-seq3 原生流程）

```
raw fq → fqfilter_v2.pl（识别 5' pattern ATTGCGCAATG：
           pattern reads 取 UMI(12-19)，internal reads 保留、无 UMI）
       → STAR 2.7.3a（WS298 专用索引 star_index_ws298_zumis）
       → Rsubread featureCounts 打基因标签
       → UMI hamming-1 错误校正去重（internal reads 单独计数）
       → dgecounts.rds（umicount / readcount / readcount_internal
                        × exon / inex / intron 全套矩阵）
```

- **zUMIs 版本**：v2.9.7e，`~/src/zUMIs`，自带独立环境 `~/src/zUMIs/zUMIs-env`
- **提交脚本**：
  - 项目批量：`/share/home/zhllab_student/yzx/260727_singleworm/scripts/zumis-test-cpu02.sbatch`（单样本测试）/ `zumis-batch-fat.sbatch`（785 全量）
  - 通用版：`260727_singleworm/scripts/zumis-generic.sbatch <sample> <R1> <R2> <outdir>`
- **yaml 验证配置**（2026-08-03/04 在 day11_CF_4 与 N2_Adult-1/2 上验证通过）：
  - `cDNA(23-150) + UMI(12-19)`，`find_pattern: ATTGCGCAATG`
  - `barcode_num: 1 / automatic: no / BarcodeBinning: 0 / strand: 0 / Ham_Dist: 1`
  - `additional_STAR_params: "--clip3pAdapterSeq CTGTCTCTTATACACATCT"`
- **关键补丁（勿丢）**：
  1. `~/src/zUMIs/fqfilter_v2.pl:268`：无 `BC()` 配置时 BC tag 为空 → R fread 成 NA
     → `ScanBamParam tagFilter` 拒收 NA → Counting 崩溃。补丁将空 BC 赋常量 `CELL1`。
     原文件备份 `fqfilter_v2.pl.bak.20260803`。
  2. `zUMIs-BCdetection.R`：`BarcodeBinning: 0` 时也写 `kept_barcodes_binned.txt`
     （Counting 无条件读该文件）。
  3. yaml 必须 `BarcodeBinning: 0`：单 barcode + binning=1 会让 `BCbin` 的 setnames 崩溃。
- **资源注意**：单样本 16c/90G；两个作业不要挤同一节点
  （samtools sort 内存叠加会 OOM：`couldn't allocate memory for bam_mem`）。
- **产出**：`<outdir>/zUMIs_output/expression/<sample>.dgecounts.rds` + loom + `stats/`

## 流程四：umite（Smart-seq3 UMI 提取 + 校正去重，2026-08-05 验证）

```
raw fq → umiextract（识别 5' pattern ATTGCGCAATG + 8bp UMI，--fuzzy_umi 容错；
           UMI 移入 read 头、剪掉 pattern/UMI 序列；internal reads 保留）
       → STAR 2.7.11b（未排序 BAM）
       → samtools sort -n（按 read name 排序）
       → umicount（--mm_count_primary 多重比对计主比对；
                    --UMI_correct directional hamming 校正去重；
                    exon/intron × UMI/internal 分段计数）
```

- **umite 版本**：0.1.1（pip），文档 https://github.com/leoforster/umite
- **通用脚本**：`/share/home/zhllab_student/singleworm_smartseq3/umite/scripts/run_umite.sh`
  `bash run_umite.sh <sample> <R1> <R2> <star_index> <gtf> <outdir> [threads]`
- **验证结果**（N2_Adult-1/2，2026-08-05，输出 `.../F26A040005318_CAEofpsyT_0615/05umite/`）：
  - UMI 检出率 ~48%；UE(UMI exon) ~11.5%、RE(internal exon) ~50%、D(UMI 重复) ~30%
  - umite_UMI 总量 1.15M/1.26M（与 HTSeq+umi_tools 的 1.40M/1.31M 同量级）
  - 与另三方法相关性全部 ρ≥0.95（vs htseq_UMI Spearman 0.964，vs zUMIs_reads 0.974）
  - 对比文件：`05umite/comparison/`（四方法合并 csv + 散点图 + summary）
- **关键补丁（勿丢）**：
  1. `umiextract.py` 剥 DNBSEQ 读名 `/1` `/2` 后缀——原版对 R1/R2 读名做全等校验，
     华大读名（`.../1`、`.../2`）直接报 `readname mismatch`；且 umicount 靠 read name
     配对，必须同名。备份 `umiextract.py.bak.*`。
     （`smartseq3` conda env 和 `singleworm_smartseq3/umite/env` 两处均已打）
  2. umicount 硬性要求 GTF exon 行有 `exon_id` 属性（Ensembl 有、NCBI 没有）：
     NCBI GTF 用 `singleworm_smartseq3/umite/scripts/gtf_add_exon_id.py` 转换
     （已生成 `Reference_genome/GCF_000002985.6/genomic.umite.gtf`）；
     run_umite.sh 检测到缺 exon_id 会自动转换。
- **产出**：`<outdir>/counts/umite.{UE,UI,RE,RI,D}.tsv`（样本 × 基因；
  UE=UMI exon 去重计数为主结果；D=UMI 重复 reads 用于 QC）
- **注意**：umite UMI 计数经错误校正坍缩，量级 ~1/8 于 zUMIs 未校正 distinct UMI 计数，
  属方法学差异；计数文件名固定，不同样本用不同 outdir。

## 四种流程一句话对比

| | 公司 bulk | 自建 umi_tools | zUMIs | umite |
|---|---|---|---|---|
| UMI 去重 | 无 | UMI + 比对位置 | UMI 序列 hamming-1 校正 | UMI directional hamming 校正 |
| 污染去除 | 无 | E. coli + rRNA（bowtie2） | 无（定量后按基因集过滤） | 无 |
| Smart-seq3 reads 分段 | 无 | 无（全部当 UMI reads） | 5' pattern / internal 分开定量 | UMI / internal 分开定量（exon/intron 再细分） |
| 多重比对 | 丢弃 | 唯一比对后计数 | 主比对保留计数 | 主比对保留计数（--mm_count_primary） |
| 注释 / Geneid | NCBI GCF（CELE_ 名） | WormBase WS298（WBGene） | WormBase WS298（WBGene） | 随 GTF（NCBI 则 CELE_ 名，需 exon_id） |
| 产出 | read counts | raw + UMI counts | UMI / reads / internal 全套 | UE/UI/RE/RI/D 五套 |

## ID 桥接（跨流程对比必备）

公司结果的 `CELE_` ID → WBGene 桥接文件：
`/share/home/zhllab_student/yzx/qc_check/cele2wb.tsv`
（来源：NCBI gene_info 的 LocusTag ↔ dbXrefs WormBase，覆盖 99.9%，46,890 条；
生成脚本基于本地 `~/rnaseq_project/singleworm_260727/ref/cel_gene_info.gz`）

## 对比分析脚本

- `/share/home/zhllab_student/yzx/260727_singleworm/scripts/compare_generic.R`
  （featureCounts + zUMIs + 可选 HTSeq 三方对比，输出 summary.txt / 逐基因 csv / 散点图）
- `.../F26A040005318_CAEofpsyT_0615/05umite/comparison/compare_4methods.py`
  （featureCounts + zUMIs + HTSeq + umite 四方对比；本地副本
  `~/rnaseq_project/smartseq3_umite_N2/`）

---

## 方法选择推荐（2026-08-06，基于 N2_Adult-1/2 四方法实测对比）

**结论：zUMIs 为主流程，umite 为正交验证。**

单线虫 Smart-seq3 的特点决定了这个选择：

- 单虫起始量极低 → PCR 重复严重（本批数据 UMI 重复率 ~30%），UMI 去重有价值；
- 但 UMI reads 只占 ~11.5%、internal reads 占 ~50%，只用 UMI 计数会丢掉 4/5 信号；
- 四方法定量相关性 ρ≥0.95，方法间分歧几乎都在 1-3 counts 的边缘基因。

| 方法 | 定位 | 理由 |
|---|---|---|
| **zUMIs** | **主流程**（DE、时序分析） | Smart-seq3 原生；pattern/internal 分开再合并，信号利用最全；检出基因最多（20.9k）；785 全量批量已验证。用 inex（exon+intron）矩阵最合理 |
| **umite** | **验证流程** | UMI directional 校正最严格，独有假阳性基因最少（24 个）；与 HTSeq+umi_tools 同量级互证。用 UE+RE 对关键结论复现一次即可 |
| HTSeq+umi_tools | 特定样本备用 | 唯一去 E. coli/rRNA 污染的流程，但 cutadapt 丢 ~40% 短读长、流程步骤最多最易碎 |
| featureCounts | 仅快速粗看 | 无 UMI 概念，PCR 偏好直接进计数，绝对定量虚高；不建议单独下结论 |

检出对比详情：`yzx/smartseq3_umite_N2/detection_summary.txt`
（四方法共同检出 17,509/21,536 = 81.3%；分歧集中在 1-3 counts 低表达基因；
阈值 ≥10 时四方法检出收敛到 10.6k-13.8k）。
