# zumis — zUMIs 2.9.7e Smart-seq3 流程（独立环境）

```
raw fq → fqfilter（5' pattern ATTGCGCAATG / UMI(12-19) 识别）
       → STAR 2.7.3a → Rsubread 打基因标签 → UMI hamming-1 去重
       → dgecounts.rds（umicount / readcount / readcount_internal）
```

## 环境

`zUMIs/`：zUMIs 源码 + 自带独立环境 `zUMIs-env/`（STAR 2.7.3a、samtools 1.9、
pigz 2.3.4、R 3.6.3 及全部 R 依赖），从 `~/src/zUMIs` 整体复制后
用 `conda-unpack` 重定位，已验证（R 包加载、STAR/samtools 均正常）。

**包含关键补丁**（勿覆盖回未打补丁的版本）：

1. `fqfilter_v2.pl:268`：无 `BC()` 配置时 BC tag 为空 → 下游 R 读成 NA 导致
   Counting 崩溃；补丁把空 BC 赋常量 `CELL1`（原件 `fqfilter_v2.pl.bak.20260803`）。
2. `zUMIs-BCdetection.R`：`BarcodeBinning: 0` 时也写 `kept_barcodes_binned.txt`。
3. yaml 必须 `BarcodeBinning: 0`（单 barcode + binning=1 会让 BCbin 崩溃）。

## 运行

```bash
sbatch scripts/run_zumis.sbatch <sample> <R1.fq.gz> <R2.fq.gz> <STAR_INDEX> <GTF> <outdir>
```

例（C. elegans WS298，索引必须用 zUMIs 自带 STAR 2.7.3a 构建的 `star_index_ws298_zumis`）：

```bash
sbatch scripts/run_zumis.sbatch N2_Adult-1 \
  /share/home/zhllab_student/yzx/swRNAseq_test_202606/F26A040005318_CAEofpsyT_0615/01raw_data/N2_Adult-1_1.fq.gz \
  /share/home/zhllab_student/yzx/swRNAseq_test_202606/F26A040005318_CAEofpsyT_0615/01raw_data/N2_Adult-1_2.fq.gz \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/star_index_ws298_zumis \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf \
  ./output/N2_Adult-1
```

## 注意

- 单样本 16c/90G；**两个作业不要挤同一节点**（samtools sort 内存叠加 OOM）。
- STAR 索引必须兼容 2.7.3a（2.7.11b 建的索引不能用）。
- 产出 `<outdir>/zUMIs_output/expression/<sample>.dgecounts.rds`；
  已有该文件的样本会被自动跳过（断点续跑）。
