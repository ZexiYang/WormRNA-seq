# bulk — 公司标准 bulk 流程（独立环境）

```
raw fq → fastp（接头/质量修剪）→ STAR → featureCounts（read 计数, 不去重）
```

## 环境

`env/`：fastp 0.24.0、STAR 2.7.11b、subread(featureCounts) 2.0.8、samtools 1.21
（bioconda 独立安装，与原 aligning 环境同版本，不共享任何文件）

## 运行

```bash
bash scripts/run_bulk.sh <SAMPLE> <R1.fq.gz> <R2.fq.gz> <STAR_INDEX> <GTF> <OUTDIR> [THREADS]
```

例（C. elegans WS298）：

```bash
bash scripts/run_bulk.sh N2_Adult-1 N2_Adult-1_1.fq.gz N2_Adult-1_2.fq.gz \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/star_index_ws298 \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf \
  ./out/N2_Adult-1 8
```

## 产出

`OUTDIR/{clean,align,counts,logs}`，计数在 `counts/<sample>_read.count`。

## 注意

- 公司原版用的是 NCBI GCF_000002985.6 注释（Geneid 为 CELE_ 序列名）；
  与 WBGene 结果对比需用 `yzx/qc_check/cele2wb.tsv` 桥接。
- 计算请提交 Slurm 作业跑，不要在 login node 直接运行（login node 只装环境）。
