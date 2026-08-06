# umi_tools — 自建 Smart-seq3 流程（独立环境）

```
raw fq → fastp QC（fastqc 可选）→ umi_tools extract(8nt UMI) → cutadapt
       → bowtie2 vs E. coli（去细菌）→ bowtie2 vs rRNA（去 rRNA）
       → STAR → HTSeq(raw) → umi_tools dedup → HTSeq(UMI)
       → MultiQC 汇总报告
```

## 环境

`env/`：fastp 0.24.0 + multiqc（QC 环节）、fastqc 0.12.1（可选，自带 openjdk）、
umi_tools 1.1.6、htseq 2.0.9、cutadapt 5.2、bowtie2 2.5.4、STAR 2.7.11b、samtools 1.21

与原流程的差异：工具路径全部指向本环境；QC 步骤由 FastQC 改为 fastp 报告
（无 java 依赖），末尾用 MultiQC 汇总全流程日志为一个 html；
如需 FastQC 可在脚本中取消注释（环境自带 openjdk）。

## 运行

```bash
bash scripts/run_umi_tools.sh <SAMPLE> <R1.fq.gz> <R2.fq.gz> <STAR_IDX> <GTF> \
  <BT2_ECOLI_PREFIX> <BT2_RRNA_PREFIX> <OUTDIR> [THREADS]
```

例（C. elegans WS298）：

```bash
bash scripts/run_umi_tools.sh N2_Adult-1 N2_Adult-1_1.fq.gz N2_Adult-1_2.fq.gz \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/star_index_ws298 \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf \
  /share/home/zhllab_student/rnaseq_project/ref/bowtie2_ecoli/ecoli \
  /share/home/zhllab_student/rnaseq_project/ref/bowtie2_rrna/rrna \
  ./out/N2_Adult-1 8
```

## 产出

`OUTDIR/counts/<sample>.htseq_raw.tsv`（read 计数）与
`OUTDIR/counts/<sample>.htseq_umi.tsv`（UMI 分子计数），
中间 fastq 自动清理，bam 保留在 `star/`。

## 注意

计算请提交 Slurm 作业跑，不要在 login node 直接运行。
