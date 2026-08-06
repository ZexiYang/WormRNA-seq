# umite — Smart-seq3 UMI 提取 + UMI 校正去重定量（独立环境）

```
raw fq → umiextract（识别 5' pattern ATTGCGCAATG + 8bp UMI，fuzzy 容错；
           UMI 移入 read 头、剪掉 pattern/UMI 序列）
       → STAR（未排序 BAM）
       → samtools sort -n（按 read name 排序，配对相邻）
       → umicount（exon/intron 分段的 UMI 与 internal read 计数；
                    --UMI_correct 做 directional hamming 校正去重）
```

## 环境

`env/`：python 3.9、umite 0.1.1（pip wheel 安装）、HTSeq 2.0.9、regex、
RapidFuzz 3.13.0、STAR 2.7.11b、samtools 1.21
（conda-forge/bioconda 独立安装，不共享任何文件）

**环境内补丁（勿丢）**：`env/lib/python3.9/site-packages/umite/umiextract.py`
剥除 DNBSEQ 读名 `/1` `/2` 后缀（原版对 R1/R2 读名全等校验，华大读名会报
`readname mismatch`；且 umicount 靠 read name 配对，必须同名）。
原文件备份同目录 `umiextract.py.bak.*`。

## 运行

```bash
bash scripts/run_umite.sh <SAMPLE> <R1.fq.gz> <R2.fq.gz> <STAR_INDEX> <GTF> <OUTDIR> [THREADS]
```

例（C. elegans WS298）：

```bash
bash scripts/run_umite.sh N2_Adult-1 N2_Adult-1_1.fq.gz N2_Adult-1_2.fq.gz \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/star_index_ws298 \
  /share/home/zhllab_student/rnaseq_project/ref/WS298/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf \
  ./out/N2_Adult-1 8
```

## 产出

`OUTDIR/{umiextract,align,counts,logs}`。计数在 `counts/`：

- `umite.UE.tsv` / `umite.UI.tsv` — UMI 去重后 exon / intron 计数（主结果用 UE）
- `umite.RE.tsv` / `umite.RI.tsv` — internal reads（无 UMI）exon / intron 计数
- `umite.D.tsv` — UMI 重复 reads 计数（QC 用）
- 矩阵为样本 × 基因；前几列是 `_unmapped` / `_multimapping` / `_no_feature` / `_ambiguous` 类别统计

## 注意

- umicount 硬性要求 GTF 的 exon 行有 `exon_id` 属性（Ensembl 有，NCBI 没有）；
  脚本检测到缺失会自动用 `scripts/gtf_add_exon_id.py` 转换（同时补 `gene_name`）。
- GTF 染色体名必须与 STAR 索引/BAM 一致（NCBI `NC_*` 与 WormBase `I/II/...` 不通用）。
- umite 的 UMI 计数经过错误校正坍缩，量级与 HTSeq+umi_tools 相当（~1/8 于 zUMIs
  未校正的 distinct UMI 计数），属方法学差异而非数据问题。
- umicount 计数文件名固定（`umite.*.tsv`），不同样本请用不同 OUTDIR。
- umiextract 的 `-d/--output_dir` 只认**已存在的绝对路径**（相对路径会误报
  "Folder does not exist"）；run_umite.sh 内部已自动转绝对路径。
- 计算请提交 Slurm 作业跑，不要在 login node 直接运行（login node 只装环境）。
