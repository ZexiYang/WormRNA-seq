#!/bin/bash
# run_umi_tools.sh — 自建 Smart-seq3 流程（umi_tools + HTSeq），独立环境版
# 步骤与 ~/rnaseq_project/scripts/01_run_sample.sh 完全一致,
# 区别仅在于所有工具（含 fastqc 依赖的 java）都来自本流程自己的环境,
# 不再跨环境借用（原版 fastqc 借用了 aligning 环境的 java）。
#
# 用法: bash run_umi_tools.sh <SAMPLE> <R1.fq.gz> <R2.fq.gz> <STAR_IDX> <GTF> \
#         <BT2_ECOLI_PREFIX> <BT2_RRNA_PREFIX> <OUTDIR> [THREADS]
set -euo pipefail

SAMPLE=$1; R1=$2; R2=$3; STAR_IDX=$4; GTF=$5; BT2_ECOLI=$6; BT2_RRNA=$7; OUT=$8
THREADS=${9:-8}
ENV="$(cd "$(dirname "${BASH_SOURCE[0]}")/../env" && pwd)"   # 本流程独立环境

CUTADAPT="$ENV/bin/cutadapt"
BOWTIE2="$ENV/bin/bowtie2"
STAR="$ENV/bin/STAR"
SAMTOOLS="$ENV/bin/samtools"
FASTP="$ENV/bin/fastp"
MULTIQC="$ENV/bin/multiqc"
UMI_TOOLS="$ENV/bin/umi_tools"
HTSEQ="$ENV/bin/htseq-count"
for t in "$CUTADAPT" "$BOWTIE2" "$STAR" "$SAMTOOLS" "$FASTP" "$MULTIQC" "$UMI_TOOLS" "$HTSEQ"; do
  [ -x "$t" ] || { echo "[FAIL] 缺少工具: $t"; exit 1; }
done
# fastqc 留作可选（需要时取消下方注释即可, 环境自带 openjdk）
export PATH="$ENV/bin:$PATH"

mkdir -p "$OUT"/{qc,umi,trim,ecoli,rrna,star,counts,logs}
LOG="$OUT/logs"
echo "[$(date '+%F %T')] START $SAMPLE (threads=$THREADS)"

# 1. fastp QC（输出写 /dev/null, 只要报告; 替代原 FastQC 步骤） -----------------
"$FASTP" -i "$R1" -I "$R2" -o /dev/null -O /dev/null \
    -h "$OUT/qc/${SAMPLE}_fastp.html" -j "$OUT/qc/${SAMPLE}_fastp.json" \
    -w "$THREADS" > "$LOG/fastp.log" 2>&1 || echo "WARN: fastp failed (non-fatal)"
# 可选: 如需 FastQC 报告, 取消下一行注释
# "$ENV/bin/fastqc" -t "$THREADS" -o "$OUT/qc" "$R1" "$R2" > "$LOG/fastqc.log" 2>&1

# 2. UMI extract (Smart-seq3: 8 nt UMI at 5' of R1) ---------------------------
"$UMI_TOOLS" extract --bc-pattern=NNNNNNNN --ignore-read-pair-suffixes \
    -I "$R1" -S "$OUT/umi/${SAMPLE}_umi_R1.fastq.gz" \
    --read2-in="$R2" --read2-out="$OUT/umi/${SAMPLE}_umi_R2.fastq.gz" \
    --log="$LOG/umi_extract.log"

# 3. cutadapt: adapters + quality ----------------------------------------------
# MGI/DNBSEQ standard adapters (data is DNBSEQ PE150)
"$CUTADAPT" -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
            -A AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
            -q 20 -m 25 -j "$THREADS" \
            -o "$OUT/trim/${SAMPLE}_trim_R1.fastq.gz" \
            -p "$OUT/trim/${SAMPLE}_trim_R2.fastq.gz" \
            "$OUT/umi/${SAMPLE}_umi_R1.fastq.gz" \
            "$OUT/umi/${SAMPLE}_umi_R2.fastq.gz" \
            > "$LOG/cutadapt.log" 2>&1

# 4. bowtie2 vs E. coli — keep only pairs that do NOT align -------------------
"$BOWTIE2" -x "$BT2_ECOLI" -p "$THREADS" --very-sensitive \
    -1 "$OUT/trim/${SAMPLE}_trim_R1.fastq.gz" \
    -2 "$OUT/trim/${SAMPLE}_trim_R2.fastq.gz" \
    --un-conc-gz "$OUT/ecoli/${SAMPLE}_noecoli_%.fastq.gz" \
    -S /dev/null 2> "$LOG/bowtie2_ecoli.log"

# 5. bowtie2 vs C. elegans rRNA — keep only pairs that do NOT align -----------
"$BOWTIE2" -x "$BT2_RRNA" -p "$THREADS" --very-sensitive \
    -1 "$OUT/ecoli/${SAMPLE}_noecoli_1.fastq.gz" \
    -2 "$OUT/ecoli/${SAMPLE}_noecoli_2.fastq.gz" \
    --un-conc-gz "$OUT/rrna/${SAMPLE}_clean_%.fastq.gz" \
    -S /dev/null 2> "$LOG/bowtie2_rrna.log"

# 6. STAR vs WBcel235 ----------------------------------------------------------
"$STAR" --runThreadN "$THREADS" \
    --genomeDir "$STAR_IDX" \
    --readFilesIn "$OUT/rrna/${SAMPLE}_clean_1.fastq.gz" "$OUT/rrna/${SAMPLE}_clean_2.fastq.gz" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outSAMattributes NH HI NM MD AS \
    --quantMode GeneCounts \
    --outFileNamePrefix "$OUT/star/${SAMPLE}_" \
    > "$LOG/star.log" 2>&1

BAM="$OUT/star/${SAMPLE}_Aligned.sortedByCoord.out.bam"
"$SAMTOOLS" index -@ "$THREADS" "$BAM"

# 7. HTSeq-count: raw reads (no UMI collapse) -----------------------------------
"$HTSEQ" -f bam -r pos -s no -t exon -i gene_id \
    "$BAM" "$GTF" > "$OUT/counts/${SAMPLE}.htseq_raw.tsv" \
    2> "$LOG/htseq_raw.log"

# 8. UMI dedup -> HTSeq-count: molecule counts ----------------------------------
"$UMI_TOOLS" dedup -I "$BAM" -S "$OUT/star/${SAMPLE}_dedup.bam" \
    --log="$LOG/umi_dedup.log"
"$SAMTOOLS" index -@ "$THREADS" "$OUT/star/${SAMPLE}_dedup.bam"
"$HTSEQ" -f bam -r pos -s no -t exon -i gene_id \
    "$OUT/star/${SAMPLE}_dedup.bam" "$GTF" \
    > "$OUT/counts/${SAMPLE}.htseq_umi.tsv" \
    2> "$LOG/htseq_umi.log"

# cleanup intermediate fastqs (keep bam + counts + logs) -----------------------
rm -rf "$OUT/umi" "$OUT/trim" "$OUT/ecoli" "$OUT/rrna"

# 9. MultiQC 汇总（fastp/cutadapt/bowtie2/STAR/umi_tools 报告 → 一个 html） -----
"$MULTIQC" "$OUT/qc" "$LOG" "$OUT/star" \
    -o "$OUT/qc" -n "${SAMPLE}_multiqc_report" \
    > "$LOG/multiqc.log" 2>&1 || echo "WARN: multiqc failed (non-fatal)"

echo "[$(date '+%F %T')] DONE $SAMPLE"
