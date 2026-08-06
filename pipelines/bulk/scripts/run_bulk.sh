#!/bin/bash
# run_bulk.sh — 公司标准 bulk 流程（fastp → STAR → featureCounts），独立环境版
# 用法: bash run_bulk.sh <SAMPLE> <R1.fq.gz> <R2.fq.gz> <STAR_INDEX> <GTF> <OUTDIR> [THREADS]
# 例:   bash run_bulk.sh N2_Adult-1 R1.fq.gz R2.fq.gz \
#         /path/star_index /path/annotation.gtf ./out 8
set -euo pipefail

SAMPLE=$1; R1=$2; R2=$3; STAR_IDX=$4; GTF=$5; OUT=$6; THREADS=${7:-8}
ENV="$(cd "$(dirname "${BASH_SOURCE[0]}")/../env" && pwd)"   # 本流程独立环境

FASTP="$ENV/bin/fastp"
STAR="$ENV/bin/STAR"
SAMTOOLS="$ENV/bin/samtools"
FEATURECOUNTS="$ENV/bin/featureCounts"
for t in "$FASTP" "$STAR" "$SAMTOOLS" "$FEATURECOUNTS"; do
  [ -x "$t" ] || { echo "[FAIL] 缺少工具: $t"; exit 1; }
done

mkdir -p "$OUT"/{clean,align,counts,logs}
echo "[$(date '+%F %T')] START $SAMPLE (threads=$THREADS)"

# 1. fastp
"$FASTP" -i "$R1" -I "$R2" \
  -o "$OUT/clean/${SAMPLE}_1.clean.fq.gz" -O "$OUT/clean/${SAMPLE}_2.clean.fq.gz" \
  -h "$OUT/logs/${SAMPLE}_fastp.html" -j "$OUT/logs/${SAMPLE}_fastp.json" \
  --detect_adapter_for_pe -q 15 -u 40 --cut_right \
  --cut_right_window_size 4 --cut_right_mean_quality 20 \
  --length_required 25 --poly_g_min_len 10 \
  -w "$THREADS" > "$OUT/logs/fastp.log" 2>&1

# 2. STAR
"$STAR" --runThreadN "$THREADS" --genomeDir "$STAR_IDX" \
  --readFilesIn "$OUT/clean/${SAMPLE}_1.clean.fq.gz" "$OUT/clean/${SAMPLE}_2.clean.fq.gz" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT/align/${SAMPLE}_" \
  --outSAMtype BAM SortedByCoordinate > "$OUT/logs/star.log" 2>&1
"$SAMTOOLS" index -@ "$THREADS" "$OUT/align/${SAMPLE}_Aligned.sortedByCoord.out.bam"

# 3. featureCounts
"$FEATURECOUNTS" -T "$THREADS" -a "$GTF" \
  -o "$OUT/counts/${SAMPLE}_read.count" \
  "$OUT/align/${SAMPLE}_Aligned.sortedByCoord.out.bam" > "$OUT/logs/featurecounts.log" 2>&1

echo "[$(date '+%F %T')] DONE $SAMPLE -> $OUT"
