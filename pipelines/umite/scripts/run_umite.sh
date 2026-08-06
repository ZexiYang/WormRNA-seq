#!/bin/bash
# run_umite.sh — umite Smart-seq3 流程（umiextract → STAR → samtools sort -n → umicount），独立环境版
# 用法: bash run_umite.sh <SAMPLE> <R1.fq.gz> <R2.fq.gz> <STAR_INDEX> <GTF> <OUTDIR> [THREADS]
# 例:   bash run_umite.sh N2_Adult-1 R1.fq.gz R2.fq.gz \
#         /path/star_index /path/annotation.gtf ./out/N2_Adult-1 8
# 注意: GTF 若无 exon_id 属性（如 NCBI GTF），脚本会自动调用 gtf_add_exon_id.py 转换后用转换版
set -euo pipefail

SAMPLE=$1; R1=$2; R2=$3; STAR_IDX=$4; GTF=$5; OUT=$6; THREADS=${7:-8}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV="$(cd "$HERE/../env" && pwd)"   # 本流程独立环境
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"          # umiextract 的 -d 只认绝对路径，统一转绝对

UMIEXTRACT="$ENV/bin/umiextract"
UMICOUNT="$ENV/bin/umicount"
STAR="$ENV/bin/STAR"
SAMTOOLS="$ENV/bin/samtools"
PYTHON="$ENV/bin/python"
for t in "$UMIEXTRACT" "$UMICOUNT" "$STAR" "$SAMTOOLS" "$PYTHON"; do
  [ -x "$t" ] || { echo "[FAIL] 缺少工具: $t"; exit 1; }
done

mkdir -p "$OUT"/{umiextract,align,counts,logs}
echo "[$(date '+%F %T')] START $SAMPLE (threads=$THREADS)"

# 0. GTF 检查：umicount 要求 exon 行有 exon_id 属性，缺则自动转换
if zcat -f "$GTF" | grep -m1 $'\texon\t' | grep -q 'exon_id'; then
  GTF_U="$GTF"
else
  GTF_U="$OUT/logs/$(basename "${GTF%.gtf}").umite.gtf"
  echo "[$(date '+%F %T')] GTF 缺 exon_id，转换: $GTF -> $GTF_U"
  "$PYTHON" "$HERE/gtf_add_exon_id.py" "$GTF" "$GTF_U" > "$OUT/logs/gtf_convert.log" 2>&1
fi

# 1. umiextract（Smart-seq3 8bp UMI，fuzzy 容错识别；UMI 移入 read 头并剪掉 5' pattern）
"$UMIEXTRACT" -1 "$R1" -2 "$R2" \
  -d "$OUT/umiextract" -c 1 --umilen 8 --fuzzy_umi \
  -l "$OUT/logs/umiextract.log"

B1=$(basename "$R1"); B1=${B1%.fq.gz}; B1=${B1%.fastq.gz}
B2=$(basename "$R2"); B2=${B2%.fq.gz}; B2=${B2%.fastq.gz}
R1U="$OUT/umiextract/${B1}_umiextract.fastq.gz"
R2U="$OUT/umiextract/${B2}_umiextract.fastq.gz"
for f in "$R1U" "$R2U"; do
  [ -s "$f" ] || { echo "[FAIL] umiextract 未产出: $f"; exit 1; }
done

# 2. STAR（输出未排序 BAM，UMI 在 read 头中）→ samtools 按 read name 排序
"$STAR" --runThreadN "$THREADS" --genomeDir "$STAR_IDX" \
  --readFilesIn "$R1U" "$R2U" \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$OUT/align/${SAMPLE}_" > "$OUT/logs/star.log" 2>&1
"$SAMTOOLS" sort -n -@ "$THREADS" \
  -o "$OUT/align/${SAMPLE}.namesorted.bam" \
  "$OUT/align/${SAMPLE}_Aligned.out.bam" 2> "$OUT/logs/samtools_sort.log"
rm -f "$OUT/align/${SAMPLE}_Aligned.out.bam"

# 3. umicount（多重比对计主比对；UMI directional hamming 校正去重）
"$UMICOUNT" \
  --bams "$OUT/align/${SAMPLE}.namesorted.bam" \
  --gtf "$GTF_U" \
  --mm_count_primary --UMI_correct \
  -d "$OUT/counts" -c 1 \
  -l "$OUT/logs/umicount.log"

echo "[$(date '+%F %T')] DONE $SAMPLE -> $OUT/counts (umite.UE/UI/RE/RI/D.tsv)"
