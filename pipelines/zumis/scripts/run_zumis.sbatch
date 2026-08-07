#!/bin/bash
#SBATCH --job-name=zumis-sws3
#SBATCH --partition=compute_cpu
#SBATCH --exclude=cpu01,cpu02
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=90G
#SBATCH --output=/share/home/zhllab_student/singleworm_smartseq3/zumis/logs/zumis.%j.log
# =============================================================================
# Smart-seq3 zUMIs —— 独立环境版（singleworm_smartseq3/zumis/）
# 用法: sbatch run_zumis.sbatch <sample> <R1.fq.gz> <R2.fq.gz> <STAR_INDEX> <GTF> <outdir>
#   例: sbatch run_zumis.sbatch N2_Adult-1 R1.fq.gz R2.fq.gz \
#         /share/home/zhllab_student/rnaseq_project/ref/WS298/star_index_ws298_zumis \
#         /share/home/zhllab_student/rnaseq_project/ref/WS298/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf \
#         ./output/N2_Adult-1
# 注意: STAR 索引必须用 zUMIs 自带 STAR 2.7.3a（或兼容版本）构建;
#       单样本作业不要两个挤同一节点（samtools sort 内存叠加会 OOM）。
# =============================================================================
set -euo pipefail
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY all_proxy

SAMPLE=${1:?用法: sbatch run_zumis.sbatch <sample> <R1> <R2> <STAR_INDEX> <GTF> <outdir>}
R1=${2:?}; R2=${3:?}; STAR_IDX=${4:?}; GTF=${5:?}; OUT=${6:?}
THREADS=16
MEMLIMIT=80

ZUMISDIR=/share/home/zhllab_student/singleworm_smartseq3/zumis/zUMIs
YAML=$OUT/../yamls/$SAMPLE.zUMIs.yaml

for f in "$R1" "$R2" "$STAR_IDX/SAindex" "$GTF"; do
  [ -s "$f" ] || { echo "[FAIL] 文件不存在: $f"; exit 1; }
done
mkdir -p "$OUT" "$(dirname "$YAML")"

if [ -s "$OUT/zUMIs_output/expression/${SAMPLE}.dgecounts.rds" ]; then
  echo "[SKIP] $SAMPLE 已完成 ($OUT)"; exit 0
fi

# ---- Smart-seq3 验证配置 (2026-08-03/04 验证; 必须 ASCII-only) ----
cat > "$YAML" <<EOF
project: $SAMPLE
sequence_files:
  file1:
    name: $R1
    base_definition:
      - cDNA(23-150)
      - UMI(12-19)
    find_pattern: ATTGCGCAATG
  file2:
    name: $R2
    base_definition:
      - cDNA(1-150)
reference:
  STAR_index: $STAR_IDX
  GTF_file: $GTF
  exon_extension: no
  extension_length: 0
  scaffold_length_min: 0
  additional_files:
  additional_STAR_params: "--clip3pAdapterSeq CTGTCTCTTATACACATCT"
out_dir: $OUT
num_threads: $THREADS
mem_limit: $MEMLIMIT
filter_cutoffs:
  BC_filter:
    num_bases: 3
    phred: 20
  UMI_filter:
    num_bases: 3
    phred: 20
barcodes:
  barcode_num: 1
  barcode_file: null
  barcode_sharing: null
  automatic: no
  BarcodeBinning: 0
  nReadsperCell: 100
  demultiplex: no
counting_opts:
  introns: yes
  intronProb: no
  downsampling: 0
  strand: 0
  Ham_Dist: 1
  velocyto: no
  primaryHit: yes
  multi_overlap: no
  fraction_overlap: 0
  twoPass: no
make_stats: yes
which_Stage: Filtering
EOF

echo "[INFO] $(date) 开始 zUMIs: $SAMPLE (threads=$THREADS, mem_limit=${MEMLIMIT}G)"
"$ZUMISDIR/zUMIs.sh" -c -y "$YAML"
echo "[DONE] $(date) $SAMPLE -> $OUT"
