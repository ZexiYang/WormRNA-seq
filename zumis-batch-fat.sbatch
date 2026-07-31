#!/bin/bash
#SBATCH --job-name=zumis-batch
#SBATCH --partition=compute_fat
#SBATCH --array=1-785%2
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --output=/share/home/zhllab_student/yzx/260727_singleworm/logs/zumis-batch.%A_%a.log
# =============================================================================
# Smart-seq3 zUMIs —— fat 节点全量批量版 (测试通过后再提交!)
# 用法: sbatch zumis-batch-fat.sbatch
#   --array=1-785%2: 785 个样本每个一个 array 任务, 最多 2 个并发
#                    (fat01+fat02 各占一个 128 核任务)
# 断点续跑: 已生成 dgecounts.rds 的样本自动 [SKIP], 重复提交安全
# 样本清单: scripts/samples.txt (新数据用 gen_samples_txt.sh 重新生成,
#           并把上面 --array=1-785 改成新的行数)
# =============================================================================
set -euo pipefail
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY all_proxy

THREADS=128
MEMLIMIT=480

BASE=/share/home/zhllab_student/yzx/260727_singleworm
REF=/share/home/zhllab_student/rnaseq_project/ref/WS298

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BASE/scripts/samples.txt")
if [ -z "$SAMPLE" ]; then
  echo "[SKIP] array id $SLURM_ARRAY_TASK_ID 无对应样本"; exit 0
fi

R1=$BASE/01raw_data/$SAMPLE/${SAMPLE}_1.fq.gz
R2=$BASE/01raw_data/$SAMPLE/${SAMPLE}_2.fq.gz
OUT=$BASE/06zumis/$SAMPLE
YAML=$BASE/06zumis/yamls/$SAMPLE.zUMIs.yaml

for f in "$R1" "$R2"; do
  [ -s "$f" ] || { echo "[FAIL] 文件不存在: $f"; exit 1; }
done
mkdir -p "$OUT" "$BASE/06zumis/yamls"

if [ -s "$OUT/zUMIs_output/expression/${SAMPLE}.dgecounts.rds" ]; then
  echo "[SKIP] $SAMPLE 已完成"; exit 0
fi

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
  STAR_index: $REF/star_index_ws298_zumis
  GTF_file: $REF/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf
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
    num_bases: 1
    phred: 20
  UMI_filter:
    num_bases: 1
    phred: 20
barcodes:
  barcode_num: null
  barcode_file: null
  barcode_sharing: null
  automatic: yes
  BarcodeBinning: 1
  nReadsperCell: 100
  demultiplex: no
counting_opts:
  introns: yes
  intronProb: no
  downsampling: 0
  strand: 1
  Ham_Dist: 1
  velocyto: no
  primaryHit: yes
  multi_overlap: no
  fraction_overlap: 0
  twoPass: no
make_stats: yes
which_Stage: Filtering
EOF

echo "[INFO] $(date) 开始 zUMIs: $SAMPLE (array=$SLURM_ARRAY_TASK_ID, node=$(hostname))"
"$HOME/src/zUMIs/zUMIs.sh" -c -y "$YAML"
echo "[DONE] $(date) $SAMPLE"
