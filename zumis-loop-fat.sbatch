#!/bin/bash
#SBATCH --job-name=zumis-loop
#SBATCH --partition=compute_fat
#SBATCH --time=14-00:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --output=/share/home/zhllab_student/yzx/260727_singleworm/logs/zumis-loop.%j.log
# =============================================================================
# Smart-seq3 zUMIs —— fat 节点循环版 (单个作业内遍历样本清单)
# 用法:
#   sbatch zumis-loop-fat.sbatch            # 跑 samples.txt 全部样本
#   sbatch zumis-loop-fat.sbatch 1 392      # 只跑第 1~392 行
#   sbatch zumis-loop-fat.sbatch 393 785    # 只跑第 393~785 行
# 推荐: 提交两次(两个范围各一个), fat01/fat02 各占一个, 时间减半:
#   sbatch zumis-loop-fat.sbatch 1 392
#   sbatch zumis-loop-fat.sbatch 393 785
# 断点续跑: 已生成 dgecounts.rds 的样本自动 [SKIP]; 作业被杀后重新提交即可
# 与 array 版的区别: array 是每个样本一个 Slurm 任务(天然并行+单任务隔离),
#   循环版是一个作业独占节点依次处理(逻辑简单, 但单个样本崩溃会中断后续,
#   已用 || 兜底: 单样本失败记日志并继续下一个)
# =============================================================================
set -uo pipefail   # 注意: 不开 -e, 单样本失败不中断整个循环
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY all_proxy

THREADS=128
MEMLIMIT=480
BASE=/share/home/zhllab_student/yzx/260727_singleworm
REF=/share/home/zhllab_student/rnaseq_project/ref/WS298
LIST=$BASE/scripts/samples.txt

TOTAL=$(wc -l < "$LIST")
START=${1:-1}
END=${2:-$TOTAL}
echo "[INFO] $(date) 节点 $(hostname), 处理 samples.txt 第 $START~$END 行 (共 $TOTAL)"

for i in $(seq "$START" "$END"); do
  SAMPLE=$(sed -n "${i}p" "$LIST")
  [ -n "$SAMPLE" ] || { echo "[WARN] 第 $i 行为空, 跳过"; continue; }

  R1=$BASE/01raw_data/$SAMPLE/${SAMPLE}_1.fq.gz
  R2=$BASE/01raw_data/$SAMPLE/${SAMPLE}_2.fq.gz
  OUT=$BASE/06zumis/$SAMPLE
  YAML=$BASE/06zumis/yamls/$SAMPLE.zUMIs.yaml

  if [ -s "$OUT/zUMIs_output/expression/${SAMPLE}.dgecounts.rds" ]; then
    echo "[SKIP] ($i/$END) $SAMPLE 已完成"; continue
  fi
  if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
    echo "[FAIL] ($i/$END) $SAMPLE 缺少原始 fq, 跳过"; continue
  fi
  mkdir -p "$OUT" "$BASE/06zumis/yamls"

  cat > "$YAML" <<EOF
project: $SAMPLE
sequence_files:
  file1:
    name: $R1
    base_definition:
      - BC(1-11)
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
  barcode_num: 1
  barcode_file: null
  barcode_sharing: null
  automatic: no
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

  echo "[RUN] ($i/$END) $(date) $SAMPLE"
  if "$HOME/src/zUMIs/zUMIs.sh" -c -y "$YAML"; then
    echo "[OK] ($i/$END) $(date) $SAMPLE"
  else
    echo "[FAIL] ($i/$END) $(date) $SAMPLE (继续下一个)"
  fi
done
echo "[ALL DONE] $(date) 范围 $START~$END 处理完毕"
