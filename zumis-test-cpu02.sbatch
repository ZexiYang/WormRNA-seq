#!/bin/bash
#SBATCH --job-name=zumis-test
#SBATCH --partition=compute_cpu
#SBATCH --exclude=cpu01,cpu02   # 2026-08-01: cpu01/cpu02 /share 挂载故障, 排除; 不指定节点让调度器自选
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=90G
#SBATCH --output=/share/home/zhllab_student/yzx/260727_singleworm/logs/zumis-test.%j.log
# =============================================================================
# Smart-seq3 zUMIs —— cpu02 单样本测试
# 用法: sbatch zumis-test-cpu02.sbatch <sample_name>
#   例: sbatch zumis-test-cpu02.sbatch day11_CF_4
# 输入: 01raw_data/<sample>/<sample>_1.fq.gz, _2.fq.gz (150bp PE, Smart-seq3)
# 输出: 06zumis/<sample>/  (zUMIs 全部输出, 含 UMI counts)
#       06zumis/yamls/<sample>.zUMIs.yaml
# 与 fat 批量版的唯一区别: 线程数 16 (32 会 OOM: ~30 路并行 STAR 超 90G) / mem_limit 80G
# =============================================================================
set -euo pipefail
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY all_proxy

SAMPLE=${1:?用法: sbatch zumis-test-cpu02.sbatch <sample_name>}
THREADS=16
MEMLIMIT=80   # zUMIs mem_limit (GB), 供 STAR 排序使用; 线虫基因组小, 80G 足够

BASE=/share/home/zhllab_student/yzx/260727_singleworm
REF=/share/home/zhllab_student/rnaseq_project/ref/WS298
R1=$BASE/01raw_data/$SAMPLE/${SAMPLE}_1.fq.gz
R2=$BASE/01raw_data/$SAMPLE/${SAMPLE}_2.fq.gz
OUT=$BASE/06zumis/$SAMPLE
YAML=$BASE/06zumis/yamls/$SAMPLE.zUMIs.yaml

# ---- 输入检查 ----
for f in "$R1" "$R2" "$REF/star_index_ws298_zumis/SAindex" \
         "$REF/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf"; do
  [ -s "$f" ] || { echo "[FAIL] 文件不存在: $f"; exit 1; }
done
mkdir -p "$OUT" "$BASE/06zumis/yamls"

# 完成标记: 已生成最终 UMI counts 则跳过
if [ -s "$OUT/zUMIs_output/expression/${SAMPLE}.dgecounts.rds" ]; then
  echo "[SKIP] $SAMPLE 已完成 ($OUT)"; exit 0
fi

# ---- 生成 yaml (Smart-seq3 标准配置; 必须 ASCII-only, 集群 locale 为 C) ----
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

echo "[INFO] $(date) 开始 zUMIs: $SAMPLE (threads=$THREADS, mem_limit=${MEMLIMIT}G)"
"$HOME/src/zUMIs/zUMIs.sh" -c -y "$YAML"
echo "[DONE] $(date) $SAMPLE -> $OUT"
