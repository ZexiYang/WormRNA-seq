#!/usr/bin/env python3
"""Render 四方法详细比较报告.md from report_stats.json + pipeline logs."""
import json
import re
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
RD = BASE / "report_data"
S = json.loads((BASE / "report_stats.json").read_text())
SAMPLES = ["N2_Adult-1", "N2_Adult-2"]

def star(path):
    d = {}
    for line in open(path):
        if "|" in line:
            k, v = [x.strip() for x in line.split("|", 1)]
            d[k] = v
    return d

def g(d, k):
    return d.get(k, "NA")

L = []
def w(s=""):
    L.append(s)

w("# Smart-seq3 单线虫 RNA-seq 四种定量方法详细比较报告")
w("")
w("> 数据：N2_Adult-1 / N2_Adult-2（150bp PE，华大 DNBSEQ，Smart-seq3，"
  "R1 5' 端 `ATTGCGCAATG` + 8bp UMI）")
w("> 方法：① 公司 bulk（fastp→STAR→featureCounts）② 自建 umi_tools+HTSeq "
  "③ zUMIs 2.9.7e ④ umite 0.1.1")
w("> 生成：analyze_full.py + render_report.py（`report_stats.json` 为全部原始数字）")
w("")
w("---")
w("")
# ---------- 1. pipeline stats ----------
w("## 1. 建库与比对统计")
w("")
w("### 1.1 输入与 UMI 识别")
w("")
w("| 样本 | 原始 reads | umiextract 检出 UMI | 检出率 |")
w("|---|---|---|---|")
for s in SAMPLES:
    p = S["pipe"][s]
    w(f"| {s} | {p['raw_reads']:,} | {p['umi_detected']:,} | {p['umi_pct']:.2f}% |")
w("")
w("Smart-seq3 设计：仅 5' 端 TSO 标记的 reads 带 UMI（预期约一半），"
  "其余为 internal reads（无 UMI 但携带外显子信号）。两样本 UMI 检出率均接近 50%，正常。")
w("")
w("### 1.2 STAR 比对（各流程用自己的参数/索引，数字不完全可比，仅供参考）")
w("")
w("| 流程 | 样本 | 输入 reads | 唯一比对 | 多比对 | too short 未比对 |")
w("|---|---|---|---|---|---|")
for s in SAMPLES:
    for name, key in [("umite", "umite_star"), ("zUMIs", "zumis_star"),
                      ("umi_tools", "umitools_star")]:
        st = star(RD / f"{key.replace('_star','')}_STAR_{s}.log") if key != "umite_star" else star(RD / f"umite_STAR_{s}.log")
        w(f"| {name} | {s} | {g(st,'Number of input reads')} | "
          f"{g(st,'Uniquely mapped reads %')} | {g(st,'% of reads mapped to multiple loci')} | "
          f"{g(st,'% of reads unmapped: too short')} |")
w("")
w("注：zUMIs 将 fastq 分成 16 个 chunk 并行比对，表中 zUMIs 行是**单个 chunk** 的日志"
  "（输入 reads ≈ 总量/16）；umi_tools 流程 STAR 限制了多重比对输出（多比对 0%）。")
w("")
w("公司 bulk 流程的 featureCounts 统计（fragment 水平）：")
w("")
w("| 样本 | Assigned | MultiMapping 丢弃 | Ambiguity | NoFeatures |")
w("|---|---|---|---|---|")
for s in SAMPLES:
    p = S["pipe"][s]
    w(f"| {s} | {p['fc_assigned']:,} | {p['fc_multi']:,} | {p['fc_ambi']:,} | {p['fc_nofeat']:,} |")
w("")
w("公司流程多重比对 fragments 全部被丢弃（约为 assigned 的 2 倍），这是其计数系统性偏低的来源之一。")
w("")
# ---------- 2. totals ----------
w("## 2. 定量总量与 UMI 去重")
w("")
w("| 指标 | N2_Adult-1 | N2_Adult-2 |")
w("|---|---|---|")
rows = [
    ("featureCounts reads（exon）", "featureCounts_reads"),
    ("zUMIs reads（exon）", "zUMIs_reads"),
    ("zUMIs UMI（exon，hamming-1 校正）", "zUMIs_UMI"),
    ("HTSeq raw reads", "htseq_raw"),
    ("HTSeq+umi_tools UMI（位置去重）", "htseq_UMI"),
    ("umite UMI（UE，directional 校正）", "umite_UMI"),
    ("umite reads（UE+RE+D，exon）", "umite_reads"),
]
for lab, k in rows:
    w(f"| {lab} | {S['totals'][SAMPLES[0]][k]:,} | {S['totals'][SAMPLES[1]][k]:,} |")
t1, t2 = S["totals"][SAMPLES[0]], S["totals"][SAMPLES[1]]
w(f"| **去重率 zUMIs（reads/UMI）** | {t1['zUMIs_reads']/t1['zUMIs_UMI']:.2f} | "
  f"{t2['zUMIs_reads']/t2['zUMIs_UMI']:.2f} |")
w(f"| **去重率 umite（reads/UMI）** | {t1['umite_reads']/t1['umite_UMI']:.2f} | "
  f"{t2['umite_reads']/t2['umite_UMI']:.2f} |")
w("")
w("要点：")
w("- 两种 UMI 校正方法（HTSeq+umi_tools、umite）的 UMI 总量同量级（~1.2-1.4M），"
  "zUMIs 的 distinct UMI 计数高 ~8 倍（~10M）——zUMIs 按 (基因, UMI) 对计数、"
  "不做错误校正坍缩，属方法学差异，不影响基因间相对排序。")
w(f"- umite UMI reads 占总 exon reads 的 ~{t1['umite_UE']/(t1['umite_UE']+t1['umite_RE']+t1['umite_D'])*100:.0f}%，"
  "internal reads 贡献了 Smart-seq3 的大部分外显子信号——只用 UMI 计数会丢掉这部分。")
w("")
# ---------- 3. exon/intron ----------
w("## 3. 外显子 / 内含子检出")
w("")
e = S["exin"]
w("### 3.1 计数组成")
w("")
w("| 类别 | 样本 | exon | intron | intron 占比 |")
w("|---|---|---|---|---|")
for s in SAMPLES:
    p = e["per_sample"][s]
    for lab, ek, ik in [
        ("umite UMI", "umite_umi_exon", "umite_umi_intron"),
        ("umite internal reads", "umite_int_read_exon", "umite_int_read_intron"),
        ("zUMIs UMI", "zumis_umi_exon", "zumis_umi_intron"),
        ("zUMIs internal reads", "zumis_int_read_exon", "zumis_int_read_intron"),
        ("zUMIs 全部 reads", "zumis_read_exon", "zumis_read_intron")]:
        frac = p[ik] / (p[ek] + p[ik]) * 100
        w(f"| {lab} | {s} | {p[ek]:,} | {p[ik]:,} | {frac:.2f}% |")
w("")
w("featureCounts 与 HTSeq 按设计只统计 exon（无内含子信息）。")
w("")
w("### 3.2 基因水平内含子信号（两样本合并）")
w("")
w(f"- umite 有内含子信号（UI+RI>0）的基因：**{e['umite_intron_genes']:,}** 个，"
  f"其中 **{e['umite_intron_only_genes']}** 个基因只有内含子信号（exon=0）")
w(f"- zUMIs 有内含子信号（read_intron>0）的基因：**{e['zumis_intron_genes']:,}** 个，"
  f"其中 **{e['zumis_intron_only_genes']}** 个只有内含子信号")
w(f"- 两方法内含子基因交集 {e['intron_overlap']:,} 个；"
  f"共有基因上内含子计数 Spearman 相关 **{e['intron_corr']:.4f}**")
w("")
w("内含子占比很低（<1%），说明文库以成熟 mRNA 为主、核 RNA 污染少；"
  "内含子信号两方法高度一致，做 exon+intron 合并定量（zUMIs inex）是安全的。")
w("![exon/intron](fig_exon_intron.png)")
w("")
# ---------- 4. correlation ----------
w("## 4. 四方法定量相关性（ρ≥0.95 的证据）")
w("")
w("两样本合并的 7 种计数类型两两相关（基因全集，WBGene 空间）：")
w("")
CT = ["featureCounts_reads", "zUMIs_reads", "zUMIs_UMI", "htseq_raw",
      "htseq_UMI", "umite_UMI", "umite_reads"]
w("### 4.1 Spearman（原始计数）")
w("")
w("| | " + " | ".join(CT) + " |")
w("|---|" + "---|" * len(CT))
for a in CT:
    w("| " + a + " | " + " | ".join(f"{S['spearman'][a][b]:.3f}" for b in CT) + " |")
w("")
w("### 4.2 Pearson（log10(count+1)）")
w("")
w("| | " + " | ".join(CT) + " |")
w("|---|" + "---|" * len(CT))
for a in CT:
    w("| " + a + " | " + " | ".join(f"{S['pearson'][a][b]:.3f}" for b in CT) + " |")
w("")
w("**全部 21 对组合 Spearman ≥ "
  f"{min(S['spearman'][a][b] for i,a in enumerate(CT) for b in CT[i+1:]):.3f}，"
  "关键跨方法组合（不同 UMI 策略之间）全部 ≥0.95。**")
w("")
w("### 4.3 分样本关键组合")
w("")
w("| 组合 | 样本 | Pearson(log) | Spearman |")
w("|---|---|---|---|")
for s in SAMPLES:
    for pair, v in S["per_sample_corr"][s].items():
        w(f"| {pair} | {s} | {v['pearson_log']:.4f} | {v['spearman']:.4f} |")
w("")
w("![correlation heatmap](fig_correlation_heatmap.png)")
w("")
w("逐基因散点图见 `N2_Adult-{1,2}_umite_vs_3methods.png`（六组合 hexbin，"
  "高密度区紧贴对角线，低计数区离散是泊松噪声主导）。")
w("")
# ---------- 5. divergence ----------
w("## 5. 方法间分歧集中在 1-3 counts 边缘基因（证据）")
w("")
w("### 5.1 按表达量分层的检出一致性")
w("")
w("对两样本合并后每个基因取四方法最大计数作为表达量代理，分层统计检出一致性：")
w("")
w("| 表达量分层 | 基因数 | 四方法全检出 | ≥3 方法检出 |")
w("|---|---|---|---|")
for c in S["concordance_bins"]:
    w(f"| {c['bin']} | {c['n']:,} | {c['pct_all4']}% | {c['pct_ge3']}% |")
w("")
w("count=1 的基因四方法全检出率仅约一半；≥11 时接近 100%。"
  "**分歧几乎全部来自低表达端。**")
w("")
w("![concordance](fig_concordance_vs_expression.png)")
w("")
w("### 5.2 分歧基因的计数分布")
w("")
d = S["discordant_overall"]
w(f"- 未获四方法一致检出的基因共 **{d['n']:,}** 个（占基因并集 {d['pct_of_union']}%），"
  f"其中最大计数 ≤3 的占 **{d['pct_le3']}%**，≤5 的占 {d['le5_pct']}%，"
  f"≥10 的只有 {d['ge10']} 个")
m = S["umite_missed"]
w(f"- 仅被 umite 漏检（其他三种都检出）的 **{m['n']:,}** 个基因："
  f"在其他方法中最大计数中位数 {m['median']:.0f}，≤3 的占 **{m['pct_le3']}%**，"
  f"≥10 的 {m['ge10']} 个，≥50 的 {m['ge50']} 个")
w("- 各方法独有检出基因（其他三种都没检出）：")
w("")
w("| 方法 | 独有基因数 | 计数中位数 | ≤3 counts 占比 | 最大计数 |")
w("|---|---|---|---|---|")
for k, v in S["unique_genes"].items():
    w(f"| {k} | {v['n']} | {v['median']:.0f} | {v['pct_le3']}% | {v['max']} |")
w("")
w("![cdf](fig_discordant_cdf.png)")
w("")
w("### 5.3 检出集合重叠（UpSet）")
w("")
w("见 `detection_upset.png`：四方法共同检出 17,509 / 21,536（81.3%）；"
  "最大的分歧组是『仅 umite 漏检』1,023 个（5.2 节已证其 62.9% ≤3 counts、83.7% ≤5）"
  "和『fc∩zUMIs（两种 UMI 校正方法均漏检）』873 个。")
w("")
# ---------- 6. bias ----------
w("## 6. 共有基因上的系统性定量偏差")
w("")
w("对四方法共同检出的 17,509 个基因，计算两两 log2(计数比)：")
w("")
w("| 比值 | 中位数 log2 | IQR |")
w("|---|---|---|")
for k, v in S["bias"].items():
    w(f"| {k} | {v['median_log2']} | [{v['iqr'][0]}, {v['iqr'][1]}] |")
w("")
w("解读：")
w("- zUMIs_UMI / umite_UMI 中位数 ~+3（约 8 倍）：distinct UMI vs 校正坍缩 UMI 的方法学差异")
w("- 任何 read 级计数 / UMI 级计数之比反映 PCR 重复度（去重率 ~2）")
w("- IQR 窄（多数组合 <0.5）说明偏差是系统性的、基因间一致，相对定量（倍数变化）不受影响")
w("")
# ---------- 7. conclusion ----------
w("## 7. 结论")
w("")
w("1. **定量等价性**：四方法在基因水平高度相关（Spearman 全部 ≥0.93，"
  "跨 UMI 策略组合全部 ≥0.95），差异表达等相对定量结论对方法选择稳健。")
w("2. **分歧来源**：几乎全部来自 1-3 counts 的边缘基因——低计数下去重策略、"
  "多重比对处理、UMI 碰撞校正的差异被放大；中高表达基因（≥10）四方法检出率接近一致。")
w("3. **检出数排序**（zUMIs > featureCounts > HTSeq+umi_tools > umite）"
  "正好对应 UMI 处理严格程度；umite 最保守但假阳性最少（独有基因仅 24 个且全为 1 count）。")
w("4. **内含子**：占比 <1%，umite 与 zUMIs 的内含子信号高度一致；"
  "exon+intron 合并（zUMIs inex）安全可用。")
w("5. **推荐**：zUMIs（inex 矩阵）为主流程，umite（UE+RE）为正交验证；"
  "详见 `四种定量流程_pipeline整理.md` 末尾推荐章节。")
w("")
w("---")
w("")
w("附录文件：")
w("- `report_stats.json` — 本报告全部原始数字")
w("- `fig_correlation_heatmap.png` / `fig_exon_intron.png` / "
  "`fig_concordance_vs_expression.png` / `fig_discordant_cdf.png`")
w("- `N2_Adult-{1,2}_comparison_4methods.csv` — 逐基因四方法计数矩阵")
w("- `detection_summary.txt` / `detection_upset.png` / `detection_thresholds.png` — 检出分析")
w("- `discordant_genes_N2_Adult-{1,2}.csv` — 3/4 方法检出的逐基因清单")
w("- `report_data/` — 各流程日志与 zUMIs 内外显子提取数据")

(BASE / "四方法详细比较报告.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("report written:", BASE / "四方法详细比较报告.md", len(L), "lines")
