#!/usr/bin/env python3
"""Detection-level comparison of 4 quantification methods on N2_Adult-1/2.

Inputs: N2_Adult-{1,2}_comparison_4methods.csv (per-gene counts, WBGene IDs)
Methods compared (detection = count > 0):
  featureCounts_reads, zUMIs_UMI, htseq_UMI, umite_UMI
Outputs: detection_summary.txt, detection_upset.png, detection_thresholds.png,
         discordant_genes_<sample>.csv
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
SAMPLES = ["N2_Adult-1", "N2_Adult-2"]
METHODS = ["featureCounts_reads", "zUMIs_UMI", "htseq_UMI", "umite_UMI"]
SHORT = {"featureCounts_reads": "featureCounts", "zUMIs_UMI": "zUMIs",
         "htseq_UMI": "HTSeq+umi_tools", "umite_UMI": "umite"}
THRESHOLDS = [1, 2, 3, 5, 10, 20, 50]

dfs = {s: pd.read_csv(BASE / f"{s}_comparison_4methods.csv").set_index("gene")
       for s in SAMPLES}

lines = []

def log(msg=""):
    print(msg)
    lines.append(msg)

# ---------- 1. detection sets & pairwise Jaccard (pooled: detected in >=1 sample) ----------
detect = {}  # method -> set of genes (pooled across samples)
for m in METHODS:
    detect[m] = set().union(*[set(dfs[s].index[dfs[s][m] > 0]) for s in SAMPLES])

log("=" * 70)
log("1. 检出基因数（count>0，两样本并集 / 各样本）")
for m in METHODS:
    per = ", ".join(f"{s}: {(dfs[s][m] > 0).sum()}" for s in SAMPLES)
    log(f"  {SHORT[m]:18s} 并集 {len(detect[m]):6d}   ({per})")

log("")
log("2. 两两 Jaccard 相似度（检出集合，两样本并集）")
header = " " * 20 + "".join(f"{SHORT[m]:>18s}" for m in METHODS)
log(header)
for a in METHODS:
    row = f"{SHORT[a]:20s}"
    for b in METHODS:
        j = len(detect[a] & detect[b]) / len(detect[a] | detect[b])
        row += f"{j:18.4f}"
    log(row)

# ---------- 3. UpSet (pooled) ----------
combos = {}
for r in range(1, 5):
    for combo in combinations(METHODS, r):
        inset = set.intersection(*[detect[m] for m in combo])
        others = set.union(*[detect[m] for m in METHODS if m not in combo]) \
            if len(combo) < 4 else set()
        combos[combo] = inset - others

log("")
log("3. 检出集合交集拆分（UpSet，两样本并集；'独有' = 仅该组合检出）")
for combo, genes in sorted(combos.items(), key=lambda kv: -len(kv[1])):
    label = " ∩ ".join(SHORT[m] for m in combo)
    if len(combo) == 1:
        label = f"仅 {SHORT[combo[0]]} 独有"
    log(f"  {label:55s} {len(genes):6d}")

all4 = set.intersection(*detect.values())
any_ = set.union(*detect.values())
log(f"\n  四方法共同检出: {len(all4)} / 总检出 {len(any_)} "
    f"({len(all4)/len(any_)*100:.1f}%)")

# ---------- 4. method-unique genes: expression level in the detecting method ----------
log("")
log("4. 各方法独有检出基因的计数特征（两样本合并计数）")
merged = sum(dfs[s] for s in SAMPLES)
for m in METHODS:
    uniq = set(detect[m])
    for other in METHODS:
        if other != m:
            uniq -= detect[other]
    if uniq:
        vals = merged.loc[sorted(uniq), m]
        log(f"  {SHORT[m]:18s} 独有 {len(vals):5d} 个: 计数中位数 {vals.median():.0f}, "
            f"均值 {vals.mean():.1f}, ≥10 counts 的 {int((vals >= 10).sum())} 个")

# ---------- 5. genes missed by exactly one method ----------
log("")
log("5. 恰好被一种方法漏检的基因（其他三种都检出）")
for miss in METHODS:
    others = [m for m in METHODS if m != miss]
    genes = set.intersection(*[detect[o] for o in others]) - detect[miss]
    if genes:
        vals = merged.loc[sorted(genes), others].max(axis=1)
        log(f"  仅 {SHORT[miss]:18s} 漏检 {len(genes):5d} 个: "
            f"在其他方法中最大计数中位数 {vals.median():.0f} "
            f"(≥10 的 {int((vals >= 10).sum())} 个)")

# ---------- 6. detection vs count threshold ----------
log("")
log("6. 不同计数阈值下的检出基因数（两样本并集）")
log("  阈值      " + "".join(f"{SHORT[m]:>18s}" for m in METHODS))
for t in THRESHOLDS:
    row = f"  ≥{t:<8d}"
    for m in METHODS:
        n = len(set().union(*[set(dfs[s].index[dfs[s][m] >= t]) for s in SAMPLES]))
        row += f"{n:18d}"
    log(row)

# ---------- plots ----------
# UpSet-style bar plot
fig, ax = plt.subplots(figsize=(12, 5))
items = sorted(combos.items(), key=lambda kv: -len(kv[1]))
labels = [" &\n".join(SHORT[m] for m in c) if len(c) > 1 else f"only\n{SHORT[c[0]]}"
          for c, _ in items]
vals = [len(g) for _, g in items]
colors = ["#d62728" if len(c) == 1 else ("#2ca02c" if len(c) == 4 else "#1f77b4")
          for c, _ in items]
ax.bar(range(len(vals)), vals, color=colors)
ax.set_xticks(range(len(vals)))
ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
ax.set_ylabel("genes")
ax.set_title("Detection overlap of 4 methods (union of 2 samples; red=method-unique, green=all-4)")
for i, v in enumerate(vals):
    ax.text(i, v + 30, str(v), ha="center", fontsize=7)
fig.tight_layout()
fig.savefig(BASE / "detection_upset.png", dpi=130)
plt.close(fig)

# threshold curves (per sample + pooled)
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
for ax, s in zip(axes, SAMPLES + ["pooled"]):
    for m in METHODS:
        ns = []
        for t in THRESHOLDS:
            if s == "pooled":
                n = len(set().union(*[set(dfs[x].index[dfs[x][m] >= t]) for x in SAMPLES]))
            else:
                n = int((dfs[s][m] >= t).sum())
            ns.append(n)
        ax.plot(THRESHOLDS, ns, marker="o", label=SHORT[m])
    ax.set_xscale("log")
    ax.set_title(s)
    ax.set_xlabel("count threshold (log)")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("detected genes")
axes[-1].legend(fontsize=8)
fig.suptitle("Detected genes vs count threshold")
fig.tight_layout()
fig.savefig(BASE / "detection_thresholds.png", dpi=130)
plt.close(fig)

# discordant gene tables (per sample): genes detected by >=3 methods but not all 4
for s in SAMPLES:
    df = dfs[s]
    ndet = (df[METHODS] > 0).sum(axis=1)
    disc = df.loc[ndet == 3, METHODS].sort_values(
        by=METHODS, ascending=False)
    disc["missed_by"] = disc.apply(
        lambda r: ",".join(SHORT[m] for m in METHODS if r[m] == 0), axis=1)
    disc.to_csv(BASE / f"discordant_genes_{s}.csv")
    log(f"\n7. {s}: 被 3/4 方法检出的基因 {len(disc)} 个 -> discordant_genes_{s}.csv")

(BASE / "detection_summary.txt").write_text("\n".join(lines) + "\n")
