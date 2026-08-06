#!/usr/bin/env python3
"""Comprehensive 4-method comparison report for N2_Adult-1/2 Smart-seq3.

Methods: 公司 bulk (fastp+STAR+featureCounts), 自建 umi_tools+HTSeq, zUMIs, umite.
Inputs (all in report_data/ or parent dir):
  ../N2_Adult-{1,2}_comparison_4methods.csv   per-gene counts (WBGene)
  ../umite.{UE,UI,RE,RI,D}.tsv                umite exon/intron x UMI/internal (CELE_)
  zumis_exin_*.csv                            zUMIs exon/intron per-gene (WBGene)
  *.log / *.txt                               pipeline QC stats
Outputs: 四方法详细比较报告.md, fig_correlation_heatmap.png, fig_exon_intron.png,
         fig_concordance_vs_expression.png, fig_discordant_cdf.png
"""
import gzip
import json
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
RD = BASE / "report_data"
GENE_INFO = Path.home() / "rnaseq_project/singleworm_260727/ref/cel_gene_info.gz"
SAMPLES = ["N2_Adult-1", "N2_Adult-2"]

R = {"samples": {}}          # results dict, rendered into markdown at the end

# ============ data loading ============
dfs = {s: pd.read_csv(BASE / f"{s}_comparison_4methods.csv").set_index("gene")
       for s in SAMPLES}
pooled = sum(dfs[s] for s in SAMPLES)   # gene x count-type, summed over 2 samples

# CELE -> WBGene bridge
cele2wb = {}
with gzip.open(GENE_INFO, "rt") as f:
    header = None
    for line in f:
        if line.startswith("#"):
            header = line.lstrip("#").rstrip("\n").split("\t")
            i_loc, i_xref = header.index("LocusTag"), header.index("dbXrefs")
            continue
        p = line.rstrip("\n").split("\t")
        m = re.search(r"WormBase:(WBGene\d+)", p[i_xref])
        if m and p[i_loc] not in ("", "-"):
            cele2wb[p[i_loc]] = m.group(1)

def load_umite(suffix):
    df = pd.read_csv(BASE / f"umite.{suffix}.tsv", sep="\t", index_col=0)
    df.index = [s.replace(".namesorted.bam", "") for s in df.index]
    return df.drop(columns=[c for c in df.columns if c.startswith("_")])

UM = {k: load_umite(k) for k in ["UE", "UI", "RE", "RI", "D"]}

def umite_wb(df, sample):
    s = df.loc[sample].rename(index=lambda c: cele2wb.get(c))
    return s[s.index.notna()].groupby(level=0).sum()

# zUMIs exon/intron
zexin = {}
for s in SAMPLES:
    z = pd.read_csv(RD / f"zumis_exin_{s}.csv").set_index("gene")
    z.columns = ["zumis_umi_exon", "zumis_umi_intron", "zumis_read_exon",
                 "zumis_read_intron", "zumis_int_read_exon", "zumis_int_read_intron"]
    zexin[s] = z

# HTSeq raw/umi per gene
htseq = {}
for s in SAMPLES:
    h = {}
    for kind in ["raw", "umi"]:
        t = pd.read_csv(RD / f"{s}.htseq_{kind}.tsv", sep="\t", header=None,
                        names=["gene", "count"], index_col=0)["count"]
        t = t[~t.index.str.startswith("__")]
        h[kind] = t
    htseq[s] = h

# ============ 1. pipeline stats from logs ============
def star_stats(path):
    d = {}
    for line in open(path):
        if "|" not in line:
            continue
        k, v = [x.strip() for x in line.split("|", 1)]
        d[k] = v
    return d

def pct(x):
    return float(x.rstrip("%"))

pipe = {}
for s in SAMPLES:
    ue_log = (RD / "umiextract.log").read_text()
    m = re.search(rf"{s}_1\.fq\.gz: (\d+) reads, (\d+) with UMI \(([\d.]+)%\)", ue_log)
    n_in, n_umi, p_umi = int(m.group(1)), int(m.group(2)), float(m.group(3))

    st_u = star_stats(RD / f"umite_STAR_{s}.log")
    st_z = star_stats(RD / f"zumis_STAR_{s}.log")
    st_t = star_stats(RD / f"umitools_STAR_{s}.log")

    fc = pd.read_csv(RD / f"fc_summary_{s}.txt", sep="\t", index_col=0).iloc[:, 0]
    n_frag = int(fc.sum())

    # company pipeline input fragments = assigned + unassigned (pairs)
    pipe[s] = {
        "raw_reads": n_in,
        "umi_detected": n_umi, "umi_pct": p_umi,
        "umite_star": st_u, "zumis_star": st_z, "umitools_star": st_t,
        "fc_assigned": int(fc["Assigned"]),
        "fc_multi": int(fc["Unassigned_MultiMapping"]),
        "fc_ambi": int(fc["Unassigned_Ambiguity"]),
        "fc_nofeat": int(fc["Unassigned_NoFeatures"]),
    }
R["pipe"] = pipe

# ============ 2. totals & dedup ============
tot = {}
for s in SAMPLES:
    df = dfs[s]
    tot[s] = {
        "featureCounts_reads": int(df["featureCounts_reads"].sum()),
        "zUMIs_reads": int(df["zUMIs_reads"].sum()),
        "zUMIs_UMI": int(df["zUMIs_UMI"].sum()),
        "htseq_raw": int(htseq[s]["raw"].sum()),
        "htseq_UMI": int(htseq[s]["umi"].sum()),
        "umite_UMI": int(df["umite_UMI"].sum()),
        "umite_reads": int(df["umite_reads"].sum()),
        "umite_UE": int(UM["UE"].loc[s].sum()), "umite_UI": int(UM["UI"].loc[s].sum()),
        "umite_RE": int(UM["RE"].loc[s].sum()), "umite_RI": int(UM["RI"].loc[s].sum()),
        "umite_D": int(UM["D"].loc[s].sum()),
    }
R["totals"] = tot

# ============ 3. exon / intron ============
exin = {}
for s in SAMPLES:
    z = zexin[s]
    ue, ui, re_, ri = (UM[k].loc[s].sum() for k in ["UE", "UI", "RE", "RI"])
    exin[s] = {
        "umite_umi_exon": int(ue), "umite_umi_intron": int(ui),
        "umite_read_exon": int(re_ + ue + UM["D"].loc[s].sum()),
        "umite_int_read_exon": int(re_), "umite_int_read_intron": int(ri),
        "zumis_umi_exon": int(z["zumis_umi_exon"].sum()),
        "zumis_umi_intron": int(z["zumis_umi_intron"].sum()),
        "zumis_read_exon": int(z["zumis_read_exon"].sum()),
        "zumis_read_intron": int(z["zumis_read_intron"].sum()),
        "zumis_int_read_exon": int(z["zumis_int_read_exon"].sum()),
        "zumis_int_read_intron": int(z["zumis_int_read_intron"].sum()),
    }

# gene-level intron detection (pooled)
umite_ui_wb = umite_wb(UM["UI"] + UM["RI"], SAMPLES[0]).add(
    umite_wb(UM["UI"] + UM["RI"], SAMPLES[1]), fill_value=0)
umite_ue_wb = umite_wb(UM["UE"] + UM["RE"] + UM["D"], SAMPLES[0]).add(
    umite_wb(UM["UE"] + UM["RE"] + UM["D"], SAMPLES[1]), fill_value=0)
z_intron = zexin[SAMPLES[0]][["zumis_umi_intron", "zumis_read_intron"]].add(
    zexin[SAMPLES[1]][["zumis_umi_intron", "zumis_read_intron"]], fill_value=0)
z_exon = zexin[SAMPLES[0]][["zumis_umi_exon", "zumis_read_exon"]].add(
    zexin[SAMPLES[1]][["zumis_umi_exon", "zumis_read_exon"]], fill_value=0)

u_in_genes = set(umite_ui_wb[umite_ui_wb > 0].index)
z_in_genes = set(z_intron[z_intron["zumis_read_intron"] > 0].index)
u_in_only = set(umite_ui_wb[(umite_ui_wb > 0)].index) - set(umite_ue_wb[umite_ue_wb > 0].index)
z_in_only = set(z_intron[z_intron["zumis_read_intron"] > 0].index) - \
            set(z_exon[z_exon["zumis_read_exon"] > 0].index)

# intron correlation umite vs zumis (common genes)
common_in = sorted(u_in_genes & z_in_genes)
ci_x = np.log10(umite_ui_wb[common_in] + 1)
ci_y = np.log10((z_intron["zumis_umi_intron"] + z_intron["zumis_read_intron"])[common_in] + 1)
exin_corr = stats.spearmanr(umite_ui_wb[common_in],
                            (z_intron["zumis_umi_intron"] + z_intron["zumis_read_intron"])[common_in])[0]

R["exin"] = {"per_sample": exin,
             "umite_intron_genes": len(u_in_genes), "zumis_intron_genes": len(z_in_genes),
             "intron_overlap": len(u_in_genes & z_in_genes),
             "umite_intron_only_genes": len(u_in_only), "zumis_intron_only_genes": len(z_in_only),
             "intron_corr": float(exin_corr)}

# ============ 4. correlation matrix ============
CTYPES = ["featureCounts_reads", "zUMIs_reads", "zUMIs_UMI",
          "htseq_raw", "htseq_UMI", "umite_UMI", "umite_reads"]
mat = pooled.copy()
mat["htseq_raw"] = (htseq[SAMPLES[0]]["raw"].add(htseq[SAMPLES[1]]["raw"], fill_value=0)
                    .reindex(mat.index).fillna(0).astype(int))
mat["htseq_UMI"] = (htseq[SAMPLES[0]]["umi"].add(htseq[SAMPLES[1]]["umi"], fill_value=0)
                    .reindex(mat.index).fillna(0).astype(int))
M = mat[CTYPES]

spear = pd.DataFrame(np.eye(len(CTYPES)), index=CTYPES, columns=CTYPES)
pear = spear.copy()
for a, b in combinations(CTYPES, 2):
    spear.loc[a, b] = spear.loc[b, a] = stats.spearmanr(M[a], M[b])[0]
    pear.loc[a, b] = pear.loc[b, a] = stats.pearsonr(np.log10(M[a] + 1),
                                                     np.log10(M[b] + 1))[0]
R["spearman"] = spear.round(4).to_dict()
R["pearson"] = pear.round(4).to_dict()

# per-sample key correlations
key_pairs = [("umite_UMI", "zUMIs_UMI"), ("umite_UMI", "htseq_UMI"),
             ("umite_UMI", "featureCounts_reads"), ("umite_reads", "featureCounts_reads"),
             ("umite_reads", "zUMIs_reads"), ("zUMIs_UMI", "htseq_UMI"),
             ("htseq_UMI", "featureCounts_reads")]
per_sample_corr = {}
for s in SAMPLES:
    d = dfs[s].copy()
    d["htseq_raw"] = htseq[s]["raw"].reindex(d.index).fillna(0).astype(int)
    d["htseq_UMI"] = htseq[s]["umi"].reindex(d.index).fillna(0).astype(int)
    per_sample_corr[s] = {
        f"{a} vs {b}": {"pearson_log": round(stats.pearsonr(np.log10(d[a]+1), np.log10(d[b]+1))[0], 4),
                        "spearman": round(stats.spearmanr(d[a], d[b])[0], 4)}
        for a, b in key_pairs}
R["per_sample_corr"] = per_sample_corr

# ============ 5. divergence concentrated at 1-3 counts ============
MAIN4 = ["featureCounts_reads", "zUMIs_UMI", "htseq_UMI", "umite_UMI"]
det = {m: set(pooled.index[pooled[m] > 0]) for m in MAIN4}
all_genes = set.union(*det.values())
maxcnt = pooled[MAIN4].max(axis=1)

bins = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 10), (11, 50), (51, 10**9)]
bin_labels = ["1", "2", "3", "4-5", "6-10", "11-50", ">50"]
conc = []
for (lo, hi), lab in zip(bins, bin_labels):
    g = [g for g in all_genes if lo <= maxcnt[g] <= hi]
    if not g:
        continue
    ndet = (pooled.loc[g, MAIN4] > 0).sum(axis=1)
    conc.append({"bin": lab, "n": len(g),
                 "all4": int((ndet == 4).sum()), "ge3": int((ndet >= 3).sum()),
                 "pct_all4": round((ndet == 4).mean() * 100, 1),
                 "pct_ge3": round((ndet >= 3).mean() * 100, 1)})
R["concordance_bins"] = conc

# CDF of discordant genes' counts
umite_missed = set.intersection(*[det[m] for m in MAIN4 if m != "umite_UMI"]) - det["umite_UMI"]
missed_max = maxcnt[sorted(umite_missed)]
R["umite_missed"] = {
    "n": len(umite_missed),
    "le3": int((missed_max <= 3).sum()), "pct_le3": round((missed_max <= 3).mean() * 100, 1),
    "le5": int((missed_max <= 5).sum()), "pct_le5": round((missed_max <= 5).mean() * 100, 1),
    "ge10": int((missed_max >= 10).sum()), "ge50": int((missed_max >= 50).sum()),
    "median": float(missed_max.median()), "max": int(missed_max.max())}

uniq_stats = {}
for m in MAIN4:
    u = set(det[m])
    for o in MAIN4:
        if o != m:
            u -= det[o]
    vals = pooled.loc[sorted(u), m] if u else pd.Series(dtype=int)
    uniq_stats[m] = {"n": len(u),
                     "le3": int((vals <= 3).sum()),
                     "pct_le3": round((vals <= 3).mean() * 100, 1) if len(vals) else 0,
                     "median": float(vals.median()) if len(vals) else 0,
                     "max": int(vals.max()) if len(vals) else 0}
R["unique_genes"] = uniq_stats

# all-discordant (any gene not detected by all 4): count distribution
ndet_all = (pooled.loc[sorted(all_genes), MAIN4] > 0).sum(axis=1)
disc = ndet_all[ndet_all < 4].index
disc_max = maxcnt[disc]
R["discordant_overall"] = {
    "n": len(disc), "pct_of_union": round(len(disc) / len(all_genes) * 100, 1),
    "le3": int((disc_max <= 3).sum()),
    "pct_le3": round((disc_max <= 3).mean() * 100, 1),
    "le5_pct": round((disc_max <= 5).mean() * 100, 1),
    "ge10": int((disc_max >= 10).sum())}

# ============ 6. quantitative bias on commonly detected genes ============
common4 = sorted(set.intersection(*det.values()))
C = pooled.loc[common4, MAIN4]
bias = {}
for a, b in combinations(MAIN4, 2):
    r = np.log2((C[a] + 1) / (C[b] + 1))
    bias[f"{a} / {b}"] = {"median_log2": round(float(r.median()), 3),
                          "iqr": [round(float(r.quantile(0.25)), 3),
                                  round(float(r.quantile(0.75)), 3)]}
R["bias"] = bias

# ============ figures ============
# fig 1: correlation heatmaps
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
for ax, m_, title in [(axes[0], spear, "Spearman (raw counts)"),
                      (axes[1], pear, "Pearson (log10(count+1))")]:
    im = ax.imshow(m_.values, vmin=0.85, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(len(CTYPES))); ax.set_yticks(range(len(CTYPES)))
    ax.set_xticklabels(CTYPES, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(CTYPES, fontsize=8)
    ax.set_title(title)
    for i in range(len(CTYPES)):
        for j in range(len(CTYPES)):
            ax.text(j, i, f"{m_.values[i,j]:.3f}", ha="center", va="center",
                    fontsize=7, color="white" if m_.values[i,j] < 0.95 else "black")
fig.suptitle("Pairwise correlation of per-gene counts (pooled 2 samples, WBGene universe)")
fig.colorbar(im, ax=axes, shrink=0.8)
fig.tight_layout()
fig.savefig(BASE / "fig_correlation_heatmap.png", dpi=130)
plt.close(fig)

# fig 2: exon/intron composition
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
for ax, s in zip(axes, SAMPLES):
    cats = ["umite\nUMI", "umite\ninternal", "zUMIs\nUMI", "zUMIs\ninternal", "zUMIs\nall reads"]
    exon_v = [exin[s]["umite_umi_exon"], exin[s]["umite_int_read_exon"],
              exin[s]["zumis_umi_exon"], exin[s]["zumis_int_read_exon"],
              exin[s]["zumis_read_exon"]]
    intron_v = [exin[s]["umite_umi_intron"], exin[s]["umite_int_read_intron"],
                exin[s]["zumis_umi_intron"], exin[s]["zumis_int_read_intron"],
                exin[s]["zumis_read_intron"]]
    x = np.arange(len(cats))
    ax.bar(x, exon_v, label="exon", color="#1f77b4")
    ax.bar(x, intron_v, bottom=exon_v, label="intron", color="#ff7f0e")
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=8)
    ax.set_title(s)
    ax.set_ylabel("counts")
    for i, (e, iv) in enumerate(zip(exon_v, intron_v)):
        frac = iv / (e + iv) * 100 if e + iv else 0
        ax.text(i, e + iv, f"{frac:.1f}% int", ha="center", fontsize=7)
axes[0].legend()
fig.suptitle("Exon / intron composition by method (featureCounts & HTSeq are exon-only by design)")
fig.tight_layout()
fig.savefig(BASE / "fig_exon_intron.png", dpi=130)
plt.close(fig)

# fig 3: concordance vs expression bin
fig, ax = plt.subplots(figsize=(8.5, 4.8))
x = np.arange(len(conc))
ax.bar(x - 0.2, [c["pct_all4"] for c in conc], width=0.4, label="detected by all 4")
ax.bar(x + 0.2, [c["pct_ge3"] for c in conc], width=0.4, label="detected by >=3")
for i, c in enumerate(conc):
    ax.text(i - 0.2, c["pct_all4"] + 1, f'{c["pct_all4"]:.0f}%', ha="center", fontsize=8)
    ax.text(i + 0.2, c["pct_ge3"] + 1, f'{c["pct_ge3"]:.0f}%', ha="center", fontsize=8)
    ax.text(i, -7, f'n={c["n"]:,}', ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([c["bin"] for c in conc])
ax.set_xlabel("gene max count across the 4 methods (pooled)")
ax.set_ylabel("% of genes")
ax.set_ylim(0, 109)
ax.legend()
ax.set_title("Detection concordance rises sharply with expression level")
fig.tight_layout()
fig.savefig(BASE / "fig_concordance_vs_expression.png", dpi=130)
plt.close(fig)

# fig 4: CDF of discordant/unique gene counts
fig, ax = plt.subplots(figsize=(8, 4.8))
for vals, lab, col in [
        (disc_max, f"all discordant genes (n={len(disc_max)})", "#1f77b4"),
        (missed_max, f"missed only by umite (n={len(missed_max)})", "#d62728")]:
    xs = np.sort(vals.values)
    ys = np.arange(1, len(xs) + 1) / len(xs) * 100
    ax.step(xs, ys, where="post", label=lab, color=col)
for m, col in zip(MAIN4, ["#2ca02c", "#9467bd", "#8c564b", "#7f7f7f"]):
    u = set(det[m])
    for o in MAIN4:
        if o != m:
            u -= det[o]
    if u:
        vals = pooled.loc[sorted(u), m].values
        xs = np.sort(vals); ys = np.arange(1, len(xs) + 1) / len(xs) * 100
        ax.step(xs, ys, where="post", ls="--", color=col,
                label=f"unique to {m} (n={len(vals)})")
ax.axvline(3, color="gray", ls=":", lw=1)
ax.text(3.1, 20, "count=3", fontsize=8, color="gray")
ax.set_xscale("log")
ax.set_xlabel("gene count (max across methods for discordant; own count for unique)")
ax.set_ylabel("CDF (%)")
ax.set_title("Discordant & method-unique genes are almost all <=3 counts")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(BASE / "fig_discordant_cdf.png", dpi=130)
plt.close(fig)

# ============ save stats json ============
def default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    raise TypeError

(BASE / "report_stats.json").write_text(json.dumps(R, ensure_ascii=False, indent=1, default=default))
print("stats written; key numbers:")
print(json.dumps({k: R[k] for k in ["concordance_bins", "umite_missed", "discordant_overall", "unique_genes"]},
                 ensure_ascii=False, indent=1, default=default))
