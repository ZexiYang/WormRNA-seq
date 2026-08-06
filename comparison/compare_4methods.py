#!/usr/bin/env python3
"""4-method comparison for N2_Adult-1/2 Smart-seq3:
featureCounts (reads), zUMIs (reads/UMI), HTSeq+umi_tools (UMI), umite (reads/UMI).
umite gene IDs are NCBI locus_tag (CELE_*); bridge to WBGene via NCBI gene_info.
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
GENE_INFO = Path.home() / "rnaseq_project/singleworm_260727/ref/cel_gene_info.gz"
SAMPLES = ["N2_Adult-1", "N2_Adult-2"]

# ---------- CELE_ -> WBGene bridge ----------
cele2wb = {}
with gzip.open(GENE_INFO, "rt") as f:
    header = None
    for line in f:
        if line.startswith("#"):
            header = line.lstrip("#").rstrip("\n").split("\t")
            i_loc = header.index("LocusTag")
            i_xref = header.index("dbXrefs")
            continue
        if header is None:
            continue
        parts = line.rstrip("\n").split("\t")
        loc = parts[i_loc]
        if not loc or loc == "-":
            continue
        m = re.search(r"WormBase:(WBGene\d+)", parts[i_xref])
        if m:
            cele2wb[loc] = m.group(1)
print(f"bridge: {len(cele2wb)} CELE -> WBGene mappings")

# ---------- load umite tsvs ----------
def load_umite(suffix):
    df = pd.read_csv(BASE / f"umite.{suffix}.tsv", sep="\t", index_col=0)
    df.index = [s.replace(".namesorted.bam", "") for s in df.index]
    drop = [c for c in df.columns if c.startswith("_")]
    return df.drop(columns=drop)

UE = load_umite("UE")   # UMI exon (dedup + corrected)
RE = load_umite("RE")   # internal-read exon
D = load_umite("D")     # UMI-duplicate exon reads
umite_umi = UE
umite_reads = UE + RE + D

def to_wb(row):
    """map CELE columns of one sample row to WBGene, summed"""
    s = row.rename(index=lambda c: cele2wb.get(c))
    s = s[s.index.notna()]
    return s.groupby(level=0).sum()

summary_lines = []
for sample in SAMPLES:
    comp = pd.read_csv(BASE / f"{sample}_comparison_wb.csv").set_index("gene")
    comp.columns = [c if not c.startswith("zUMIs") else c for c in comp.columns]

    u_umi = to_wb(umite_umi.loc[sample])
    u_reads = to_wb(umite_reads.loc[sample])
    n_mapped = len(u_umi)

    df = comp.copy()
    df["umite_UMI"] = u_umi
    df["umite_reads"] = u_reads
    df = df.fillna(0).astype(int)
    df.to_csv(BASE / f"{sample}_comparison_4methods.csv")

    # stats
    def corr(a, b, how="spearman"):
        x, y = np.log10(df[a] + 1), np.log10(df[b] + 1)
        if how == "pearson":
            return stats.pearsonr(x, y)[0]
        return stats.spearmanr(df[a], df[b])[0]

    L = []
    L.append(f"sample: {sample}")
    L.append(f"CELE genes mapped to WBGene : {n_mapped} / {umite_umi.shape[1]}")
    L.append(f"totals: fc_reads={df['featureCounts_reads'].sum():,}  "
             f"zUMIs_UMI={df['zUMIs_UMI'].sum():,}  htseq_UMI={df['htseq_UMI'].sum():,}  "
             f"umite_UMI={df['umite_UMI'].sum():,}  umite_reads={df['umite_reads'].sum():,}")
    det = {c: int((df[c] > 0).sum()) for c in
           ["featureCounts_reads", "zUMIs_UMI", "htseq_UMI", "umite_UMI"]}
    L.append("detected genes: " + "  ".join(f"{k}={v:,}" for k, v in det.items()))
    pairs = [
        ("umite_UMI", "zUMIs_UMI"),
        ("umite_UMI", "htseq_UMI"),
        ("umite_UMI", "featureCounts_reads"),
        ("umite_reads", "featureCounts_reads"),
        ("umite_reads", "zUMIs_reads"),
        ("zUMIs_UMI", "htseq_UMI"),
    ]
    for a, b in pairs:
        L.append(f"  {a:15s} vs {b:20s} Pearson(log10+1)={corr(a,b,'pearson'):.4f}  "
                 f"Spearman={corr(a,b,'spearman'):.4f}")
    L.append("")
    summary_lines.extend(L)

    # scatter plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, (a, b) in zip(axes.flat, pairs):
        x, y = np.log10(df[a] + 1), np.log10(df[b] + 1)
        ax.hexbin(x, y, gridsize=60, cmap="viridis", mincnt=1, bins="log")
        ax.set_xlabel(a)
        ax.set_ylabel(b)
        ax.set_title(f"rho={stats.spearmanr(df[a], df[b])[0]:.3f}", fontsize=9)
    fig.suptitle(f"{sample}: umite vs previous methods (log10(count+1))")
    fig.tight_layout()
    fig.savefig(BASE / f"{sample}_umite_vs_3methods.png", dpi=120)
    plt.close(fig)

summary = "\n".join(summary_lines)
(BASE / "comparison_4methods_summary.txt").write_text(summary + "\n")
print(summary)
