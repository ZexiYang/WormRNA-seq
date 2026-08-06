#!/usr/bin/env python3
"""gtf_add_exon_id.py — 把 NCBI 风格 GTF 转成 umite umicount 可用格式。

umicount 硬性要求 exon 行有 exon_id 属性（NCBI GTF 没有），
gene_name 缺失时显示为空。本脚本为每个 exon 补 exon_id（transcript_id:exon_number），
为所有行补 gene_name（取 gene / locus_tag / gene_id）。
用法: python gtf_add_exon_id.py <in.gtf> <out.gtf>
"""
import re
import sys

inp, out = sys.argv[1], sys.argv[2]
n_exon = n_fix = n_gn = 0
with open(inp) as fi, open(out, "w") as fo:
    for line in fi:
        if line.startswith("#"):
            fo.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        attr = fields[8]
        attrs = dict(re.findall(r'(\S+) "([^"]*)";', attr))
        extras = []
        if fields[2] == "exon":
            n_exon += 1
            if "exon_id" not in attrs:
                tid = attrs.get("transcript_id") or attrs.get("gene_id", "NA")
                en = attrs.get("exon_number", "0")
                extras.append('exon_id "%s:%s";' % (tid, en))
                n_fix += 1
        if "gene_name" not in attrs:
            gn = attrs.get("gene") or attrs.get("locus_tag") or attrs.get("gene_id", "")
            extras.append('gene_name "%s";' % gn)
            n_gn += 1
        if extras:
            fields[8] = attr.rstrip() + " " + " ".join(extras)
        fo.write("\t".join(fields) + "\n")
print("exons=%d, exon_id added=%d, gene_name added=%d" % (n_exon, n_fix, n_gn))
