#!/usr/bin/env python3
"""Annotate genome-wide significant SNPs against the NCBI RefSeq GTF (ARS-UCD1.2).

For every significant SNP the script reports
  * class   : exonic / intronic / upstream (<=5 kb) / downstream (<=5 kb) / intergenic
  * genes   : genes whose transcript/exon records overlap the SNP (or the nearest transcript)
  * nearest_TSS_gene, nearest_TSS_tx, bp_upstream_of_TSS : signed distance to the nearest
    transcription start site (positive = upstream of the TSS, strand-aware)

Transcript/exon records are used rather than whole-gene spans, because repeated tRNA genes
share one gene_id and give meaningless spans.

Run from the repository root:  python scripts/annotate_variants.py
"""
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

FIRTH_FILE = "results/association/splotchy_gwas_firth_forced.Splotchy.glm.firth"
SNP_FILE = "results/association/significant_snps.tsv"
GTF = "data/reference/ARS-UCD1.2.ncbiRefSeq.gtf.gz"
OUT = "results/association/significant_snps_annotated_v3.tsv"
CHR22_ACCESSION = "NC_037349.1"   # RefSeq accession of ARS-UCD1.2 chromosome 22
ALPHA = 0.05

if not Path(SNP_FILE).exists():     # build the significant-SNP table from the primary GWAS
    d = pd.read_csv(FIRTH_FILE, sep=r"\s+")
    d = d[(d.TEST == "ADD") & d.P.notna()]
    thr = ALPHA / len(d)
    sig = d[d.P < thr].rename(columns={"#CHROM": "CHR"})[["CHR", "POS", "ID", "P"]]
    sig.to_csv(SNP_FILE, sep="\t", index=False)
    print(f"Bonferroni threshold {thr:.3g}; {len(sig)} significant SNPs -> {SNP_FILE}")

snps = pd.read_csv(SNP_FILE, sep="\t")
tx, ex = [], []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line[0] == "#":
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] != CHR22_ACCESSION or p[2] not in ("transcript", "exon"):
            continue
        g = re.search(r'gene_name "([^"]+)"', p[8])
        t = re.search(r'transcript_id "([^"]+)"', p[8])
        rec = (int(p[3]), int(p[4]), p[6], g.group(1) if g else ".", t.group(1) if t else ".")
        (tx if p[2] == "transcript" else ex).append(rec)

rows = []
for _, s in snps.iterrows():
    pos = int(s["POS"])
    oex = [x for x in ex if x[0] <= pos <= x[1]]
    otx = [x for x in tx if x[0] <= pos <= x[1]]
    if oex:
        cls, genes = "exonic", {x[3] for x in oex}
    elif otx:
        cls, genes = "intronic", {x[3] for x in otx}
    else:
        def dist(x):
            return x[0] - pos if pos < x[0] else pos - x[1]
        n = min(tx, key=dist)
        before = pos < n[0]
        up = (before and n[2] == "+") or (not before and n[2] == "-")
        cls = ("upstream" if up else "downstream") if dist(n) <= 5000 else "intergenic"
        genes = {n[3]}

    def tss(x):                       # transcription start site of a transcript
        return x[0] if x[2] == "+" else x[1]

    def signed(x):                    # positive = SNP lies upstream of the TSS
        return tss(x) - pos if x[2] == "+" else pos - tss(x)

    t = min(tx, key=lambda x: abs(signed(x)))
    rows.append({**s.to_dict(), "class": cls, "genes": ";".join(sorted(genes)),
                 "nearest_TSS_gene": t[3], "nearest_TSS_tx": t[4], "bp_upstream_of_TSS": signed(t)})

out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False)
print(out["class"].value_counts().to_string())
print(f"wrote {OUT}")
