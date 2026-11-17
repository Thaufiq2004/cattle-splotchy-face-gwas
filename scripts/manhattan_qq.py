#!/usr/bin/env python3
"""Manhattan and QQ plots (and genomic inflation lambda) for the primary Firth GWAS.

Run from the repository root:  python scripts/manhattan_qq.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

F = "results/association/splotchy_gwas_firth_forced.Splotchy.glm.firth"
d = pd.read_csv(F, sep=r"\s+")
d = d[(d.TEST == "ADD") & d.P.notna() & (d.P > 0)].copy()
d["CHR"] = pd.to_numeric(d["#CHROM"], errors="coerce")
d = d.dropna(subset=["CHR"]).sort_values(["CHR", "POS"])
d["lp"] = -np.log10(d.P)

# genomic inflation factor from squared Wald Z statistics
chi2 = d.Z_STAT.astype(float) ** 2
lam = np.nanmedian(chi2) / stats.chi2.ppf(0.5, 1)
print("SNPs:", len(d), " lambda:", round(lam, 3))

# Manhattan plot
off, ticks, x = 0, [], []
for c, g in d.groupby("CHR"):
    x.append(g.POS + off)
    ticks.append((c, off + g.POS.max() / 2))
    off += g.POS.max()
d["x"] = pd.concat(x)
fig, ax = plt.subplots(figsize=(13, 5))
cols = ["#1f77b4", "#7f7f7f"]
for i, (c, g) in enumerate(d.groupby("CHR")):
    ax.scatter(g.x, g.lp, s=6, color=cols[i % 2])
ax.axhline(-np.log10(0.05 / len(d)), ls="--", color="red")
ax.set_xticks([t[1] for t in ticks])
ax.set_xticklabels([int(t[0]) for t in ticks], fontsize=7)
ax.set_xlabel("Chromosome")
ax.set_ylabel("-log10(P)")
ax.set_title("Splotchy face GWAS (Firth, PC1-3, n=122)")
plt.tight_layout()
plt.savefig("results/plots/splotchy_manhattan.png", dpi=200)

# QQ plot
o = np.sort(d.lp.values)[::-1]
e = -np.log10((np.arange(1, len(o) + 1) - 0.5) / len(o))
plt.figure(figsize=(5, 5))
plt.scatter(e, o, s=6)
plt.plot([0, e.max()], [0, e.max()], "r--")
plt.xlabel("Expected -log10(P)")
plt.ylabel("Observed -log10(P)")
plt.title(f"QQ plot (lambda = {lam:.2f})")
plt.tight_layout()
plt.savefig("results/plots/splotchy_qq.png", dpi=200)
print("Saved results/plots/splotchy_manhattan.png and results/plots/splotchy_qq.png")
