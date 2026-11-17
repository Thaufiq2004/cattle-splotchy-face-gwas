#!/usr/bin/env python3
"""Chromosome 22 association plot (report Figure 3). Run from the repository root."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

F = "results/association/splotchy_gwas_firth_forced.Splotchy.glm.firth"
d = pd.read_csv(F, sep=r"\s+")
d = d[(d.TEST == "ADD") & d.P.notna() & (d.P > 0)].copy()
n_tests = len(d)                                   # genome-wide number of tests for Bonferroni
c = d[d["#CHROM"] == 22].copy()
c["lp"] = -np.log10(c.P)

plt.figure(figsize=(12, 5))
plt.scatter(c.POS / 1e6, c.lp, s=8)
plt.axhline(-np.log10(0.05 / n_tests), ls="--", color="red")
plt.xlabel("Chr22 position (Mb)")
plt.ylabel("-log10(P)")
plt.title("Splotchy face GWAS, chromosome 22")
plt.tight_layout()
plt.savefig("results/plots/splotchy_chr22_fixed.png", dpi=200)
print("chr22 SNPs plotted:", len(c))
