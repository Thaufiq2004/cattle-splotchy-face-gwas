#!/usr/bin/env python3
"""Regional association plot around MITF (chr22: 31.4-32.5 Mb). Run from the repository root."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

F = "results/association/splotchy_gwas_firth_forced.Splotchy.glm.firth"
d = pd.read_csv(F, sep=r"\s+")
d = d[(d["TEST"] == "ADD") & d["P"].notna() & (d["P"] > 0)].copy()
n_tests = len(d)                                   # genome-wide number of tests for Bonferroni

d = d[(d["#CHROM"] == 22) & (d["POS"] >= 31_400_000) & (d["POS"] <= 32_500_000)].copy()
d["-log10P"] = -np.log10(d["P"])

plt.figure(figsize=(12, 5))
plt.scatter(d["POS"] / 1e6, d["-log10P"], s=18, alpha=0.75)
plt.axhline(-np.log10(0.05 / n_tests), linestyle="--", label="Bonferroni threshold")
plt.xlabel("Chromosome 22 position (Mb)")
plt.ylabel("-log10(P)")
plt.title("Regional GWAS: Cattle Splotchy Face")
plt.legend()
plt.tight_layout()
plt.savefig("results/plots/splotchy_regional_association.png", dpi=200)
print("Regional SNPs plotted:", len(d))
