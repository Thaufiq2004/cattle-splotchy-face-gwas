#!/usr/bin/env python3
"""Recompute the headline numbers quoted in the analysis report directly from the PLINK output.

Run from the repository root after the pipeline:  python scripts/verify_report_numbers.py
Writes results/report_numbers.txt (also printed). Sensitivity analyses whose result files are
missing are skipped, so the script also works on the trimmed result set stored in this repository.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

A = Path("results/association")
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def load(path):
    d = pd.read_csv(path, sep=r"\s+")
    d = d[d.TEST == "ADD"].copy()
    d["chi2"] = d.Z_STAT.astype(float) ** 2
    return d


def lam(chi2):
    return float(np.nanmedian(chi2) / stats.chi2.ppf(0.5, 1))


# ---- primary analysis: Firth, PC1-PC3, n = 122 -------------------------------------------
p = load(A / "splotchy_gwas_firth_forced.Splotchy.glm.firth")
n = int(p.P.notna().sum())
thr = 0.05 / n
sig = p[p.P < thr]
say("PRIMARY ANALYSIS (Firth, PC1-PC3)")
say(f"  SNPs tested                     : {n}")
say(f"  Bonferroni threshold            : {thr:.3g}  (-log10 = {-np.log10(thr):.2f})")
say(f"  Significant SNPs                : {len(sig)}  on chr {[int(c) for c in sorted(sig['#CHROM'].unique())]}, "
    f"{sig.POS.min():,}-{sig.POS.max():,}")
top = p.loc[p.P.idxmin()]
say(f"  Lowest P                        : {top.P:.2g} at chr{int(top['#CHROM'])}:{int(top.POS):,}")
say(f"  Best -log10P outside chr22      : {(-np.log10(p[p['#CHROM'] != 22].P)).max():.2f}")
say(f"  OR range, significant SNPs      : {sig.OR.min():.1f} - {sig.OR.max():.1f}")
l_all, l_ex22 = lam(p.chi2), lam(p[p["#CHROM"] != 22].chi2)
say(f"  lambda (all / excluding chr22)  : {l_all:.3f} / {l_ex22:.3f}")
gc_p = stats.chi2.sf(p.chi2 / l_ex22, 1)
say(f"  After genomic control           : {(gc_p < thr).sum()} significant, lowest P {gc_p.min():.2g}")
for pos in (31651379, 31651404, 31655675):
    r = p[p.POS == pos]
    if len(r):
        say(f"  chr22:{pos:,}                P = {r.P.iloc[0]:.2g}, OR = {r.OR.iloc[0]:.1f}")

# ---- unrelated subset ----------------------------------------------------------------------
f = A / "splotchy_unrel.Splotchy.glm.firth"
if f.exists():
    u = load(f)
    say("\nRELATIVES REMOVED (KING cutoff 0.177)")
    say(f"  lambda                          : {lam(u.chi2):.3f}")
    say(f"  Significant SNPs                : {(u.P < thr).sum()}, lowest P {u.P.min():.2g}")

# ---- optional sensitivity analyses ---------------------------------------------------------
sens = [("Logistic, no covariates (124 animals)", "splotchy_gwas.Splotchy.glm.logistic.hybrid"),
        ("Logistic + PC1-PC3 (124 animals)", "splotchy_gwas_pca3.Splotchy.glm.logistic.hybrid"),
        ("Duplicate removed, PC1-PC3 (122)", "splotchy_gwas_no15_23_pca3.Splotchy.glm.logistic.hybrid")]
rows = []
for label, fn in sens:
    if (A / fn).exists():
        d = load(A / fn)
        d = d[d.P.notna()]
        t = d.loc[d.P.idxmin()]
        g = {pos: d[d.POS == pos].P.iloc[0] for pos in (31651379, 31651404)}
        rows.append(f"  {label:<40} top {int(t.POS):,} P={t.P:.2g}; "
                    f"31,651,379 P={g[31651379]:.2g}; 31,651,404 P={g[31651404]:.2g}")
if rows:
    say("\nSENSITIVITY ANALYSES")
    for r in rows:
        say(r)

# ---- genotype-by-phenotype tables at the two lead SNPs -------------------------------------
raw = Path("results/qc/mitf_candidate_genotypes.raw")
pheno = Path("data/raw/splotchy_phenos.txt")
if raw.exists() and pheno.exists():
    g = pd.read_csv(raw, sep=r"\s+").merge(pd.read_csv(pheno, sep=r"\s+"), on=["FID", "IID"])
    say("\nGENOTYPE TABLES (counted allele = A; 1 = heterozygous, 2 = homozygous A)")
    for c in [c for c in g.columns if c.startswith("22:")]:
        t = pd.crosstab(g.Splotchy.map({1: "case", 0: "control"}), g[c])
        say(f"  {c}")
        say("    " + t.to_string().replace("\n", "\n    "))

Path("results/report_numbers.txt").write_text("\n".join(lines) + "\n")
