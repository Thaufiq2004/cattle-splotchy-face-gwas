#!/usr/bin/env bash
# Full analysis pipeline for the splotchy-face GWAS (PLINK 2 + Python).
# Run from the repository root:   bash scripts/run_pipeline.sh
#
# Requirements
#   * PLINK v2.0.0-a.7.12 (or newer) on your PATH, or set PLINK2=/path/to/plink2
#   * Python 3 with the packages in requirements.txt
#   * data/raw/splotchy_GWAS_geno_set.{bed,bim,fam} and data/raw/splotchy_phenos.txt (see data/raw/README.md)
#   * data/reference/ARS-UCD1.2.ncbiRefSeq.gtf.gz              (see data/reference/README.md)
#
# The PLINK options below are the ones recorded in the original run logs (results/logs/).
# Set RUN_SENSITIVITY=1 to also re-run the logistic sensitivity analyses of report section 3.3.
set -euo pipefail

PLINK2="${PLINK2:-plink2}"
PY="${PYTHON:-python3}"
CHR="--chr-set 29"                      # cattle: 29 autosomes
B=data/processed
Q=results/qc
A=results/association
mkdir -p "$B" "$Q" "$A" results/plots

# ---- 1. Quality control (sample missingness, variant missingness, MAF) -----------------------
"$PLINK2" --bfile data/raw/splotchy_GWAS_geno_set $CHR \
  --mind 0.05 --geno 0.05 --maf 0.01 --make-bed --out $B/splotchy_qc

# ---- 2. Unique variant IDs (364 IDs are duplicated in the raw data) and PLINK phenotype file --
$PY scripts/prepare_inputs.py unique-ids $B/splotchy_qc $B/splotchy_qc_unique
$PY scripts/prepare_inputs.py phenotype data/raw/splotchy_phenos.txt $B/splotchy_pheno_plink.txt
$PY scripts/prepare_inputs.py remove-list $Q/remove_15_23.txt     # duplicated genotype, discordant phenotype

# ---- 3. Duplicate / relatedness detection (all 124 animals) ----------------------------------
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --indep-pairwise 50 5 0.2 --out $Q/splotchy_prune
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --extract $Q/splotchy_prune.prune.in \
  --make-king-table --out $Q/splotchy_king      # animals 15 and 23: kinship ~0.5, IBS0 = 0

# ---- 4. Population structure after removing animals 15 and 23 (n = 122) ----------------------
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
  --indep-pairwise 50 5 0.2 --out $Q/splotchy_prune_no15_23
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
  --extract $Q/splotchy_prune_no15_23.prune.in --pca 10 --out $Q/splotchy_pca_no15_23
$PY scripts/prepare_inputs.py pcs $Q/splotchy_pca_no15_23.eigenvec $B/splotchy_pcs_no15_23.txt --n 3

# ---- 5. PRIMARY ASSOCIATION: Firth logistic regression, additive model, PC1-PC3 --------------
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
  --pheno $B/splotchy_pheno_plink.txt --pheno-name Splotchy \
  --covar $B/splotchy_pcs_no15_23.txt --covar-name PC1,PC2,PC3 \
  --glm firth hide-covar --out $A/splotchy_gwas_firth_forced

# ---- 6. Clumping of significant signals (P1 = 0.05 / 35,048) ---------------------------------
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
  --clump $A/splotchy_gwas_firth_forced.Splotchy.glm.firth \
  --clump-p1 1.4266e-6 --clump-p2 0.001 --clump-r2 0.2 --clump-kb 1000 --out $A/splotchy_clump

# ---- 7. Robustness: relatives removed (KING cutoff 0.177), PCA and GWAS recomputed ------------
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
  --extract $Q/splotchy_prune.prune.in --king-cutoff 0.177 --out $Q/splotchy_unrel
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR \
  --extract $Q/splotchy_prune.prune.in --keep $Q/splotchy_unrel.king.cutoff.in.id \
  --pca 5 --out $Q/splotchy_unrel_pca
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --keep $Q/splotchy_unrel.king.cutoff.in.id \
  --pheno $B/splotchy_pheno_plink.txt --pheno-name Splotchy \
  --covar $Q/splotchy_unrel_pca.eigenvec --covar-name PC1,PC2,PC3 \
  --glm firth hide-covar --out $A/splotchy_unrel

# ---- 8. Genotypes of the candidate MITF SNPs (22:31,651,379 and 22:31,651,404) ---------------
printf '%s\n' 22:31651379:G:A:28851 22:31651404:G:A:28852 > $Q/mitf_candidate_ids.txt
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
  --extract $Q/mitf_candidate_ids.txt --export A --out $Q/mitf_candidate_genotypes
"$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
  --ld 22:31651379:G:A:28851 22:31651404:G:A:28852 --out $Q/mitf_candidate_ld

# ---- 9. Optional logistic sensitivity analyses (report section 3.3) ---------------------------
if [[ "${RUN_SENSITIVITY:-0}" == "1" ]]; then
  "$PLINK2" --bfile $B/splotchy_qc $CHR --pheno $B/splotchy_pheno_plink.txt --pheno-name Splotchy \
    --glm hide-covar allow-no-covars --out $A/splotchy_gwas                       # no covariates, n = 124
  "$PLINK2" --bfile $B/splotchy_qc_unique $CHR --extract $Q/splotchy_prune.prune.in --pca 10 \
    --out $Q/splotchy_pca
  $PY scripts/prepare_inputs.py pcs $Q/splotchy_pca.eigenvec $B/splotchy_pcs.txt --n 3
  "$PLINK2" --bfile $B/splotchy_qc_unique $CHR --pheno $B/splotchy_pheno_plink.txt --pheno-name Splotchy \
    --covar $B/splotchy_pcs.txt --covar-name PC1,PC2,PC3 --glm hide-covar \
    --out $A/splotchy_gwas_pca3                                                   # PC1-3, n = 124
  "$PLINK2" --bfile $B/splotchy_qc_unique $CHR --remove $Q/remove_15_23.txt \
    --pheno $B/splotchy_pheno_plink.txt --pheno-name Splotchy \
    --covar $B/splotchy_pcs_no15_23.txt --covar-name PC1,PC2,PC3 --glm hide-covar \
    --out $A/splotchy_gwas_no15_23_pca3                                           # PC1-3, n = 122
fi

# ---- 10. Annotation, figures and a check of the numbers quoted in the report -----------------
$PY scripts/annotate_variants.py
$PY scripts/manhattan_qq.py
$PY scripts/chr22_plot.py
$PY scripts/regional_plot.py
$PY scripts/verify_report_numbers.py
echo "Done. Key outputs are in results/."
