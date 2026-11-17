# Acknowledgements and citations

Citing a source credits and documents where the data, methods and software came from. It does not, by itself,
mean that any organisation or author has reviewed or endorsed this repository.

## 1. Data and original study (please cite whenever the data or results are used)

| Role | Citation |
|---|---|
| Source study | Jivanji S, Wilkinson E, Tang L, Tiplady KM, Yeates A, Harland C, Gray C, Couldrey C, Worth G, Gamache I, Desjardins J, Tabares JAA, Yamanaka N, McNaughton L, Brennan L, Cloutier M-P, Cowan M, Ellison R, Fransen T, Monehan T, Spelman RJ, Snell RG, Charlier C, Yamanaka Y, Garrick D, Mort R, Littlejohn MD. Structural and epistatic regulatory variants cause hallmark white spotting in cattle. *Science Advances* 11(46):eadt5913 (2025). doi:[10.1126/sciadv.adt5913](https://doi.org/10.1126/sciadv.adt5913) (open access, CC BY 4.0) |
| Dataset | Data from: Structural and epistatic regulatory variants cause hallmark white spotting in cattle. Dryad. doi:[10.5061/dryad.fxpnvx15q](https://doi.org/10.5061/dryad.fxpnvx15q) (version of 8 Oct 2025). Check the "Cite" box on the Dryad page for the exact author list and year. |
| Earlier GWAS that first nominated the MITF SNP | Jivanji S, et al. Genome-wide association analysis reveals QTL and candidate mutations involved in white spotting in cattle. *Genet Sel Evol* 51:62 (2019). doi:[10.1186/s12711-019-0506-2](https://doi.org/10.1186/s12711-019-0506-2) |
| Variant knowledge base (independent confirmation of the MITF variant at Chr22:31,651,379) | OMIA - Online Mendelian Inheritance in Animals, University of Sydney. Coat colour, white spotting in *Bos taurus* (OMIA:000214-9913). <https://omia.org/OMIA000214/9913/> (accessed 18 Nov 2026) |

**Organisations and groups named in the source study** (acknowledged for data generation and funding): Livestock Improvement
Corporation (LIC), Hamilton, New Zealand; AL Rae Centre for Genetics and Breeding, Massey University; Unit of Animal
Genomics, GIGA and Faculty of Veterinary Medicine, University of Liège; McGill Integrated Core of Animal Modeling and
Goodman Cancer Institute, McGill University; Lancaster University; University of Liverpool; University of Auckland. The
original work was funded by the New Zealand Ministry of Business, Innovation and Employment Endeavour Fund (programme
grant CONT-57639-ENDRP-LIC) and LIC.

> Transparency note taken from the paper's competing-interests statement: LIC has a patent application related to genetic
> testing of the *MITF* regulatory variant discussed there. This repository is a non-commercial research reanalysis.

## 2. Software and methods

| Used for | Citation |
|---|---|
| PLINK 2 (QC, PCA, KING, GWAS, clumping, LD) | Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 4:7 (2015). doi:10.1186/s13742-015-0047-8 - software: <https://www.cog-genomics.org/plink/2.0/> |
| KING kinship / duplicate detection | Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen W-M. Robust relationship inference in genome-wide association studies. *Bioinformatics* 26(22):2867-2873 (2010). doi:10.1093/bioinformatics/btq559 |
| Firth bias-reduced logistic regression | Firth D. Bias reduction of maximum likelihood estimates. *Biometrika* 80(1):27-38 (1993). doi:10.1093/biomet/80.1.27 |
| Genomic control | Devlin B, Roeder K. Genomic control for association studies. *Biometrics* 55(4):997-1004 (1999). doi:10.1111/j.0006-341X.1999.00997.x |
| Reference genome ARS-UCD1.2 | Rosen BD, et al. De novo assembly of the cattle reference genome with single-molecule sequencing. *GigaScience* 9(3):giaa021 (2020). doi:10.1093/gigascience/giaa021 |
| NCBI RefSeq annotation (GTF release 2025-05-22) | O'Leary NA, et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Res* 44(D1):D733-D745 (2016). doi:10.1093/nar/gkv1189 |
| Ensembl Variant Effect Predictor | McLaren W, et al. The Ensembl Variant Effect Predictor. *Genome Biol* 17:122 (2016). doi:10.1186/s13059-016-0974-4 |
| SciPy / NumPy / pandas / Matplotlib | Virtanen P, et al. *Nat Methods* 17:261-272 (2020), doi:10.1038/s41592-019-0686-2; Harris CR, et al. *Nature* 585:357-362 (2020), doi:10.1038/s41586-020-2649-2; McKinney W. *Proc 9th Python in Science Conf* 56-61 (2010), doi:10.25080/Majora-92bf1922-00a; Hunter JD. *Comput Sci Eng* 9(3):90-95 (2007), doi:10.1109/MCSE.2007.55 |

Machine-readable versions: [`references.bib`](references.bib) and [`../CITATION.cff`](../CITATION.cff).

## 3. Suggested wording for papers or posters

> We reanalysed the genotype and phenotype data of Jivanji et al. (2025; Dryad doi:10.5061/dryad.fxpnvx15q) using PLINK 2
> (Chang et al., 2015) with Firth logistic regression, KING relatedness estimation (Manichaikul et al., 2010) and genomic
> control (Devlin & Roeder, 1999). Code and results: Thaufiq (2026), https://github.com/YOUR-USERNAME/cattle-splotchy-face-gwas.
