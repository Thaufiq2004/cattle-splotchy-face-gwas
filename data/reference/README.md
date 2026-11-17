# Reference annotation (not committed)

`scripts/annotate_variants.py` expects:

    data/reference/ARS-UCD1.2.ncbiRefSeq.gtf.gz

* NCBI RefSeq gene annotation for the cattle assembly **ARS-UCD1.2** (RefSeq assembly GCF_002263795.x), GTF format.
* The file used in the original run has `ncbiRefSeq.2025-05-22` in its source column and uses RefSeq accessions as
  sequence names (chromosome 22 = `NC_037349.1`). SHA-256 is in `../CHECKSUMS.sha256`.
* TODO (author): record here the exact download URL and date you used for this file (for example the UCSC Genome Browser
  download area for the ARS-UCD1.2 / bosTau9 RefSeq annotation) so others can fetch the identical release.

Cite: Rosen *et al.* 2020 (ARS-UCD1.2) and O'Leary *et al.* 2016 (RefSeq) - see `docs/ACKNOWLEDGEMENTS_AND_CITATIONS.md`.
