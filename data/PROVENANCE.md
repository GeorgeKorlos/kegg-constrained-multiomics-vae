# PROVENANCE.md

## DepMap Release

* **Release**: DepMap Public 25Q2
* **Download date**: 2026-05-26
* **Source**: https://depmap.org/portal/download/all/
* **Citation**: DepMap, Broad (2025). DepMap Public 25Q2. Dataset. depmap.org

## File 1 — OmicsExpressionProteinCodingGenesTPMLogp1.csv

* **SHA256**: e0326e16eb23bea1be980fce315acb36b224dedd7af6b47e0ba37e7747dbcc47
* **Shape**: 1684 × 19205
* **Notes**: Row index = DepMap IDs (ACH-XXXXXX). Column headers = HGNC symbols
  with Entrez IDs appended (e.g. TSPAN6 (7105)) — strip Entrez suffix before
  HGNC → UniProt mapping in export step.

## File 2 — CCLE_metabolomics_20190502.csv

* **Release**: CCLE 2019
* **Citation**: Li et al. The landscape of cancer cell line metabolism.
  Nature Medicine 25, 850-860 (2019).
* **SHA256**: 7c1d24aa575f4c58a29019026b5df8e6d1142a56925aba32ff3f1d1d5a7fd0ac
* **Shape**: 928 × 226 (225 metabolites + CCLE_ID column)
* **Notes**: Pre-imputed at source (clean, no NaNs expected). Values are
  log-transformed. Column headers are metabolite names, not HMDB IDs or
  KEGG compound IDs — name → KEGG compound ID translation required in week 2.

## Paired Dataset

* **Join key**: DepMap ID
* **Join type**: inner
* **Paired N**: 912
* **Status**: PASS (above 800 floor)

## File 3 — kegg_module_gene_membership.tsv
* Release: KEGG 118.0+/05-29
* Fetch date: 2026-05-29
* Source endpoints: /link/hsa/module + /list/hsa (HGNC symbol map)
* SHA256: 77741ce0fdbf86bc2dedfeb2632633fdf196bab112def386ba4f094578d68e82
* Rows: 1075
* Unique modules: 131
* Unique HGNC symbols: 832
* Symbol-map misses: 1426 (pseudogenes, uncharacterized loci)
* Join misses: 3 (deprecated KEGG entries)

## File 4 — kegg_module_compound_membership.tsv
* SHA256: a85688084f904667ae3667a099a5648267e20eb295f68e1737383b4fffbdafe3
* Source endpoint: /link/compound/module
* Rows: 3537
* Unique modules: 507
* Unique compounds: 2031

## File 5 — kegg_module_metadata.tsv
* SHA256: ee15b2c7a0bcca391f30ef87329c27fe6945c6a0591cb800dfd78a6a29909aaf
* Source endpoint: /list/module
* Total modules: 573