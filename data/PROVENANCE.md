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
* SHA256: bf72dbfaf31ea29bd0072c63999027e791a26a2079b9fb4e03990e208464317e
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

## File 6 — metabolite_kegg_mapping.csv

* Generated: 2026-05-30
* Method: KEGG /find/compound exact synonym match (kegg_exact)
        + manual review and unmapped-set sweep (kegg_manual)
        + PubChem /rest/pug/compound/name fallback (pubchem_only)
* SHA256: 94c1b96a872114a85546c9fe859647f1cf908e7e2fd7bf2b38218b231c8ba2d6
* Coverage:
    kegg_exact:   81/225  (36.0%)
    kegg_manual:  23/225  (10.2%)
    pubchem_only: 16/225  (7.1%)
    unmapped:    105/225  (46.7%)
* KEGG-mapped total: 104/225 (46.2%)
* Cross-ref total:   120/225 (53.3%)
* Unmapped breakdown:
    lipid species — KEGG class-level only:    89 (LC-MS C##:# LPC/LPE/PC/SM/DAG/CE/TAG)
    isobaric mixture — no single KEGG ID:     13 (e.g., DHAP/glyceraldehyde 3P)
    no KEGG entry:                             3
* Notes:
    - Unmapped metabolites are retained in training per D-013; they appear
      as zero-rows in the metabolite-module membership matrix (Day 4).
    - Coverage limit is structural — LC-MS lipid species below KEGG's
      compound-level resolution. See D-012, OPEN-007.