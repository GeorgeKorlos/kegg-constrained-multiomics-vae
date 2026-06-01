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

## File 7 — gene_module_matrix.npy
* Generated: 2026-05-30
* Shape: (17384, 235)
* dtype: float32
* SHA256: 1463998ae812e6c139c6b4062a13e9cefe936a858096bb0bbab82aafe239f2af
* Notes: row order matches data/processed/gene_names.csv;
        column order matches data/processed/module_ids.csv;
        1006 gene-module links across 130 CCLE-detected modules.

## File 8 — metabolite_module_matrix.npy
* Generated: 2026-05-30
* Shape: (225, 235)
* dtype: float32
* SHA256: 91ce140473933a4443135255dff55a83bd59f9ab9c0069602a1aa4cb95d4afa4
* Notes: row order matches data/processed/metabolite_names.csv;
        column order matches data/processed/module_ids.csv;
        250 metabolite-module links across 154 CCLE-detected modules;
        171 unmapped+module-orphan metabolites have all-zero rows (D-013).

## File 9 — module_ids.csv
* Generated: 2026-05-30
* Rows: 235
* SHA256: 17c97f54083fb7fe1be97ab306f1ce52e5cbd1afdbe3f25ec9f4f424477143c2
* Notes: lexically sorted union of CCLE-detected KEGG modules
        (130 gene-side + 154 metabolite-side, 49 in both).
        Column order for both membership matrices.

## File 10 — pathway_metabolite_membership.csv

* Generated: 2026-05-31
* Source: `scripts/05b_pathway_validation_set.py`
* Method: KEGG /link/compound/pathway (reference map IDs, species-agnostic)
        filtered to four target pathway maps and CCLE kegg_exact + kegg_manual
* SHA256: 461bb1bc4e946e9a14254d4b55910e8bd41bc428f42566959acda813bc37e3a7
* Rows: 28
* Per-pathway counts:
    map00010 (Glycolysis): 2 (descriptive only — CCLE panel limit)
    map00020 (TCA cycle): 6
    map00260 (Glycine/Serine/Threonine): 10
    map00670 (One-carbon pool by folate): 10
* Per D-015 — pathway-map granularity for Section 4 Step 3 validation.    
  
## File 11 — block_sizes.csv

* Generated: 2026-05-31
* Source: `scripts/05a_compute_block_size.py`
* Shape: 49 rows × 5 columns
* Columns: module_id, module_name, s_k (gene count), s_k_metabolites, block_size
* SHA256: e51f28cad25be12217dd6fb81613b23bc677685a90eb0cf72ab455c7462e5f4f
* Per D-011 (revised) and D-014. Latent block allocation table:
    K = 49 modules
    Σ b_k = 128
    min(b_k) = 2
    max(b_k) = 4
    median(b_k) = 3
    Modules at floor (b_k = 2): 25
* Notes: gene-count-weighted allocation. s_k_metabolites is informational
        (allocation is gene-weighted because gene branch has substantially
        more KEGG-annotated entities: 778 vs 54).

## Processed Training Tensors (derived — gitignored, reproducible)

* Files: transcriptomics.npy, metabolomics.npy, transcriptomics_scaler.pkl,
         metabolomics_scaler.pkl
* Reproduction: `python src/data/transform.py`
* Inputs determining output bytes (all pinned):
    - data/raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv (SHA e0326e16…)
    - data/raw/CCLE_metabolomics_20190502.csv                (SHA 7c1d24aa…)
    - data/splits/split_v1.json (SHA 85292785…) — scaler fits on TRAIN rows only (D-016)
    - src/data/transform.py @ commit 1d9d49e630c9902d97919179c6418f27d307889a
* Scaler: fit on 720 train rows only. Supersedes prior fit-on-full-898 artifacts (D-016).
* Row order: inner-join on DepMap_ID; stability asserted in build_paired_dataset
  against data/processed/sample_ids.csv.