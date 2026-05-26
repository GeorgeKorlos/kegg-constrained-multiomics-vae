# protocol/data_source.md

---

## Source 1 — Transcriptomics

### Dataset Identity

* **Dataset name**: CCLE Omics Expression — Protein Coding Genes TPM Log1p
* **Source**: https://depmap.org/portal/download/all/

### File

* **Filename**: OmicsExpressionProteinCodingGenesTPMLogp1.csv
* **Format**: CSV
* **Role in pipeline**: transcriptomics modality input to VAE encoder. Gene-level expression
  values pre-transformed as log1p(TPM) at source. Rows = cell lines (DepMap IDs),
  columns = HGNC gene symbols. No additional log transformation applied in pipeline.
  HGNC symbols translated to canonical UniProt accessions at export time in `src/export/`.

### Dataset Characteristics

* **Type**: bulk RNA-seq, gene-level summarized
* **Scale**: ~1,019 cell lines (pre-pairing)
* **Taxonomic scope**: Homo sapiens (cancer cell lines, CCLE)
* **Identifier space**: HGNC gene symbols (column headers)
* **Value space**: log1p(TPM) — pre-transformed at source, do not re-transform
* **Preprocessing constraints**:
  - Retain protein-coding genes only (enforced by source file selection)
  - HGNC → UniProt translation applied at write time in `src/export/`, not during
    preprocessing

---

## Source 2 — Metabolomics

### Dataset Identity

* **Dataset name**: CCLE Omics Metabolites CCLE 2019
* **Source**: https://depmap.org/portal/download/all/

### File

* **Filename**: CCLE_metabolomics_20190502.csv
* **Format**: CSV
* **Role in pipeline**: metabolomics modality input to VAE encoder. LC-MS targeted panel.
  Rows = cell lines (DepMap IDs), columns = metabolite identifiers. Missing values present
  — imputation strategy and detection rate threshold decided after day 3 QC (OPEN-003,
  OPEN-004). Metabolite identifiers mapped to KEGG compound IDs; PubChem CID and
  ChEBI stored as cross-references in HDF5 output.

### Dataset Characteristics

* **Type**: LC-MS targeted metabolomics
* **Scale**: ~900 cell lines (pre-pairing); ~225 metabolites
* **Taxonomic scope**: Homo sapiens (cancer cell lines, CCLE)
* **Preprocessing constraints**:
  - Detection rate threshold: tentatively 50%, confirmed after day 3 QC (OPEN-004)
  - Missing value imputation: strategy TBD after day 3 QC (OPEN-003)
  - Metabolite identifiers mapped to KEGG compound IDs for KEGG constraint layer

---

## Paired Dataset

* **Join key**: DepMap ID
* **Join type**: inner — only cell lines present in both modalities retained
* **Expected paired N**: ~900
* **Minimum acceptable paired N**: 800 (D-005)
* **Action if below 800**: stop, do not proceed to QC, re-evaluate data source

---

## Source 3 — KEGG Module Membership

### Dataset Identity

* **Dataset name**: KEGG Module Database
* **Release**: 111.0 (April 2025) — locked, see D-003
* **Source**: https://rest.kegg.jp/

### Files (serialized from API, saved to `data/raw/`)

* **Filename**: kegg_module_gene_membership.tsv
* **Format**: TSV — columns: [module_id, hgnc_symbol]
* **Fetch method**: KEGG REST API — `/link/hsa/module`, `/list/module`
* **Role in pipeline**: defines gene-to-module membership for the latent block partition.

* **Filename**: kegg_module_compound_membership.tsv
* **Format**: TSV — columns: [module_id, kegg_compound_id]
* **Fetch method**: KEGG REST API — `/link/compound/module`
* **Role in pipeline**: defines compound-to-module membership for the metabolomics
  branch of the KEGG constraint.

### Dataset Characteristics

* **Type**: pathway/module gene and compound membership tables
* **Scale**: ~450 KEGG modules (human); subset with detected genes/metabolites used
* **Identifier space**: KEGG module IDs; HGNC symbols (gene branch); KEGG compound
  IDs (metabolite branch)
* **Preprocessing constraints**:
  - Restrict to human modules (hsa prefix)
  - Drop modules with zero detected genes in transcriptomics data
  - Drop modules with zero detected metabolites in metabolomics data
  - Single-modality modules retained but flagged in QC report

---

## Source 4 — UniProt ID Mapping

### Dataset Identity

* **Dataset name**: UniProt ID Mapping — HGNC to UniProt canonical accessions
* **Source**: https://rest.uniprot.org/idmapping/

### File (saved output of mapping job, `data/raw/`)

* **Filename**: uniprot_hgnc_mapping.tsv
* **Format**: TSV — columns: [hgnc_symbol, uniprot_accession, status]
* **Fetch method**: UniProt ID Mapping API — programmatic submission, results saved
  immediately on retrieval
* **Role in pipeline**: primary HGNC → UniProt translation table used in `src/export/`.
  Filters: Swiss-Prot reviewed only, Homo sapiens only, canonical isoform only.
  Unmapped symbols passed to STRING alias reconciliation (Source 5). Remaining
  unmapped symbols dropped and logged to `/protein_metadata/dropped_ids` in HDF5.

### Dataset Characteristics

* **Type**: identifier mapping table
* **Scale**: ~18,000 HGNC input symbols → ~16,500 UniProt canonical accessions
  expected (~91.7% coverage)
* **Filters applied**:
  - Swiss-Prot reviewed only (TrEMBL excluded)
  - Homo sapiens only (taxon 9606)
  - Canonical isoform only

---

## Source 5 — STRING Alias File

### Dataset Identity

* **Dataset name**: STRING Protein Aliases
* **Source**: https://string-db.org/cgi/download?sessionId=&species_text=Homo+sapiens

### File

* **Filename**: 9606.protein.aliases.v[version].txt.gz
* **Format**: TSV (compressed) — columns: [string_id, alias, source]
* **Role in pipeline**: fallback reconciliation for HGNC symbols not resolved by the
  UniProt ID Mapping API. Resolves updated or deprecated HGNC symbols. Remaining
  unmapped symbols after STRING reconciliation are dropped and logged.

### Dataset Characteristics

* **Type**: protein identifier alias table
* **Taxonomic scope**: Homo sapiens (taxon 9606)
* **Identifier space**: STRING protein IDs → alias symbols (HGNC, Ensembl, Entrez,
  synonyms)
* **Preprocessing constraints**:
  - Filter to alias sources relevant to HGNC reconciliation
  - Apply same Swiss-Prot reviewed / canonical isoform filters as Source 4

---

## Known Limitations

* LC-MS targeted panel covers ~225 metabolites. Modules without any detected
  metabolite contribute constraint signal through the transcriptomics branch only.
* Samples are cancer cell lines, not primary tissue. Appropriate for P5's DTI task
  (DAVIS/KIBA) but not generalizable to primary tissue contexts without retraining.
* No source-provided checksums from DepMap portal. Local SHA256 is the sole
  integrity record for Sources 1 and 2.
* KEGG REST API has no versioned download endpoint. Membership tables are pinned
  by SHA256 after serialization. Re-fetching against a different release is a breaking
  change for P3 and P5.