# DECISIONS.md

## [D-001] Output Artifact Schema

* **Options considered**: 64-dim, 128-dim, 256-dim embeddings; CSV, NPZ, HDF5 formats; float16, float32, float64 dtypes
* **Choice**: 128-dim, HDF5, float32, L2-normalized
* **Rationale**: Dimensionality is set by the downstream consumer's input contract, not
optimized here. 128 is sufficient to represent KEGG module-level structure (~450 modules)
with room for inter-module signal, while keeping downstream graph construction tractable.
256 would increase edge computation cost with no established benefit at this embedding
granularity.

**HDF5**: the only format that supports heterogeneous arrays (protein embeddings, metabolite
embeddings, metadata, cross-references) in a single file with typed datasets and atomic reads.
CSV cannot represent float32 arrays without precision loss. NPZ lacks native metadata and
identifier storage.

**float32**: standard for learned embeddings. float16 introduces rounding error that compounds
across L2 normalization and downstream cosine similarity operations. float64 doubles storage
with no benefit for neural embedding precision.

**L2 normalization**: required for cosine similarity in downstream graph construction to be
equivalent to dot product. Applied per embedding at write time.

Downstream impact: the output schema is a contract. Any change to dimensionality, format,
or identifier fields breaks downstream data loaders. Do not modify without notifying all
downstream consumers.

---

## [D-002] Identifier Namespaces

* **Options considered**: HGNC symbols, Ensembl gene IDs, Entrez gene IDs, UniProt accessions
(proteins); HMDB IDs, PubChem CIDs, ChEBI IDs, KEGG compound IDs (metabolites)
* **Choice**: UniProt canonical accessions (proteins); KEGG compound IDs with PubChem CID
and ChEBI cross-refs stored (metabolites)
* **Rationale**:

**Proteins — UniProt**: downstream benchmarks use UniProt as the protein identifier space.
Aligning output to UniProt eliminates a join step in dependent pipelines and removes a class
of identifier mismatch bugs. HGNC symbols are non-unique across isoforms and have irregular
versioning. Ensembl and Entrez require additional mapping to reach UniProt; using UniProt
directly removes that dependency.

Swiss-Prot reviewed only, human only, canonical isoform only. These filters enforce a
deterministic one-to-one mapping. TrEMBL entries are excluded because annotation quality
is unreviewed and isoform handling is inconsistent. Multi-isoform genes map to the canonical
accession; isoform-specific accessions are dropped and logged.

Translation from HGNC to UniProt happens at write time in `src/export/`. Internal pipeline
operates on HGNC throughout. Downstream consumers receive UniProt accessions and never
see HGNC symbols. Mapping method: UniProt ID Mapping API → STRING alias reconciliation
for unmapped residuals.

**Metabolites — KEGG compound IDs**: KEGG compound IDs are the native identifier in
KEGG modules, which define the latent space structure. Using any other namespace as primary
would require a translation at the KEGG constraint layer. PubChem CID and ChEBI stored as
cross-references for downstream interoperability.

---

## [D-003] KEGG Version

* **Options considered**: 111.0 (April 2025), 118.0 (May 2026)
* **Choice**: 118.0 (May 2026)
* **Rationale**: 111.0 was the planned version at project scoping but the live KEGG API
  is on 118.0 at time of implementation. Downstream projects have not yet consumed any
  KEGG data. Updating now costs nothing; updating after downstream projects have built
  against a version would require retraining and rebuilding dependent pipelines.
  118.0 is the correct version to lock across the portfolio.

Do not update past 118.0 without: (1) re-deriving module membership, (2) retraining
the VAE, (3) re-exporting the HDF5 artifact, (4) notifying all downstream consumers.

---

## [D-004] Data Source

* **Options considered**: CCLE via DepMap (OmicsExpressionProteinCodingGenesTPMLogp1 +
OmicsMetabolites), TCGA (transcriptomics only, no paired metabolomics), GTEx (no
metabolomics, healthy tissue only), custom institutional datasets (not publicly available)
* **Choice**: CCLE via DepMap — paired transcriptomics and LC-MS metabolomics
* **Rationale**: The only publicly available dataset providing paired transcriptomics and targeted
metabolomics at this sample scale (~900 paired samples). TCGA and GTEx are excluded because
they lack metabolomics entirely. No alternative public dataset reaches comparable paired N
without institutional access.

Known limitations: LC-MS targeted panel covers ~225 metabolites, not untargeted coverage.
Samples are cancer cell lines, not primary tissue — embeddings reflect a cancer-cell-line
biology. This is appropriate for the intended downstream task but limits generalizability
claims to primary tissue or organismal contexts.

DepMap release: 25Q2 for the transcriptomics, CCLE 2019 for metabolomics for the metabolomics. SHA256: 

---

## [D-005] Paired Sample Floor

* **Options considered**: 700, 800, 900
* **Choice**: 800 minimum acceptable paired N
* **Rationale**: Expected paired overlap is ~900 based on current DepMap release metadata.
A floor of 800 allows for ~10% sample loss from QC filtering (expression outliers, missing
metabolite coverage) before the dataset becomes underpowered for VAE training at the planned
latent dimensionality. Below 800, the latent space partition into K KEGG module blocks risks
insufficient samples-per-block for gradient signal. If paired N falls below 800 after QC,
stop and re-evaluate the data source before proceeding.

---

## [D-006] Metabolite Missing Value Imputation

* **Options considered**: zero-fill, KNN imputation, minimum value imputation,
  drop samples with >X% missing
* **Choice**: no action required
* **Rationale**: Source file (CCLE_metabolomics_20190502.csv) is pre-imputed and
  clean at source. Confirmed no NaNs in loaded dataframe. Imputation is not a
  pipeline responsibility.

---

## [D-007] Metabolite Detection Rate Threshold

* **Options considered**: 50%, 70%, 80% detection rate thresholds
* **Choice**: no action required
* **Rationale**: Source file is pre-filtered and clean at source. No metabolites
  need to be dropped on detection rate grounds. Threshold decision is moot.

---

## [OPEN-001] KEGG Constraint Mechanism

* **Options under consideration**: regularization loss term penalizing cross-block covariance;
hard binary masking of decoder weights by KEGG module membership; soft assignment weights
on latent blocks learned during training
* **Working hypothesis**: soft assignment — latent space partitioned into K blocks, each
corresponding to a KEGG module, with overlap handled via learned soft weights
* **Decision gate**: preregistration document before any model code
* **Fallback**: hard masking with pathway priority if soft assignment does not stabilize by week 7

---

## [OPEN-002] PyTorch CUDA Target

* **Options**: CPU-only, CUDA 11.8, CUDA 12.1, CUDA 12.4
* **Decision gate**: confirm hardware before locking torch wheel variant
* **Risk**: installing torch==2.3.1 without specifying CUDA version silently installs the CPU
build on some systems. Verify with `torch.cuda.is_available()` in 00_check_environment.py.

