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
* **Confirmed post-QC paired N**: 898 (above 800 floor, PASS)

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

## [D-008] Low-Variance Gene Threshold

* **Options considered**: 0.01 (987 genes, 5.1%), 0.05 (1821 genes, 9.5%),
  0.10 (2328 genes, 12.1%), 0.20 (3004 genes, 15.6%), 0.50 (7625 genes, 39.7%)
* **Choice**: variance < 0.05 — drop 1821 genes, retain 17,384
* **Rationale**: The variance distribution shows a dense cluster of near-zero
  variance genes below 0.05. Above this threshold the removal becomes a judgment
  call rather than a clear uninformative boundary. 0.01 is too conservative —
  genes with variance 0.01–0.05 have essentially flat expression across 1684
  cell lines and contribute no gradient signal to the VAE encoder. 0.50 removes
  39.7% of the feature space, which is too aggressive. 0.05 captures the natural
  low-variance cluster without excessive feature loss.

---

## [D-009] Expression Outlier Threshold — Transcriptomics

* **Options considered**: 2 SD, 3 SD, IQR-based
* **Choice**: 3 SD from mean total expression per sample
* **Rationale**: 15 samples flagged out of 1684 (0.9%). Distribution of sample
  total expression is tight (mean 51346, std 3731) — flagged samples are genuine
  outliers, not artifacts of a skewed distribution. 3 SD is a standard conservative
  threshold that avoids over-removing samples while catching true anomalies.

---

## [D-010] Expression Outlier Threshold — Metabolomics

* **Options considered**: 2 SD, 3 SD, IQR-based
* **Choice**: 3 SD from mean total metabolite level per sample
* **Rationale**: 11 samples flagged out of 928 (1.2%). Metabolomics sample totals
  are extremely tight (mean 1322, std 7.4) — the data is heavily normalized at
  source. Flagged samples are clear anomalies. No overlap with transcriptomics
  outliers — independent quality issues in each modality.
  
---

## [D-011] Latent block size allocation

* Options considered:
  - Equal-size blocks (b = 128/K) — fails when K > 64
  - Proportional allocation with min floor — chosen
  - Increase latent_dim to 512 — rejected, breaks P3 contract
* Choice: b_k = max(2, round(128 * s_k / Σ s_k)), adjusted to sum to 128
* Rationale: K is at most 131 pre-CCLE-filter, likely 100-120 post-filter.
  Equal blocks force b < 2 which is not a meaningful sub-vector. Proportional
  allocation aligns latent capacity with module size, preserves the 128-dim
  contract, and requires no change to the loss formulation — Frobenius norm
  of cross-block covariance is shape-agnostic.
* Block size table to be generated and logged at Day 4 alongside K.
* Revised retention: only "both" modules (49) retained for block partition.
  Single-modality modules (gene-only 81, metabolite-only 105) excluded
  from the latent block partition to keep b_k ≥ 2. They are documented
  in the membership matrices but are not assigned to latent dimensions.
* **Confirmed at Day 4:** K = 235 (49 modules in both modalities, 81 gene-only,
  105 metabolite-only). This is substantially larger than the 100–120 estimate
  at preregistration time.
* **Block allocation problem:** with K = 235, minimum total dim if every block
  gets b_k ≥ 2 is 470 > 128. The min-floor of 2 is incompatible with K = 235
  at latent_dim = 128.
* **Decision deferred to Day 5:** the choice between (a) dropping the min-floor
  (allowing b_k = 1 for small modules), or (b) retaining only the 49 "both"
  modules for the latent block partition while keeping single-modality modules
  in the membership matrices for reporting, is logged for Day 5 resolution
  alongside the coverage report. The cleaner architectural path is (b) — see
  Day 5 deliverable.
* **Resolved at Day 5:** K = 49 (both-modality only). Block sizes computed
  from gene counts per module, min floor b_k = 2, Σ b_k = 128. min(b_k) = 2,
  max(b_k) = 4, median = 3. 25 of 49 modules at the floor.
  Per-module table in `data/processed/block_sizes.csv` (SHA e51f28cad25be12217dd6fb81613b23bc677685a90eb0cf72ab455c7462e5f4f).
  See `reports/kegg_coverage.md` for full statistics.

---


## [D-012] Metabolite KEGG coverage limit

* **Observed**: 104/225 (46.2%) of CCLE metabolites map to a KEGG compound ID
  (81 via /find/compound exact synonym match, 23 via manual review and
  unmapped-set sweep). 16 additional metabolites have PubChem cross-references
  only. 105 are fully unmapped.
* **Cause**: ~85% of unmapped entries are LC-MS-resolved lipid species
  (C##:# LPC/LPE/PC/SM/DAG/CE/TAG). KEGG catalogs lipid classes, not
  species-level resolution. ~12% are isobaric mixtures (e.g.,
  DHAP/glyceraldehyde 3P, UDP-galactose/UDP-glucose) the instrument cannot
  separate into single compounds.
* **Consequence**: ~54% of CCLE metabolites are unconstrained by L_KEGG —
  they contribute to reconstruction loss but receive no signal from the
  KEGG regularization term.
* **Mitigation**: none feasible. Lipid species are not in KEGG at this
  resolution. Isobaric mixtures cannot be assigned to a single compound.
* **Implication for preregistration**: Section 4 Step 3 (pathway activity
  validation) coverage is bounded by the 104 mapped metabolites. All four
  named pathways (glycolysis, serine biosynthesis, TCA, one-carbon) fall
  within the mapped subset. Validation set is sufficient. Section 6
  updated to flag coverage constraint explicitly.
* Additionally observed: of the 104 KEGG-mapped metabolites, only 54 
  belong to any KEGG module. The other 50 have valid KEGG compound IDs 
  but are not catalogued in any module (e.g., NADP, 1-methylnicotinamide, 
  2-hydroxyglutarate). Effective KEGG-constrained subset is 54/225 
  (24.0%). The 50 module-orphan metabolites are treated identically to 
  the 105 fully-unmapped metabolites per D-013: retained in training, 
  contribute to reconstruction loss, receive zero-rows in the 
  metabolite-module matrix.
* **Refined at Day 4:** of the 104 KEGG-mapped metabolites, only 54 are
  members of any KEGG module. The other 50 have valid KEGG compound IDs
  but no module annotation (e.g., NADP, 2-hydroxyglutarate, 4-pyridoxate).
  Effective KEGG-constrained metabolite subset is 54/225 (24.0%), not 46.2%.
  The 50 module-orphan metabolites are treated identically to the 105
  fully-unmapped metabolites per D-013.

---

## [D-013] Unmapped metabolites retained in training, unconstrained by KEGG

* **Context**: 105 of 225 CCLE metabolites have no KEGG compound ID
  (see D-012).
* **Options considered**:
  - (a) Retain in training — pass through encoder/decoder, contribute to
        reconstruction loss, receive zero-rows in metabolite-module
        membership matrix so L_KEGG does not act on them.
  - (b) Drop from input — train on 104 mapped metabolites only.
* **Choice**: (a) — retain all 225, mark unmapped as zero-rows in
  membership matrix.
* **Rationale**:
  - Output contract is locked at 225 metabolite embeddings (P3-facing).
    Dropping unmapped would break the contract and require portfolio
    renegotiation.
  - The architecture handles unconstrained features by construction
    (zero rows in soft assignment matrix produce no KEGG gradient on
    those entities).
  - Reconstruction signal still produces usable embeddings for downstream.
  - Drift risk on the unconstrained subset is logged in OPEN-007.
* **Implication for evaluation**: Section 4 Step 2 (pathway coherence)
  runs over the mapped subset only. Unmapped metabolites are not part of
  the in-module vs out-of-module cosine comparison.

---

---

## [D-014] Single-modality modules excluded from latent partition

* **Context:** Day 4 found 235 KEGG modules with at least one CCLE entity.
  Of these, only 49 have both gene and metabolite annotations. 81 are
  gene-only, 105 are metabolite-only.
* **Options considered:**
  - (a) K = 235 (full union) — accept b_k = 1 for many blocks.
  - (b) K = 49 (both-modality only) — single-modality modules documented
        in the membership matrices but excluded from latent partition.
* **Choice:** (b)
* **Rationale:**
  - The KEGG constraint's theoretical justification is cross-modality
    coherence. Single-modality blocks cannot be evaluated for cross-block
    KEGG-aligned variance.
  - The loss term L_KEGG = Σ_{k≠l} ||Cov(z_k, z_l)||²_F is degenerate
    when blocks have no anchor in one modality.
  - K = 49 with b_k ≥ 2 fits cleanly within latent_dim = 128.
  - Section 4 Step 2 (cross-modality coherence test) requires both-modality
    annotation to be meaningful.
* **Implications:**
  - Section 4 Step 2 evaluates over 49 modules.
  - 186 single-modality modules are retained in `gene_module_matrix.npy`
    and `metabolite_module_matrix.npy` (235 columns) for downstream use
    (P3 may consume the full membership table independently).
  - Latent partition obtained at model load time by slicing matrices to
    49 columns where both matrices have non-zero column sums.

---

## [D-015] Pathway activity score validation at pathway-map granularity

* **Context:** KEGG module-level coverage of the CCLE-mapped metabolite set
  is thin for several validation pathways (Section 4 Step 3): glycolysis
  module M00001 has 2 mapped metabolites, serine biosynthesis module M00020
  has 2, TCA module M00009 has 4, one-carbon module M00141 not in retained 49.
  This is a structural limit of CCLE LC-MS coverage vs KEGG module definitions
  (see D-012), not a matrix construction failure.
* **Options considered:**
  - (a) Pathway-map granularity (`map00010`, `map00020`, `map00260`, `map00670`).
  - (b) Union of related modules per pathway.
  - (c) Accept module-level thin coverage; demote Section 4 Step 3 to
        descriptive.
* **Choice:** (a)
* **Rationale:**
  - Pathway maps are KEGG-native broader units that include peripheral
    metabolites absent from tight module definitions.
  - Per-pathway coverage rises from 2-4 to 6-10 metabolites (3 of 4 pathways),
    enabling meaningful Pearson and Spearman correlations.
  - Decouples the structural prior (module-level, K=49, D-014) from the
    validation metric (pathway-map-level). Fine prior + broad validation
    is a stronger methodological setup than same-granularity for both.
  - Does not modify the latent partition or L_KEGG. Affects only
    Section 4 Step 3 computation.
* **Implementation:** `scripts/05b_pathway_validation_set.py` fetches
  `/link/compound/pathway` (reference map IDs, species-agnostic), filters
  to four named pathway maps and the 104 KEGG-mapped CCLE metabolites,
  saves `data/processed/pathway_metabolite_membership.csv`.
* **Per-pathway counts:**
  - map00010 (Glycolysis): 2 (descriptive only — CCLE panel limit)
  - map00020 (TCA cycle): 6
  - map00260 (Glycine/Serine/Threonine): 10
  - map00670 (One-carbon pool by folate): 10
* **Glycolysis caveat:** map00010 retains 2 mapped metabolites
  (3-phosphoglycerate, PEP) because:
  - Phosphorylated hexose intermediates appear as isobaric mixtures in
    CCLE (D-012).
  - Pyruvate is not in the CCLE targeted LC-MS panel.
  Glycolysis pathway activity is reported as **descriptive, not inferential**.
  Section 4 Step 3 pass condition (Pearson r > 0 for ≥3 of 4 pathways) is
  unaffected.
* **Note on prefix:** Pathway-compound links from `/link/compound/pathway`
  return `map`-prefixed IDs (species-agnostic). The four IDs above use the
  `map` prefix accordingly. Human specificity enters via the CCLE-mapped
  metabolite filter on the join, not via the pathway ID prefix.

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

## [OPEN-005] Encoder Architecture Symmetry

* **Flag**: Transcriptomics and metabolomics have different intrinsic
  dimensionalities — transcriptomics requires >50 PCs for 80% variance,
  metabolomics reaches 80% within ~48 PCs. A symmetric encoder design
  may underserve the transcriptomics modality.
* **Options under consideration**: symmetric encoder (same depth and width
  for both modalities), asymmetric encoder (deeper or wider transcriptomics
  encoder), shared encoder with modality-specific input projection layers
* **Decision gate**: preregistration document before any model code

---

## [OPEN-007] Metabolomics input scaling vs. lipid variance dominance

* **Context**: ~54% of metabolites are unconstrained by L_KEGG and
  retained via reconstruction loss only (D-013). LC-MS lipid species
  typically have higher cross-sample variance than central metabolites,
  so reconstruction loss — which is variance-weighted by construction —
  may be dominated by lipid features. The KEGG term operates on the
  latent partition, but the partition is shaped by what the encoders
  push into it, which is shaped by reconstruction pressure. Risk: KEGG
  signal becomes secondary structure in a lipid-dominated latent manifold.
* **Current state**: metabolomics is per-feature standardized
  (mean 0, std 1) in `src/data/transform.py`. This partially mitigates
  the variance asymmetry at the input level but does not guarantee
  balanced gradient at the loss level.
* **Decision gate**: after the first ablation run completes, inspect
  the lipid-driven drift diagnostic (preregistration Section 7). If
  KEGG-block variance is suppressed below the unconstrained-residual
  variance, consider:
    (a) re-weight reconstruction loss to upweight KEGG-mapped metabolites
    (b) variance-rank truncation on metabolomics input — REJECTED in
        advance, breaks output contract
    (c) accept the result and report lipid-dominant latent as a finding
* **Architecture fallback**: hard masking (already in ablation matrix
  as condition C) is the canonical alternative if soft assignment
  proves dominated by unconstrained variance.
* **Do not act before the data tells you to.**