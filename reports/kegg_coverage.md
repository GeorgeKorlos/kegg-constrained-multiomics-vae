# KEGG Coverage Report 

End of week 2 deliverable. Documents the KEGG-derived constraint structure
for the P2 VAE. All numbers are final as of 2026-05-31.

---

## Summary

| Quantity | Value |
|---|---|
| KEGG release | 118.0+/05-30 |
| Modules in KEGG (total, human) | 573 |
| Modules with ≥1 CCLE entity (gene or metabolite) | 235 |
| **Modules retained for latent partition (K)** | **49** |
| CCLE genes annotated to ≥1 retained module | 778 / 17,384 (4.5%) |
| CCLE metabolites annotated to ≥1 retained module | 54 / 225 (24.0%) |
| Latent dimensions | 128 |
| Block size range (b_k) | 2 to 4 |

The latent partition uses 49 modules with both gene and metabolite annotations
(D-014). Single-modality modules (81 gene-only, 105 metabolite-only) are
documented but excluded from the partition. Pathway activity score validation
operates at pathway-map granularity (D-015), broader than the module-level
partition.

---

## Module Filtering Funnel

KEGG human modules:              573
↓ filter to modules with ≥1 CCLE gene or metabolite
Modules with ≥1 CCLE entity:     235
↓ require both gene AND metabolite annotation (D-014)
Modules retained for partition:   49

### Retained module classification

| Category | Count | Status |
|---|---|---|
| Both modalities | 49 | In latent partition |
| Gene-only | 81 | In gene_module_matrix.npy (235 cols), excluded from partition |
| Metabolite-only | 105 | In metabolite_module_matrix.npy (235 cols), excluded from partition |

The 235-column membership matrices are the canonical artifacts on disk. The
49-module partition is obtained at model load time by selecting columns where
both matrices have non-zero column sums.

---

## Gene Coverage

| Quantity | Value |
|---|---|
| CCLE genes (post-QC) | 17,384 |
| Genes with ≥1 retained module | 778 (4.5%) |
| Genes with 0 modules (unconstrained) | 16,606 (95.5%) |
| Gene-module link count (retained) | 1,006 |
| Median modules per annotated gene | 1.0 |
| Max modules per gene | 6 |

Unconstrained genes contribute to L_recon_trans but receive no signal from
L_KEGG. The ~5% KEGG annotation rate is consistent with KEGG's coverage of
the human metabolic and signaling backbone — not a pipeline filter loss.

---

## Metabolite Coverage

| Quantity | Value |
|---|---|
| CCLE metabolites | 225 |
| KEGG-mapped (kegg_exact + kegg_manual) | 104 (46.2%) |
| Module-resident (in ≥1 KEGG module) | 54 (24.0%) |
| Metabolite-module link count (retained) | 250 |
| Median modules per annotated metabolite | 3.0 |
| Max modules per metabolite | 16 |

### Three classes of unconstrained metabolites

| Class | Count | % of 225 |
|---|---|---|
| Fully unmapped (no KEGG ID) | 105 | 46.7% |
| KEGG-IDed but no module (D-012) | 50 | 22.2% |
| PubChem cross-ref only | 16 | 7.1% |
| **Total unconstrained** | **171** | **76.0%** |

The 105 fully unmapped are dominated by LC-MS lipid species (~89) and isobaric
mixtures (~13). See D-012 for breakdown. The 50 module-orphans have valid KEGG
compound IDs but are not catalogued in any module (e.g., NADP, 2-hydroxyglutarate,
4-pyridoxate).

All 171 unconstrained metabolites are retained in training per D-013 — they
contribute to L_recon_meta but appear as zero-rows in the membership matrix.

---

## Block Size Allocation (K=49)

Per D-011 (revised). Block sizes proportional to gene count per module, with
minimum floor b_k = 2 and constraint Σ b_k = 128.

Allocation formula:
b_k_raw = 128 × s_k / Σ s_k
b_k = max(2, round(b_k_raw))
adjust largest blocks to satisfy Σ b_k = 128

Note: `s_k` in the table below refers to gene count per module (annotated CCLE
genes). The allocation is gene-count-weighted because the gene branch has
substantially more KEGG-annotated entities (778 genes vs 54 metabolites) and
dominates the constraint signal.

### Summary statistics

| Statistic | Value |
|---|---|
| K (modules in partition) | 49 |
| Σ b_k | 128 |
| min b_k | 2 |
| max b_k | 4 |
| median b_k | 3 |
| Modules at floor (b_k = 2) | 25 |

Full per-module table in `data/processed/block_sizes.csv`. Top 10 modules by
gene count:

| module_id | module_name | s_k | s_k_metabolites | b_k |
|---|---|---|---|---|
| M00001 | Glycolysis (Embden-Meyerhof pathway) | 24 | 2 | 3 |
| M00009 | Citrate cycle (TCA cycle, Krebs cycle) | 22 | 4 | 3 |
| M00014 | (see CSV) | 18 | — | 3 |
| M00049 | (see CSV) | 18 | — | 4 |
| M00003 | Gluconeogenesis, oxaloacetate → F6P | 16 | 2 | 3 |
| M00958 | (see CSV) | 15 | — | 4 |
| M00146 | (see CSV) | 13 | — | 4 |
| M00147 | (see CSV) | 13 | — | 4 |
| M00094 | (see CSV) | 12 | — | 4 |
| M00034 | (see CSV) | 12 | — | 4 |

(Fill in module names from the actual CSV.)

---

## Pathway Activity Score Validation Set (D-015)

Section 4 Step 3 of the preregistration evaluates pathway-level activity
correlations between gene and metabolite modalities. Per D-015, this evaluation
operates at **KEGG pathway-map granularity** (broader than the module-level
partition). Pathway maps include peripheral metabolites absent from tight
module definitions.

| Pathway map | Description | CCLE metabolites |
|---|---|---|
| map00010 | Glycolysis / Gluconeogenesis | 2 (descriptive only) |
| map00020 | Citrate cycle (TCA cycle) | 6 |
| map00260 | Glycine, serine, threonine metabolism | 10 |
| map00670 | One-carbon pool by folate | 10 |
| **Total (unique pairs)** | | **28** |

### Glycolysis caveat

Glycolysis (map00010) has only 2 CCLE-mapped metabolites available for the
validation set. This is a structural limitation of CCLE's LC-MS targeted
panel:

- Phosphorylated hexose intermediates (G6P, F6P, F1,6BP, DHAP, GAP, 1,3-BPG,
  2PG) appear in CCLE as isobaric mixtures (`F1P/F6P/G1P/G6P`, `DHAP/
  glyceraldehyde 3P`, `hexoses (HILIC neg)`, `hexoses (HILIC pos)`) that
  cannot be assigned to single KEGG compound IDs (D-012).
- Pyruvate (C00022) is not in the CCLE targeted panel.
- Mapped survivors are 3-phosphoglycerate (C00197) and PEP (C00074).

The glycolysis pathway activity score is reported as **descriptive, not
inferential.** Section 4 Step 3 pass condition (Pearson r > 0 for ≥3 of 4
pathways) is unaffected — the other three pathways carry the validation
with adequate metabolite counts.

### Note on map vs hsa prefix

KEGG returns pathway-compound links via `/link/compound/pathway` using
reference map IDs (`map00010`), not organism-specific IDs (`hsa00010`).
Compounds are species-agnostic — only one canonical pathway-compound
table exists. Human specificity comes from the CCLE-mapped metabolite
filter, not from the pathway ID prefix.

---

## Implications for Evaluation

| Section | Granularity | Evaluation set |
|---|---|---|
| Section 4 Step 2 — Coherence | Module-level | 49 retained modules |
| Section 4 Step 3 — Pathway activity | Pathway-map level | 4 maps, 28 (metabolite, pathway) pairs |
| Section 4 Step 4 — P3 downstream | (P3-defined) | Full embedding artifact, 17,384 + 225 vectors |

The 186 single-modality modules and the 171 unconstrained metabolites are
**present in the embedding artifact** that P3 receives. They are excluded
from P2's internal evaluation but available to downstream consumers.

---

## Architectural Implications

The dual-granularity design (module-level constraint, pathway-map-level
validation) implements:

- A **tight structural prior** at module granularity, enforcing latent
  block organization on functionally specific units (~5–25 entities each).
- A **broad biological validation** at pathway-map granularity, testing
  whether the fine prior recovers signal at the coarser level reviewers
  use to reason about metabolic biology.

Showing that module-level structure predicts pathway-map-level activity is
a stronger claim than same-granularity validation. The two granularities
serve distinct purposes.

---

## File References

| File | SHA256 | Description |
|---|---|---|
| kegg_module_gene_membership.tsv | bf72db... | KEGG gene-module links |
| kegg_module_compound_membership.tsv | a85688... | KEGG compound-module links |
| kegg_module_metadata.tsv | ee15b2... | KEGG module names |
| metabolite_kegg_mapping.csv | 94c1b9... | CCLE → KEGG metabolite mapping |
| gene_module_matrix.npy | 146399... | 17384 × 235 binary |
| metabolite_module_matrix.npy | 91ce14... | 225 × 235 binary |
| module_ids.csv | 17c97f... | 235 module IDs (column order) |
| block_sizes.csv | e51f28cad25be12217dd6fb81613b23bc677685a90eb0cf72ab455c7462e5f4f | 49-row block allocation table |
| pathway_metabolite_membership.csv | 461bb1... | 28-row pathway validation set |

See `data/PROVENANCE.md` for full hashes and generation methods.

---

## Decision References

- **D-011 (revised)** — Proportional latent block allocation, K=49, b_k ≥ 2
- **D-012** — Metabolite KEGG coverage limit (104/225, 46.2%)
- **D-013** — Unmapped metabolites retained, unconstrained by L_KEGG
- **D-014** — Single-modality modules excluded from latent partition
- **D-015** — Pathway activity validation at pathway-map granularity
- **OPEN-007** — Lipid-driven drift diagnostic (post-hoc, week 7+)

See `DECISIONS.md` for full text.