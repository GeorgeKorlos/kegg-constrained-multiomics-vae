## Transcriptomics

* **Input shape**: 1684 × 19205
* **Missing values**: 0
* **Gene variance threshold**: 0.05 (D-008)
* **Genes dropped (low variance)**: 1821 (9.5%)
* **Genes retained**: 17,384
* **Outlier threshold**: 3 SD total expression (D-009)
* **Samples dropped (outliers)**: 15 (0.9%)

## Metabolomics

* **Input shape**: 928 × 225
* **Missing values**: 0 (pre-imputed at source)
* **Zero-variance metabolites**: 0 — all 225 retained
* **Outlier threshold**: 3 SD total metabolite level (D-010)
* **Samples dropped (outliers)**: 11 (1.2%)

## Paired Dataset

* **Pre-QC paired N**: 912
* **Post-QC paired N**: 898
* **Outlier overlap**: 0 — independent quality issues per modality
* **Status**: PASS (above 800 floor)

## Notes

* Metabolomics column headers are metabolite names, not compound IDs —
  name → KEGG compound ID translation required in week 2
* Transcriptomics column headers had Entrez IDs appended — stripped during QC
* CCLE_ID column dropped from metabolomics (all NaN after index set on DepMap ID)

## Exploration Flags

* Metabolite-metabolite correlation matrix shows clear block structure —
  co-regulated metabolite clusters exist. Supports KEGG constraint design.
* Transcriptomics intrinsic dimensionality is high — 50 PCs explain ~56% variance.
  80% threshold not reached within 50 PCs.
* Metabolomics intrinsic dimensionality is lower — 50 PCs explain ~82% variance.
  No visible cancer type clustering in PC1/PC2.
* Modality asymmetry confirmed — symmetric encoder architecture may underserve
  transcriptomics. Flag for preregistration.
* Cross-modality sample-level mean correlation: r=0.028, p=0.398 — not significant.
  Explained by aggregation artifact, not a biological finding.
* Known gene-metabolite pairs: 2/4 detectable (LDHA/aconitate r=0.124,
  SHMT2/3-phosphoglycerate r=-0.155). Weak signal present. Supports KEGG prior.