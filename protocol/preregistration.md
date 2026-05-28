## Empirical observations informing architecture (from exploration)

* Metabolite co-regulation block structure confirmed in correlation matrix —
  KEGG module partition of latent space is empirically motivated, not only
  theoretically justified
* Transcriptomics intrinsic dimensionality substantially higher than metabolomics —
  encoder architecture must account for asymmetry
* Cross-modality signal is weak at the global level, detectable at the
  gene-metabolite pair level — joint embedding is justified but the KEGG
  structural prior is load-bearing, not decorative
* Weak cross-modality signal reduces expected downstream node initialization
  advantage — ablation design must isolate this contribution explicitly

## Planned diagnostics — cross-modality signal validation

* Pathway activity scores (glycolysis, serine biosynthesis, TCA, one-carbon
  metabolism) correlated across modalities as post-hoc validation of KEGG
  constraint. Gene score = mean expression of pathway genes, metabolite score =
  mean of pathway metabolites in panel.
* Spearman correlation as robustness check on Pearson results — less sensitive
  to outliers given the tight metabolomics value range.
* Mutual information between modalities as non-linear signal check.
* Stratified analysis by cancer lineage (lung, breast, melanoma) to test whether
  cross-modality signal emerges after removing inter-lineage heterogeneity.