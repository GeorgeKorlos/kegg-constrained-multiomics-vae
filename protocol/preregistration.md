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