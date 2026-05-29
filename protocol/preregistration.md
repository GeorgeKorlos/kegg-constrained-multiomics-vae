# P2 — KEGG-Constrained Multi-Omics VAE

---

## Section 1 — Falsifiable Claim

### Primary claim

A VAE with KEGG module membership encoded as a structural prior on the latent space
produces latent representations with measurably higher pathway-level coherence than
a VAE trained on the same data without the KEGG prior, and this difference in
coherence translates to a measurable performance advantage in the downstream P3 graph
task.

Two sub-claims, each independently falsifiable:

**Claim A — Biological coherence:** In the KEGG-constrained VAE, genes and metabolites
assigned to the same KEGG module will have higher cosine similarity in embedding space
than genes and metabolites assigned to different modules. The unconstrained baseline
will not reproduce this block structure at the same level.

**Claim B — Downstream utility:** P3's GNN, initialized with KEGG-constrained
embeddings, will outperform the same GNN initialized with unconstrained VAE embeddings
on its primary evaluation metric (held fixed in P3's own preregistration).

### What falsifies the claim

| Outcome | Verdict |
|---|---|
| Unconstrained baseline matches or beats KEGG-VAE on pathway coherence | Claim A falsified — prior adds no structure |
| KEGG-VAE coherence improves but P3 performance does not | Claim B falsified — structural prior does not transfer to task |
| Per-block variance collapses (Guard 2 fires) | Model failed to learn; result is uninterpretable, not falsifying |
| Either modality fails to beat mean predictor (Guard 3 fires) | Reconstruction failed; result is uninterpretable, not falsifying |
| λ sweep shows no monotonic relationship between constraint strength and coherence | Suggests KEGG regularization is not the active ingredient |

Guard-triggered failures are halts, not falsifications. The claim is only tested on a model
that passes all four training guards.

### What is not claimed

- That the KEGG-constrained VAE produces better reconstruction than the baseline
  (it may not — the constraint trades reconstruction accuracy for biological structure)
- That embeddings generalize to primary tissue (training data is cancer cell lines)
- That the architecture is the methodological contribution of the portfolio (it is not —
  see P3, P5)

---

## Section 2 — Loss Function and Training Guards

### 2.1 Loss definition

```
L_total = L_recon_trans + L_recon_meta + β · L_KL + λ(t) · L_KEGG
```

**Reconstruction — transcriptomics**

```
L_recon_trans = (1 / F_trans) · ||x_trans - x̂_trans||²_F
```

Where F_trans = 17,384 (genes retained post-QC). Normalized by feature count to
prevent transcriptomics from dominating the scalar loss purely by dimensionality.
MSE over the standardized input — no sigmoid, inputs are not bounded.

**Reconstruction — metabolomics**

```
L_recon_meta = (1 / F_meta) · ||x_meta - x̂_meta||²_F
```

Where F_meta = 225. Same normalization convention.

**KL divergence**

```
L_KL = -½ · (1/d) · Σ_i [1 + log σ²_i − μ²_i − σ²_i]
```

Where d = 128 (latent dimension). Normalized by d. Computed on the joint posterior
from the product-of-experts fusion of both modality encoders.

**KEGG regularization**

```
L_KEGG = z_k ∈ ℝ^(B × b_k)
```

Where z_k ∈ ℝ^(B × b) is the sub-vector of the latent corresponding to block k
(b = 128/K), Cov_batch is the empirical covariance matrix computed over the batch
dimension, and K is the number of KEGG modules retained after week 2 filtering
(value TBD — logged in `reports/kegg_coverage.md` after Day 4).

This term penalizes cross-block covariance. Minimizing it encourages each block to
capture variance independently. The KEGG module partition defines which latent
dimensions correspond to which biological modules; the regularization enforces that
they remain separable.

**Modality weighting and gradient balancing — Guard 4**

The raw loss sum weights transcriptomics 17,384/225 ≈ 77× more than metabolomics by
feature count even after per-feature normalization above. Guard 4 corrects this at
the gradient level.

Let η_t = ||∇_θ L_recon_trans||₂ / ||∇_θ L_recon_meta||₂ computed per batch over
shared encoder parameters θ. The effective per-batch loss used for the backward pass
is:

```
L_recon_eff = (1 / η_t) · L_recon_trans + L_recon_meta    [if η_t > 1]
L_recon_eff = L_recon_trans + η_t · L_recon_meta           [if η_t < 1]
```

This normalizes gradient contributions so neither modality systematically dominates
parameter updates. η_t is computed before the backward pass, not applied as a static
weight. Logged per epoch for inspection.

η_t is not a hyperparameter — it is a diagnostic and a corrective computed from
the data each batch. It is not tuned.

### 2.2 λ schedule — Guard 1

λ(t) is not epoch-indexed. It is gated on reconstruction plateau.

```
λ(t) = λ_max · min(1, t_stable / T_warm)
```

Where t_stable is the number of epochs since reconstruction loss (both modalities)
ceased improving by more than δ = 0.001 on the validation set, and T_warm is a
fixed warm-up duration (default: 10 plateau-epochs, set in `config/training.yaml`).

**Hard constraint:** λ must not exceed 0.5 · λ_max before Guard 3 passes. The KEGG
regularization is not activated until reconstruction has stabilized — the model must
first learn to reconstruct before it is forced to restructure the latent space.

λ_max is a hyperparameter swept in the ablation (Section 5).

### 2.3 Four training guards

All four guards are checked at the end of every epoch. Any failure triggers a halt
and logs the epoch, guard ID, and current metric value. A model that halted on a
guard does not enter evaluation.

**Guard 1 — λ ramp gate**

λ(t) is frozen at 0 until reconstruction plateau is confirmed. If reconstruction has
not plateaued within 50 epochs, log a warning and begin λ ramp anyway to prevent
indefinite deferral. This is a safety valve, not the intended path.

Implementation: in the training loop, track `epochs_since_recon_improvement` for
both modalities. Ramp begins when both have stalled for T_warm epochs.

Config keys (config/training.yaml):
```yaml
lambda_warmup:
  T_warm: 10
  delta: 0.001
  max_deferral_epochs: 50
```

**Guard 2 — Per-block variance**

At each epoch, compute Var(z_k) for each of the K blocks over the full validation set.
Compute the scale-normalized failure condition:

```
Failure: Var(z_k) < α · E_k[Var(z_k)]  for p consecutive epochs
```

Default: α = 0.1, p = 5. A block with variance below 10% of the mean block variance
for 5 consecutive epochs is dead. Halt training.

Config keys:
```yaml
block_variance_guard:
  alpha: 0.1
  p: 5
```

This guard fires before the model is evaluated. A block-collapsed model does not
produce interpretable embeddings regardless of reconstruction performance.

**Guard 3 — Per-modality reconstruction baseline**

Compute R² for each modality against the mean-predictor baseline (predicting the
training mean for every sample) on a held-out validation set. Evaluated after the
λ ramp has been active for at least 5 epochs.

```
Failure: R²_trans < 0  OR  R²_meta < 0
```

A negative R² means the model predicts worse than the mean. This is a categorical
failure. Pearson r thresholds were considered and rejected — they conflate ranking
with magnitude and are sensitive to distribution shape. R² vs mean predictor is
distribution-free and directly interpretable.

R² is computed on standardized inputs (training mean = 0 by construction, so the
mean predictor predicts all-zeros, and R² = 1 − MSE_model / MSE_null where
MSE_null = variance of the standardized validation targets).

**Guard 4 — Dynamic gradient norm balancing**

Logged every batch. No failure condition — Guard 4 is a corrective, not a monitor.
If η_t > 100 or η_t < 0.01 for more than 20 consecutive batches, emit a warning
(not a halt) and flag the run for inspection. Values this extreme indicate a
modality is contributing negligible gradient signal.

---

## Section 3 — Tensor-Level Data Flow

Shapes use: B = batch size, F_t = 17,384, F_m = 225, d = 128, K = TBD (post week 2),
b_k = block size for module k, Σ b_k = 128, b_k ≥ 2, allocated proportionally to module gene count
Block sizes vary by module to reflect KEGG module size. Block size table logged in reports/kegg_coverage.md.

```
INPUT
  x_trans:  (B, 17384)   standardized log1p(TPM), float32
  x_meta:   (B, 225)     standardized log-metabolite levels, float32

TRANSCRIPTOMICS ENCODER
  (B, 17384)
  → Linear + LayerNorm + GELU  [dim: h_t, TBD — pending OPEN-005]
  → Linear + LayerNorm + GELU  [dim: h_t/2]
  → Linear                     [dim: 2·128]
  → split → μ_t: (B, 128),  log_σ²_t: (B, 128)

METABOLOMICS ENCODER
  (B, 225)
  → Linear + LayerNorm + GELU  [dim: h_m, TBD — pending OPEN-005]
  → Linear                     [dim: 2·128]
  → split → μ_m: (B, 128),  log_σ²_m: (B, 128)

PRODUCT-OF-EXPERTS FUSION
  σ²_joint_i = 1 / (1/σ²_t_i + 1/σ²_m_i + 1/σ²_prior_i)
  μ_joint_i  = σ²_joint_i · (μ_t_i/σ²_t_i + μ_m_i/σ²_m_i)
  where σ²_prior = 1 (standard normal prior per dimension)

  Output: μ_joint: (B, 128),  σ²_joint: (B, 128)

REPARAMETERIZATION
  ε ~ N(0, I)  shape (B, 128)
  z = μ_joint + ε · sqrt(σ²_joint)
  z: (B, 128)

BLOCK PARTITION
  z → [z_1 | z_2 | ... | z_K]
  z_k: (B, b)  where b = 128/K
  Partition is a deterministic slice — no learned parameters at this step.
  Soft assignment weights w_{gk} (gene g → module k) are applied in the decoder,
  not at the partition step.

TRANSCRIPTOMICS DECODER
  z: (B, 128)
  Soft assignment gate per gene g:
    ĝ_g = Σ_k w_{gk} · (W_k_dec · z_k)    [weighted sum over K block decoders]
    w_{gk} from G × K soft assignment matrix, row-normalized (softmax over k)
    initialized from binary KEGG membership, learned during training
  → Linear + LayerNorm + GELU  [dim: h_t, symmetric to encoder — pending OPEN-005]
  → Linear                     [dim: 17384]
  x̂_trans: (B, 17384)

METABOLOMICS DECODER
  z: (B, 128)
  Same soft assignment structure over M × K metabolite-module membership
  → Linear + LayerNorm + GELU  [dim: h_m]
  → Linear                     [dim: 225]
  x̂_meta: (B, 225)

LOSS COMPUTATION
  L_recon_trans: scalar  (normalized MSE)
  L_recon_meta:  scalar  (normalized MSE)
  L_KL:          scalar  (on joint posterior)
  L_KEGG:        scalar  (cross-block covariance penalty)
  η_t:           scalar  (gradient norm ratio, computed pre-backward)
  L_total:       scalar

EXPORT (write time — not at training time)
  z_gene_g = pre-aggregation decoder activation for gene g
           = Σ_k w_{gk} · (W_k_dec · z_k)   [before final projection to F_t]
           shape per gene: (128,)
  → L2 normalize → (128,)

  z_metabolite_m = equivalent pre-aggregation activation for metabolite m
  → L2 normalize → (128,)

  HGNC → UniProt translation applied at write time.
  Output: HDF5 per locked schema.
```

**Open dependency:** b = 128/K requires 128 % K == 0. K is determined after Day 4.
If K does not divide 128, options are: (a) pad latent to next multiple of K with
zero-masked dimensions, (b) use unequal block sizes with the constraint that
Σ b_k = 128. Option selection deferred to Day 4 — must be resolved before model
code is written.

---

## Section 4 — Evaluation Hierarchy

Steps are sequential. A model failing any step does not proceed to later steps.
Guards (Section 2.3) are prerequisites for Step 0.

### Step 0 — Reconstruction sanity (mandatory pass/fail)

**Metric:** R²_trans and R²_meta on held-out validation set vs mean predictor.

**Pass condition:** R²_trans > 0 AND R²_meta > 0

**Failure action:** stop. The model cannot be evaluated. Diagnose: check Guard 2
(block collapse), check η_t logs (modality dominance), check λ schedule (premature
KEGG activation).

### Step 1 — Latent space health (mandatory pass/fail)

**Metrics:**
- Var(z_k) distribution across K blocks — all blocks must be alive (Guard 2 threshold)
- η_t mean and variance over training — flagged if extreme (Guard 4)
- KL per dimension — check for posterior collapse (dimensions with KL ≈ 0 throughout)

**Pass condition:** No block-dead condition. KL active in > 80% of latent dimensions.

**Failure action:** stop. Diagnose: increase β (if KL collapse), reduce λ_max
(if block collapse coincides with KEGG term activation), check decoder capacity.

### Step 2 — KEGG constraint effect (primary comparison)

**Models compared:** KEGG-constrained VAE (condition A) vs unconstrained VAE (condition B).
Both trained under identical conditions except λ = 0 in condition B.

**Metrics:**

*Pathway coherence score:* For each KEGG module k, compute mean pairwise cosine
similarity of embeddings for genes/metabolites assigned to module k. Compare
in-module vs out-of-module similarity. Report:
```
coherence_k = mean_cosine(genes in module k) − mean_cosine(genes in random set of same size)
```
Aggregate over all K modules. Higher is better. Test: paired t-test over modules
between condition A and B.

*Reconstruction MSE:* KEGG-constrained VAE is expected to have equal or slightly
higher MSE than unconstrained (the constraint is a regularizer that trades
reconstruction for structure). If constrained has substantially higher MSE (>20%),
flag as a capacity problem.

**Pass condition for Claim A:** mean coherence_k(A) > mean coherence_k(B), p < 0.05
(paired t-test, modules as units of analysis).

### Step 3 — Biological validity (secondary)

**Metrics:**

*Pathway activity score correlation:* For each of four pathways (glycolysis,
serine biosynthesis, TCA, one-carbon metabolism), compute:
- Gene pathway activity = mean standardized expression of pathway genes
- Metabolite pathway activity = mean standardized level of pathway metabolites
  present in the 225-metabolite panel
- Pearson r and Spearman ρ between gene and metabolite pathway activity scores

*Block alignment:* Check whether the latent block assigned to module k has higher
mutual information with the pathway activity score of module k than with scores of
other modules. This is a posterior check that the block partition is doing what the
KEGG prior intends.

**Pass condition:** Pearson r > 0 for at least 3 of 4 pathways in the constrained
model. No explicit threshold on block alignment — reported descriptively.

**Failure implication:** If pathway activity correlations are zero in both constrained
and unconstrained models, the signal is too weak to validate the biological claim.
This is a known risk given the weak cross-modality correlation (r = 0.028 at global
level, weak signal detectable only at pair level). The ablation (Section 5) is
designed to diagnose this scenario.

### Step 4 — Downstream task performance (tertiary, P3-dependent)

**Condition:** P3's GNN is trained with embeddings from condition A, then from
condition B. P3's primary metric is fixed in P3's own preregistration — P2 does not
define it.

**Pass condition for Claim B:** P3 metric(A) > P3 metric(B). No threshold; the
direction of the effect and its confidence interval are the finding.

**Note:** Step 4 cannot be completed until P3 is implemented. It is listed here for
completeness and because the portfolio claim depends on it. P2's internal evaluation
stops at Step 3.

---

## Section 5 — Ablation Matrix

All ablations share the same data split, optimizer, and hyperparameters except the
varied condition. A single train/val/test split (80/10/10) is fixed before any
training and held constant across all conditions.

| ID | Label | KEGG constraint | Assignment | λ | Purpose |
|----|-------|----------------|------------|---|---------|
| A | Constrained-soft | Yes | Soft (learned) | λ_max | Primary model |
| B | Unconstrained | No | — | 0 | Falsification baseline |
| C | Constrained-hard | Yes | Hard (binary KEGG) | λ_max | Soft vs hard comparison |
| D | Half-λ | Yes | Soft | λ_max / 2 | λ sensitivity — upper |
| E | Quarter-λ | Yes | Soft | λ_max / 4 | λ sensitivity — lower |

**λ_max** is selected by held-out validation: the largest value for which condition A
passes Step 0 and Step 1. This is determined empirically before the ablation is run.
The sweep for λ_max selection is over {0.01, 0.05, 0.1, 0.5, 1.0} — log-spaced.
Results of the selection sweep are reported but are not the primary ablation.

**Conditions C, D, E** are run only if condition A passes Step 2. If A fails Step 2,
the claim is falsified and ablation comparisons are moot.

**Decision rule:**
- A > B on coherence: proceed to Step 3 and Step 4
- A ≤ B on coherence: claim falsified; report result as null finding; do not run C/D/E
- A > C: soft assignment adds value over hard — supports the learned prior claim
- A ≈ C: soft assignment confers no benefit — the KEGG structure alone is sufficient;
  consider hard masking as the cleaner implementation for P3 hand-off
- Monotone decrease in coherence from A → D → E: confirms λ is the active ingredient

---

## Section 6 — Empirical Observations Informing Architecture

(From week 1 exploration — these motivated the design choices above.)

- Metabolite co-regulation block structure confirmed in correlation matrix. KEGG module
  partition of latent space is empirically motivated, not only theoretically justified.
- Transcriptomics intrinsic dimensionality substantially higher than metabolomics
  (>50 PCs for 80% variance vs ~48 PCs). Asymmetric encoder depth is under
  consideration (OPEN-005). Product-of-experts fusion is architecture-agnostic to
  encoder depth — this choice is robust to OPEN-005 resolution.
- Cross-modality signal is weak at the global level (r = 0.028, p = 0.398), detectable
  at gene-metabolite pair level (LDHA/aconitate r = 0.124, SHMT2/3-phosphoglycerate
  r = −0.155). Joint embedding is justified but the KEGG structural prior is
  load-bearing. This is why Step 3 uses pathway-level scores rather than global
  cross-modality correlation.
- Weak cross-modality signal reduces expected downstream advantage — the ablation
  design (Step 4, condition A vs B) is designed to isolate this contribution
  independently of the correlation strength.

---

## Section 7 — Planned Diagnostics

- Pathway activity scores (glycolysis, serine biosynthesis, TCA, one-carbon
  metabolism) correlated across modalities as post-hoc validation of KEGG constraint.
  Gene score = mean expression of pathway genes; metabolite score = mean of pathway
  metabolites in panel.
- Spearman ρ as robustness check on Pearson results — less sensitive to outliers given
  the tight metabolomics value range.
- Mutual information between modalities as non-linear signal check (upper bound on
  linear correlation findings).
- Stratified analysis by cancer lineage (lung, breast, melanoma) to test whether
  cross-modality signal emerges after removing inter-lineage heterogeneity.
- η_t trajectory plots over training — visual confirmation that gradient balancing
  prevents modality dominance.
- Block variance heatmap (K × epochs) — visual confirmation Guard 2 is not close
  to firing in the primary run.

---

## Open decisions at time of preregistration

| ID | Decision | Gate |
|----|----------|------|
| OPEN-005 | Asymmetric vs symmetric encoder depth | Before model code |
| — | b = 128/K divisibility check | After Day 4 (K confirmed) |
| — | λ_max selection | Before ablation runs, via validation sweep |

No model code is written until this document is committed and OPEN-005 is resolved.

---
