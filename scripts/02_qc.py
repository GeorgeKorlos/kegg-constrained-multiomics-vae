import numpy as np
import pandas as pd
from scipy import stats

### Transcriptomics

transcriptomics_path = "data/raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv"

transcriptomics = pd.read_csv(transcriptomics_path, index_col=0)
transcriptomics.index.name = "DepMap_ID"

transcriptomics.columns = transcriptomics.columns.str.replace(
    r"\s*\(\d+\)\s*", "", regex=True
)

# Gene-level stats
gene_mean = transcriptomics.mean(axis=0)
gene_var = transcriptomics.var(axis=0)
zero_fraction = (transcriptomics == 0).mean(axis=0)

cutoff_10 = gene_var.quantile(0.10)
low_var_genes = gene_var[gene_var < cutoff_10].index

# Sample-level stats
sample_sum_tx = transcriptomics.sum(axis=1)
sample_mean_tx = transcriptomics.mean(axis=1)
sample_std_tx = transcriptomics.std(axis=1)

sample_sum_z_tx = np.abs(stats.zscore(sample_sum_tx))
outlier_mask_tx = sample_sum_z_tx > 3
outlier_samples_tx = transcriptomics.loc[outlier_mask_tx]

print(f"Transcriptomics missing values: {transcriptomics.isna().sum().sum()}")
print(f"Gene mean expression:\n{gene_mean.describe()}")
print(f"Gene variance:\n{gene_var.describe()}")
print(f"10th percentile variance cutoff: {cutoff_10:.4f}")
print(f"Low variance genes (below 10th pct): {len(low_var_genes)}")
print(f"Sample total expression:\n{sample_sum_tx.describe()}")
print(f"Outlier samples (>3 SD total expression): {outlier_mask_tx.sum()}")

for threshold in [0.01, 0.05, 0.10, 0.20, 0.50]:
    n = (gene_var < threshold).sum()
    print(f"Var < {threshold}: {n} genes ({n/len(gene_var)*100:.1f}%)")

### Metabolomics

metabolomics_path = "data/raw/CCLE_metabolomics_20190502.csv"
metabolomics = pd.read_csv(metabolomics_path)
metabolomics = metabolomics.set_index("DepMap_ID")
metabolomics = metabolomics.apply(pd.to_numeric, errors="coerce")
metabolomics = metabolomics.drop(columns=["CCLE_ID"])

metabolite_mean = metabolomics.mean(axis=0)
metabolite_var = metabolomics.var(axis=0)
metabolite_min = metabolomics.min(axis=0)
metabolite_max = metabolomics.max(axis=0)

zero_var_metabolites = metabolite_var[metabolite_var == 0].index

# Sample-level stats
sample_sum_meta = metabolomics.sum(axis=1)
sample_mean_meta = metabolomics.mean(axis=1)
sample_std_meta = metabolomics.std(axis=1)

sample_sum_z_meta = np.abs(stats.zscore(sample_sum_meta))
outlier_mask_meta = sample_sum_z_meta > 3
outlier_samples_meta = metabolomics.loc[outlier_mask_meta]

print(f"Metabolomics shape after dropping CCLE_ID: {metabolomics.shape}")
print(f"Metabolomics missing values: {metabolomics.isna().sum().sum()}")
print(f"Metabolite mean:\n{metabolite_mean.describe()}")
print(f"Metabolite variance:\n{metabolite_var.describe()}")
print(f"Zero-variance metabolites: {len(zero_var_metabolites)}")
print(f"Sample total expression:\n{sample_sum_meta.describe()}")
print(f"Outlier samples (>3 SD total expression): {outlier_mask_meta.sum()}")

print(metabolomics.columns[metabolomics.isna().any()].tolist())
print(metabolomics.shape)
print(metabolomics.head())

trans_outlier_ids = transcriptomics.index[outlier_mask_tx]
meta_outlier_ids = metabolomics.index[outlier_mask_meta]
overlap = trans_outlier_ids.intersection(meta_outlier_ids)

print(f"Transcriptomics-only outliers: {len(trans_outlier_ids) - len(overlap)}")
print(f"Metabolomics-only outliers: {len(meta_outlier_ids) - len(overlap)}")
print(f"Overlap (both modalities): {len(overlap)}")
print(f"Total unique samples to drop: {len(trans_outlier_ids.union(meta_outlier_ids))}")

### Paired dataset

trans_clean = transcriptomics[~outlier_mask_tx]
meta_clean = metabolomics[~outlier_mask_meta]

paired_clean = trans_clean.join(meta_clean, how="inner")
print(f"Pre-QC paired N: 912")
print(f"Post-QC paired N: {paired_clean.shape[0]}")
print(f"PASS" if paired_clean.shape[0] >= 800 else "FAIL — below 800 floor")


transcriptomics_gene_qc = pd.DataFrame(
    {"mean": gene_mean, "variance": gene_var, "zero_fraction": zero_fraction}
)

transcriptomics_gene_qc.to_csv("reports/transcriptomics_gene_qc.csv")

transcriptomics_sample_qc = pd.DataFrame(
    {
        "total_expression": sample_sum_tx,
        "mean_expression": sample_mean_tx,
        "std_expression": sample_std_tx,
        "outlier": outlier_mask_tx,
    }
)

transcriptomics_sample_qc.to_csv("reports/transcriptomics_sample_qc.csv")


metabolomics_metabolite_qc = pd.DataFrame(
    {
        "mean": metabolite_mean,
        "variance": metabolite_var,
        "min": metabolite_min,
        "max": metabolite_max,
    }
)

metabolomics_metabolite_qc.to_csv("reports/metabolomics_metabolite_qc.csv")

metabolomics_sample_qc = pd.DataFrame(
    {
        "total_expression": sample_sum_meta,
        "mean_expression": sample_mean_meta,
        "std_expression": sample_std_meta,
        "outlier": outlier_mask_meta,
    }
)

metabolomics_sample_qc.to_csv("reports/metabolomics_sample_qc.csv")
