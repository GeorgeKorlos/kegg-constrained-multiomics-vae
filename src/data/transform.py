import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
import pickle
import os

TRANSCRIPTOMICS_PATH = "data/raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
METABOLOMICS_PATH = "data/raw/CCLE_metabolomics_20190502.csv"
PROCESSED_DIR = "data/processed/"

VAR_THRESHOLD = 0.05
OUTLIER_SD = 3


def load_and_filter_transcriptomics(path=TRANSCRIPTOMICS_PATH):
    df = pd.read_csv(path, index_col=0)
    df.index.name = "DepMap_ID"
    df.columns = df.columns.str.replace(r"\s*\(\d+\)\s*", "", regex=True)

    gene_var = df.var(axis=0)
    df = df.loc[:, gene_var >= VAR_THRESHOLD]

    sample_sum = df.sum(axis=1)
    outlier_mask = np.abs(stats.zscore(sample_sum)) > OUTLIER_SD
    df = df[~outlier_mask]

    return df


def load_and_filter_metabolomics(path=METABOLOMICS_PATH):
    df = pd.read_csv(path)
    df = df.set_index("DepMap_ID")
    df = df.drop(columns=["CCLE_ID"])

    sample_sum = df.sum(axis=1)
    outlier_mask = np.abs(stats.zscore(sample_sum)) > OUTLIER_SD
    df = df[~outlier_mask]

    return df


def build_paired_dataset(save=True):
    trans = load_and_filter_transcriptomics()
    meta = load_and_filter_metabolomics()

    paired = trans.join(meta, how="inner")
    if paired.shape[0] < 800:
        raise ValueError(f"Paired N {paired.shape[0]} below 800 floor — stop.")

    trans_paired = paired[trans.columns]
    meta_paired = paired[meta.columns]

    sample_ids = paired.index.tolist()
    gene_names = trans.columns.tolist()
    metabolite_names = meta.columns.tolist()

    trans_scaler = StandardScaler()
    meta_scaler = StandardScaler()

    trans_scaled = trans_scaler.fit_transform(trans_paired.values)
    meta_scaled = meta_scaler.fit_transform(meta_paired.values)

    if save:
        os.makedirs(PROCESSED_DIR, exist_ok=True)

        np.save(f"{PROCESSED_DIR}transcriptomics.npy", trans_scaled.astype(np.float32))
        np.save(f"{PROCESSED_DIR}metabolomics.npy", meta_scaled.astype(np.float32))

        pd.Series(sample_ids, name="sample_id").to_csv(
            f"{PROCESSED_DIR}sample_ids.csv", index=False
        )
        pd.Series(gene_names).to_csv(f"{PROCESSED_DIR}gene_names.csv", index=False)
        pd.Series(metabolite_names).to_csv(
            f"{PROCESSED_DIR}metabolite_names.csv", index=False
        )

        with open(f"{PROCESSED_DIR}transcriptomics_scaler.pkl", "wb") as f:
            pickle.dump(trans_scaler, f)
        with open(f"{PROCESSED_DIR}metabolomics_scaler.pkl", "wb") as f:
            pickle.dump(meta_scaler, f)

        print(f"Saved processed tensors to {PROCESSED_DIR}")
        print(f"Transcriptomics: {trans_scaled.shape}")
        print(f"Metabolomics: {meta_scaled.shape}")
        print(f"Paired N: {len(sample_ids)}")

    return trans_scaled.astype(np.float32), meta_scaled.astype(np.float32), sample_ids


if __name__ == "__main__":
    build_paired_dataset(save=True)
