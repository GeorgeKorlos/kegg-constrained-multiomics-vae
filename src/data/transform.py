import os
import json
import pickle

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler

TRANSCRIPTOMICS_PATH = "data/raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
METABOLOMICS_PATH = "data/raw/CCLE_metabolomics_20190502.csv"
PROCESSED_DIR = "data/processed/"
SPLIT_PATH = "data/splits/split_v1.json"

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


def _load_train_indices(split_path, expected_n):
    with open(split_path, "r") as f:
        split = json.load(f)
    if split["n"] != expected_n:
        raise ValueError(
            f"Split n={split['n']} != paired N={expected_n}. "
            f"Split was built against a different dataset — refusing to proceed."
        )
    return split["train"]


def build_paired_dataset(save=True, split_path=SPLIT_PATH):
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

    existing_ids_path = f"{PROCESSED_DIR}sample_ids.csv"
    if os.path.exists(existing_ids_path):
        prev_ids = pd.read_csv(existing_ids_path)["sample_id"].tolist()
        if prev_ids != sample_ids:
            raise RuntimeError(
                "Paired row ordering changed vs committed sample_ids.csv. "
                "Positional split indices are no longer valid. STOP — do not "
                "refit the scaler against a reordered dataset."
            )

    train_idx = _load_train_indices(split_path, expected_n=len(sample_ids))

    trans_scaler = StandardScaler()
    meta_scaler = StandardScaler()

    trans_scaler.fit(trans_paired.values[train_idx])
    meta_scaler.fit(meta_paired.values[train_idx])

    trans_scaled = trans_scaler.transform(trans_paired.values)
    meta_scaled = meta_scaler.transform(meta_paired.values)

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
        print(f"Scaler fit on TRAIN ONLY: {len(train_idx)} rows")
        # Sanity: train rows should be ~N(0,1); val/test will drift slightly.
        print(
            f"Train-row mean/std (trans): "
            f"{trans_scaled[train_idx].mean():.4f} / {trans_scaled[train_idx].std():.4f}"
        )

    return trans_scaled.astype(np.float32), meta_scaled.astype(np.float32), sample_ids


if __name__ == "__main__":
    build_paired_dataset(save=True)
