import sys
import torch
import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader

sys.path.insert(0, ".")
PROCESSED_DIR = "data/processed/"


class MultiOmicsDataset(Dataset):
    def __init__(self, processed_dir=PROCESSED_DIR):
        self.trans = np.load(f"{processed_dir}transcriptomics.npy")
        self.meta = np.load(f"{processed_dir}metabolomics.npy")
        self.sample_ids = pd.read_csv(f"{processed_dir}sample_ids.csv")[
            "sample_id"
        ].tolist()

        assert (
            self.trans.shape[0] == self.meta.shape[0]
        ), "Transcriptomics and metabolomics sample counts do not match"
        assert self.trans.shape[0] == len(
            self.sample_ids
        ), "Sample ID count does not match tensor row count"

    def __len__(self):
        return self.trans.shape[0]

    def __getitem__(self, idx):
        return (
            torch.tensor(self.trans[idx], dtype=torch.float32),
            torch.tensor(self.meta[idx], dtype=torch.float32),
        )

    @property
    def trans_dim(self):
        return self.trans.shape[1]

    @property
    def meta_dim(self):
        return self.meta.shape[1]


if __name__ == "__main__":
    import torch
    from torch.utils.data import DataLoader

    dataset = MultiOmicsDataset()

    print(f"Transcriptomics tensor shape: {dataset.trans.shape}")
    print(f"Metabolomics tensor shape: {dataset.meta.shape}")
    print(f"Sample IDs count: {len(dataset.sample_ids)}")
    print(f"Transcriptomics dtype: {dataset.trans.dtype}")
    print(f"Metabolomics dtype: {dataset.meta.dtype}")

    print(
        f"\nTranscriptomics — mean: {dataset.trans.mean():.4f}, std: {dataset.trans.std():.4f}"
    )
    print(
        f"Metabolomics — mean: {dataset.meta.mean():.4f}, std: {dataset.meta.std():.4f}"
    )

    print(f"\nSpot-check — 5 samples:")
    for i in range(5):
        t, m = dataset[i]
        print(
            f"  Sample {dataset.sample_ids[i]}: "
            f"trans={t.shape}, range=[{t.min():.2f}, {t.max():.2f}] | "
            f"meta={m.shape}, range=[{m.min():.2f}, {m.max():.2f}]"
        )

    loader = DataLoader(dataset, batch_size=32, shuffle=True)
    batch_trans, batch_meta = next(iter(loader))
    print(f"\nDataLoader batch — transcriptomics: {batch_trans.shape}")
    print(f"DataLoader batch — metabolomics: {batch_meta.shape}")

    passed = (
        dataset.trans.shape == (898, 17384)
        and dataset.meta.shape == (898, 225)
        and dataset.trans.dtype == np.float32
        and dataset.meta.dtype == np.float32
        and abs(dataset.trans.mean()) < 0.01
        and abs(dataset.meta.mean()) < 0.01
    )
    print(f"\n{'PASS' if passed else 'FAIL — check output above.'}")
