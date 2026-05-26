import hashlib
import pandas as pd

transcriptomics_path = "data/raw/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
metabolomics_path = "data/raw/CCLE_metabolomics_20190502.csv"


def get_sha256(filepath):
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


print(f"Transcriptomics SHA256: {get_sha256(transcriptomics_path)}")
print(f"Metabolomics SHA256: {get_sha256(metabolomics_path)}")

trans_df = pd.read_csv(transcriptomics_path, index_col=0)
trans_df.index.name = "DepMap_ID"
meta_df = pd.read_csv(metabolomics_path)
meta_df = meta_df.set_index("DepMap_ID")

print(f"Transcriptomics shape: {trans_df.shape}")
print(f"Metabolomics shape: {meta_df.shape}")

print("\n")
print(trans_df.head(5))
print("\n")
print(meta_df.head(5))

joined_df = trans_df.join(meta_df, how="inner")

print(f"Paired N: {joined_df.shape[0]}")
