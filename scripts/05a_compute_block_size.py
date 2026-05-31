import numpy as np
import pandas as pd

GENE_MODULE_PATH = "data/processed/gene_module_matrix.npy"
METABOLITE_MODULE_PATH = "data/processed/metabolite_module_matrix.npy"
MODULE_IDS_PATH = "data/processed/module_ids.csv"
MODULE_METADATA_PATH = "data/raw/kegg_module_metadata.tsv"
OUT_PATH = "data/processed/block_sizes.csv"

genes = np.load(GENE_MODULE_PATH)
metabolites = np.load(METABOLITE_MODULE_PATH)

ids = pd.read_csv(MODULE_IDS_PATH)
metadata = pd.read_csv(MODULE_METADATA_PATH, sep="\t")

required = {"module_id", "module_name"}
missing = required - set(metadata.columns)

if missing:
    raise ValueError(f"Missing metadata columns: {missing}")

module_list = ids.iloc[:, 0].tolist()

gene_nonzero = genes.sum(axis=0) > 0
met_nonzero = metabolites.sum(axis=0) > 0

both_mask = gene_nonzero & met_nonzero
both_indices = np.where(both_mask)[0]
both_modules = np.array(module_list)[both_mask]

target_modules = ["M00001", "M00009", "M00020", "M00141"]
retained_module_ids = list(both_modules)

print("\n=== Four-pathway coverage check ===\n")

for mod_id in target_modules:
    if mod_id not in retained_module_ids:
        print(f"{mod_id}: NOT in retained {len(retained_module_ids)} modules")
        continue

    idx = retained_module_ids.index(mod_id)
    global_idx = both_indices[idx]

    n_genes = int(genes[:, global_idx].sum())
    n_metas = int(metabolites[:, global_idx].sum())

    flag = " thin (<4 metabolites)" if n_metas < 4 else ""
    print(f"{mod_id}: genes={n_genes}, metabolites={n_metas}{flag}")

s_k = genes.sum(axis=0)
both_s_k = s_k[both_mask]

s_k_metabolites = metabolites.sum(axis=0)
both_s_k_metabolites = s_k_metabolites[both_mask]

S = both_s_k.sum()

b_k_raw = 128.0 * both_s_k / S
b_k = np.maximum(2, np.round(b_k_raw)).astype(int)

diff = b_k.sum() - 128

while diff > 0:
    candidates = np.where(b_k > 2)[0]

    if len(candidates) == 0:
        raise ValueError("Cannot reduce further without violating min-floor")

    idx = candidates[np.argmax(b_k[candidates])]
    b_k[idx] -= 1
    diff -= 1

while diff < 0:
    idx = np.argmax(b_k)
    b_k[idx] += 1
    diff += 1

module_name_lookup = dict(zip(metadata["module_id"], metadata["module_name"]))

out = pd.DataFrame(
    {
        "module_id": both_modules,
        "module_name": [module_name_lookup.get(m, "") for m in both_modules],
        "s_k": both_s_k.astype(int),
        "s_k_metabolites": both_s_k_metabolites.astype(int),
        "block_size": b_k,
    }
)

out = out.sort_values("module_id").reset_index(drop=True)

out.to_csv(OUT_PATH, index=False)

print("\n=== Block Size Summary ===")
print(f"Retained modules: {len(out)}")
print(f"Total latent dimensions: {b_k.sum()}")
print(f"Minimum block size: {b_k.min()}")
print(f"Maximum block size: {b_k.max()}")

assert b_k.min() >= 2
assert b_k.sum() == 128
