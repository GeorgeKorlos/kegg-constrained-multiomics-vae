import hashlib
import numpy as np
import pandas as pd

KEGG_GENE_PATH = "data/raw/kegg_module_gene_membership.tsv"
KEGG_COMPOUND_PATH = "data/raw/kegg_module_compound_membership.tsv"
KEGG_MODULE_META_PATH = "data/raw/kegg_module_metadata.tsv"
CCLE_MAPPING_PATH = "data/processed/metabolite_kegg_mapping.csv"
CCLE_GENES_PATH = "data/processed/gene_names.csv"
CCLE_METABOLITES_PATH = "data/processed/metabolite_names.csv"
GENE_MATRIX_OUT = "data/processed/gene_module_matrix.npy"
META_MATRIX_OUT = "data/processed/metabolite_module_matrix.npy"
MODULE_IDS_OUT = "data/processed/module_ids.csv"


def get_sha256(filepath):
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def load_inputs() -> dict:
    kegg_genes = pd.read_csv(KEGG_GENE_PATH, sep="\t")
    assert list(kegg_genes.columns) == ["module_id", "hgnc_symbol"]

    kegg_compounds = pd.read_csv(KEGG_COMPOUND_PATH, sep="\t")
    assert list(kegg_compounds.columns) == ["module_id", "kegg_compound_id"]

    kegg_module_meta = pd.read_csv(KEGG_MODULE_META_PATH, sep="\t")
    assert list(kegg_module_meta.columns) == ["module_id", "module_name"]

    ccle_mapping = pd.read_csv(CCLE_MAPPING_PATH)
    assert ccle_mapping.shape == (225, 6)

    ccle_genes = pd.read_csv(CCLE_GENES_PATH).iloc[:, 0]
    assert len(ccle_genes) == 17384, f"Expected 17384 CCLE genes, got {len(ccle_genes)}"

    ccle_metabolites = pd.read_csv(CCLE_METABOLITES_PATH).iloc[:, 0]
    assert (
        len(ccle_metabolites) == 225
    ), f"Expected 225 CCLE metabolites, got {len(ccle_metabolites)}"

    print(f"kegg_genes:        {kegg_genes.shape}")
    print(f"kegg_compounds:    {kegg_compounds.shape}")
    print(f"kegg_module_meta:  {kegg_module_meta.shape}")
    print(f"ccle_mapping:      {ccle_mapping.shape}")
    print(f"ccle_genes:        {len(ccle_genes)}")
    print(f"ccle_metabolites:  {len(ccle_metabolites)}")

    return {
        "kegg_genes": kegg_genes,
        "kegg_compounds": kegg_compounds,
        "kegg_module_meta": kegg_module_meta,
        "ccle_mapping": ccle_mapping,
        "ccle_genes": ccle_genes,
        "ccle_metabolites": ccle_metabolites,
    }


def filter_gene_module_to_ccle(
    kegg_gene_df: pd.DataFrame,
    ccle_genes: pd.Series,
) -> pd.DataFrame:
    mask = kegg_gene_df["hgnc_symbol"].isin(ccle_genes)
    filtered = kegg_gene_df[mask].copy()

    print(f"Gene-module pre-filter rows:  {len(kegg_gene_df)}")
    print(f"Gene-module post-filter rows: {len(filtered)}")
    print(f"Unique modules: {filtered['module_id'].nunique()}")
    print(f"Unique CCLE genes covered: {filtered['hgnc_symbol'].nunique()}")

    ccle_gene_set = set(ccle_genes)
    annotated = set(filtered["hgnc_symbol"])
    print(f"CCLE genes with ≥1 module:        {len(annotated)}")
    print(
        f"CCLE genes with 0 modules:        {len(ccle_gene_set - annotated)} "
        f"({(len(ccle_gene_set - annotated) / len(ccle_gene_set)) * 100:.1f}%)"
    )
    return filtered


def filter_compound_module_to_ccle(
    kegg_compound_df: pd.DataFrame,
    ccle_mapping_df: pd.DataFrame,
) -> pd.DataFrame:
    mapped = ccle_mapping_df[
        ccle_mapping_df["source"].isin(["kegg_exact", "kegg_manual"])
    ].copy()
    mapped = mapped[["ccle_name", "kegg_compound_id"]]
    assert (
        len(mapped) == 104
    ), f"Expected 104 KEGG-mapped metabolites, got {len(mapped)}"

    filtered = kegg_compound_df.merge(
        mapped,
        on="kegg_compound_id",
        how="inner",
    )

    n_before = len(filtered)
    filtered = filtered.drop_duplicates(subset=["module_id", "ccle_name"])
    n_after = len(filtered)
    if n_before != n_after:
        print(
            f"  warning: dropped {n_before - n_after} duplicate (module, metabolite) pairs"
        )

    print(f"\nCCLE mapping (kegg_exact + kegg_manual):   {len(mapped)}")
    print(f"Compound-module pre-filter rows:           {len(kegg_compound_df)}")
    print(f"Compound-module post-filter rows:          {len(filtered)}")
    print(f"Unique modules: {filtered['module_id'].nunique()}")
    print(
        f"Unique KEGG compounds in filtered table:   {filtered['kegg_compound_id'].nunique()}"
    )
    print(
        f"Unique CCLE metabolites in ≥1 module:      {filtered['ccle_name'].nunique()}"
    )

    mapped_in_modules = set(filtered["ccle_name"])
    mapped_total = set(mapped["ccle_name"])
    mapped_no_module = mapped_total - mapped_in_modules
    print(f"CCLE metabolites with KEGG ID but 0 modules: {len(mapped_no_module)}")

    if mapped_no_module:
        print(f"  examples: {sorted(mapped_no_module)[:5]}")

    return filtered


def classify_modules(
    gene_modules: set[str],
    compound_modules: set[str],
) -> tuple[set[str], set[str], set[str]]:
    both = gene_modules & compound_modules
    gene_only = gene_modules - compound_modules
    metabolite_only = compound_modules - gene_modules
    return both, gene_only, metabolite_only


def determine_module_universe(
    filtered_gene_df: pd.DataFrame,
    filtered_compound_df: pd.DataFrame,
) -> tuple[list[str], set[str], set[str], set[str]]:
    gene_modules = set(filtered_gene_df["module_id"])
    compound_modules = set(filtered_compound_df["module_id"])

    retained = sorted(gene_modules | compound_modules)
    both, gene_only, metabolite_only = classify_modules(gene_modules, compound_modules)

    K = len(retained)
    assert len(both) + len(gene_only) + len(metabolite_only) == K

    print(f"K: {K}")
    print(f"Both modalities: {len(both)}")
    print(f"Gene-only: {len(gene_only)}")
    print(f"Metabolite-only: {len(metabolite_only)}")

    return retained, both, gene_only, metabolite_only


def build_binary_matrix(
    filtered_df: pd.DataFrame,
    row_ids: list[str],
    col_ids: list[str],
    row_col: str,
    col_col: str,
) -> np.ndarray:
    row_lookup = {row_id: i for i, row_id in enumerate(row_ids)}
    col_lookup = {col_id: j for j, col_id in enumerate(col_ids)}
    matrix = np.zeros((len(row_ids), len(col_ids)), dtype=np.float32)

    for row in filtered_df.itertuples(index=False):
        row_key = getattr(row, row_col)
        col_key = getattr(row, col_col)

        if row_key not in row_lookup or col_key not in col_lookup:
            continue

        i = row_lookup[row_key]
        j = col_lookup[col_key]
        matrix[i, j] = 1.0

    return matrix


def verify_matrix(
    matrix: np.ndarray,
    expected_rows: int,
    expected_link_count: int,
    name: str,
) -> None:
    assert matrix.dtype == np.float32
    assert matrix.shape[0] == expected_rows
    assert matrix.sum() == expected_link_count

    row_sums = matrix.sum(axis=1)
    print(
        f"{name} row sums — min: {row_sums.min()}, "
        f"median: {np.median(row_sums)}, max: {row_sums.max()}, "
        f"zero-rows: {(row_sums == 0).sum()}"
    )


def save_outputs(
    gene_matrix: np.ndarray,
    meta_matrix: np.ndarray,
    module_ids: list[str],
) -> None:
    np.save(GENE_MATRIX_OUT, gene_matrix)
    np.save(META_MATRIX_OUT, meta_matrix)
    pd.Series(module_ids, name="module_id").to_csv(MODULE_IDS_OUT, index=False)

    print(f"gene_module_matrix.npy        SHA256: {get_sha256(GENE_MATRIX_OUT)}")
    print(f"metabolite_module_matrix.npy  SHA256: {get_sha256(META_MATRIX_OUT)}")
    print(f"module_ids.csv                SHA256: {get_sha256(MODULE_IDS_OUT)}")


def print_coverage_summary(
    ccle_genes: pd.Series,
    ccle_metabolites: pd.Series,
    ccle_mapping_df: pd.DataFrame,
    kegg_module_meta: pd.DataFrame,
    filtered_gene_df: pd.DataFrame,
    filtered_compound_df: pd.DataFrame,
    retained: list[str],
    both: set[str],
    gene_only: set[str],
    metabolite_only: set[str],
    gene_matrix: np.ndarray,
    meta_matrix: np.ndarray,
) -> None:
    print("\n==============================")
    print("=== K MODULE FILTERING ===")
    print("==============================")

    print(f"Modules in KEGG total: {len(kegg_module_meta)}")
    print(f"Modules with ≥1 CCLE gene: {filtered_gene_df['module_id'].nunique()}")
    print(
        f"Modules with ≥1 CCLE metabolite: {filtered_compound_df['module_id'].nunique()}"
    )
    print(f"Modules retained (union): K = {len(retained)}")
    print(f"  Both modalities:    {len(both)}")
    print(f"  Gene-only:          {len(gene_only)}")
    print(f"  Metabolite-only:    {len(metabolite_only)}")

    print("\n==============================")
    print("=== GENE COVERAGE ===")
    print("==============================")

    gene_row_sums = gene_matrix.sum(axis=1)
    annotated_gene_mask = gene_row_sums > 0
    annotated_gene_sums = gene_row_sums[annotated_gene_mask]

    print(f"CCLE genes: {len(ccle_genes)}")
    print(f"Genes with ≥1 module: {annotated_gene_mask.sum()}")
    print(
        f"Genes with 0 modules (unconstrained): "
        f"{len(ccle_genes) - annotated_gene_mask.sum()} "
        f"({(1 - annotated_gene_mask.mean()) * 100:.1f}%)"
    )
    print(f"Gene-module link count: {int(gene_matrix.sum())}")

    if len(annotated_gene_sums) > 0:
        print(f"Median modules per annotated gene: {np.median(annotated_gene_sums)}")
        print(f"Max modules per gene: {annotated_gene_sums.max()}")
    else:
        print("Median modules per annotated gene: 0")
        print("Max modules per gene: 0")

    print("\n==============================")
    print("=== METABOLITE COVERAGE ===")
    print("==============================")

    meta_row_sums = meta_matrix.sum(axis=1)
    annotated_meta_mask = meta_row_sums > 0
    annotated_meta_sums = meta_row_sums[annotated_meta_mask]

    kegg_mapped = ccle_mapping_df[
        ccle_mapping_df["source"].isin(["kegg_exact", "kegg_manual"])
    ]["ccle_name"]

    print(f"CCLE metabolites: {len(ccle_metabolites)}")
    print(f"Metabolites mapped to KEGG: {len(kegg_mapped)}")
    print(f"Metabolites with ≥1 module: {annotated_meta_mask.sum()}")
    print(
        f"Metabolites with 0 modules (unconstrained): "
        f"{len(ccle_metabolites) - annotated_meta_mask.sum()} "
        f"({(1 - annotated_meta_mask.mean()) * 100:.1f}%)"
    )
    print(f"Metabolite-module link count: {int(meta_matrix.sum())}")

    if len(annotated_meta_sums) > 0:
        print(
            f"Median modules per annotated metabolite: {np.median(annotated_meta_sums)}"
        )
        print(f"Max modules per metabolite: {annotated_meta_sums.max()}")
    else:
        print("Median modules per annotated metabolite: 0")
        print("Max modules per metabolite: 0")


def main() -> None:
    inputs = load_inputs()

    filtered_genes = filter_gene_module_to_ccle(
        inputs["kegg_genes"],
        inputs["ccle_genes"],
    )
    filtered_compounds = filter_compound_module_to_ccle(
        inputs["kegg_compounds"],
        inputs["ccle_mapping"],
    )

    retained, both, gene_only, metabolite_only = determine_module_universe(
        filtered_genes,
        filtered_compounds,
    )

    gene_matrix = build_binary_matrix(
        filtered_genes,
        inputs["ccle_genes"].tolist(),
        retained,
        row_col="hgnc_symbol",
        col_col="module_id",
    )
    meta_matrix = build_binary_matrix(
        filtered_compounds,
        inputs["ccle_metabolites"].tolist(),
        retained,
        row_col="ccle_name",
        col_col="module_id",
    )

    verify_matrix(
        gene_matrix,
        expected_rows=17384,
        expected_link_count=len(filtered_genes),
        name="gene_matrix",
    )
    verify_matrix(
        meta_matrix,
        expected_rows=225,
        expected_link_count=len(filtered_compounds),
        name="meta_matrix",
    )

    combined_col_sums = gene_matrix.sum(axis=0) + meta_matrix.sum(axis=0)
    dead_modules = (combined_col_sums == 0).sum()
    assert (
        dead_modules == 0
    ), f"{dead_modules} retained modules have zero annotations across both matrices"

    print_coverage_summary(
        ccle_genes=inputs["ccle_genes"],
        ccle_metabolites=inputs["ccle_metabolites"],
        ccle_mapping_df=inputs["ccle_mapping"],
        kegg_module_meta=inputs["kegg_module_meta"],
        filtered_gene_df=filtered_genes,
        filtered_compound_df=filtered_compounds,
        retained=retained,
        both=both,
        gene_only=gene_only,
        metabolite_only=metabolite_only,
        gene_matrix=gene_matrix,
        meta_matrix=meta_matrix,
    )

    save_outputs(gene_matrix, meta_matrix, retained)


if __name__ == "__main__":
    main()
