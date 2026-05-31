import os
import time
import hashlib
import requests
import pandas as pd

KEGG_BASE = "https://rest.kegg.jp"
RAW_DIR = "data/raw/"
OUTPUT_PATH = os.path.join("data/processed", "pathway_metabolite_membership.csv")

REQUEST_DELAY = 0.4

TARGET_PATHWAYS = {
    "map00010": "Glycolysis / Gluconeogenesis",
    "map00020": "Citrate cycle (TCA cycle)",
    "map00260": "Glycine, serine and threonine metabolism",
    "map00670": "One carbon pool by folate",
}


def get_sha256(filepath: str) -> str:
    """Chunked SHA256 file hash."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def fetch_pathway_compound_links() -> list[tuple[str, str]]:
    """GET /link/compound/pathway -> (pathway_id, compound_id)."""
    response = requests.get(f"{KEGG_BASE}/link/compound/pathway")

    if response.status_code != 200:
        raise RuntimeError(f"HTTP {response.status_code}: KEGG fetch failed")

    pairs = []

    for line in response.text.splitlines():
        if "\t" not in line:
            continue

        a, b = line.split("\t")

        pathway_id = None
        compound_id = None

        for x in (a, b):
            if x.startswith("path:"):
                pathway_id = x.replace("path:", "", 1)
            elif x.startswith("cpd:"):
                compound_id = x.replace("cpd:", "", 1)

        if pathway_id and compound_id:
            pairs.append((pathway_id, compound_id))

    time.sleep(REQUEST_DELAY)
    return pairs


def load_ccle_mappings() -> pd.DataFrame:
    """Load CCLE → KEGG compound mappings."""
    df = pd.read_csv(os.path.join("data/processed", "metabolite_kegg_mapping.csv"))

    allowed = {"kegg_exact", "kegg_manual"}
    df = df[df["source"].isin(allowed)].copy()

    df = df[["ccle_name", "kegg_compound_id"]].dropna()

    time.sleep(REQUEST_DELAY)
    return df


def normalize_kegg_ids(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace("cpd:", "", regex=False)


def build_validation_set() -> pd.DataFrame:
    pathway_links = fetch_pathway_compound_links()
    mappings = load_ccle_mappings()

    pathway_df = pd.DataFrame(pathway_links, columns=["pathway_id", "kegg_compound_id"])

    pathway_df["kegg_compound_id"] = normalize_kegg_ids(pathway_df["kegg_compound_id"])
    mappings["kegg_compound_id"] = normalize_kegg_ids(mappings["kegg_compound_id"])

    pathway_df = pathway_df[pathway_df["pathway_id"].isin(TARGET_PATHWAYS)].copy()

    print("\n=== DEBUG PATHWAYS ===")
    print("Rows:", len(pathway_df))

    overlap = set(pathway_df["kegg_compound_id"]) & set(mappings["kegg_compound_id"])
    print("Overlap size:", len(overlap))

    if len(overlap) == 0:
        raise ValueError("No overlap between KEGG pathway compounds and CCLE mappings")

    pathway_df["pathway_name"] = pathway_df["pathway_id"].map(TARGET_PATHWAYS)

    merged = pathway_df.merge(
        mappings,
        on="kegg_compound_id",
        how="inner",
    )

    merged = merged[
        ["pathway_id", "pathway_name", "kegg_compound_id", "ccle_name"]
    ].drop_duplicates()

    return merged


def print_summary(df: pd.DataFrame) -> None:
    print("\n=== PATHWAY METABOLITE COUNTS ===")

    if len(df) == 0:
        print("Empty dataframe.")
        return

    counts = (
        df.groupby(["pathway_id", "pathway_name"])
        .size()
        .reset_index(name="metabolites")
    )

    print(counts.to_string(index=False))
    print("\nTotal rows:", len(df))


def main() -> None:
    df = build_validation_set()

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

    df.to_csv(OUTPUT_PATH, index=False)

    print_summary(df)

    digest = get_sha256(OUTPUT_PATH)

    print("\n=== SHA256 ===")
    print(digest)


if __name__ == "__main__":
    main()
