import os
import time
import hashlib
import requests
import pandas as pd

KEGG_BASE = "https://rest.kegg.jp"
RAW_DIR = "data/raw/"
EXPECTED_RELEASE = "118.0"
REQUEST_DELAY = 0.4


def get_sha256(filepath):
    """SHA256 of a file, chunked read. Copied from 01_verify_data.py."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def verify_release(expected: str = EXPECTED_RELEASE) -> str:
    """GET /info/kegg, parse the release line, assert `expected` in it.
    Returns the full release line for logging. Raises on mismatch."""
    response = requests.get(f"{KEGG_BASE}/info/kegg")

    if response.status_code != 200:
        raise RuntimeError(f"HTTP {response.status_code}: failed to fetch KEGG info")

    lines = response.text.splitlines()
    release_line = next((line for line in lines if "Release" in line), None)

    if release_line is None:
        raise ValueError("No 'Release' line found in KEGG response")

    release_line = release_line.strip()

    if expected not in release_line:
        raise AssertionError(
            f"Expected '{expected}' not found in release line: {release_line}"
        )

    time.sleep(REQUEST_DELAY)
    return release_line


def fetch_gene_module_links() -> list[tuple[str, str]]:
    """GET /link/hsa/module. Returns list of (module_id, hsa_id),
    prefixes stripped (md:, hsa:)."""
    response = requests.get(f"{KEGG_BASE}/link/hsa/module")

    if response.status_code != 200:
        raise RuntimeError(f"HTTP {response.status_code}: failed gene-module fetch")

    pairs = []

    for line in response.text.splitlines():
        a, b = line.split("\t")

        module_id = None
        hsa_id = None

        for p in (a, b):
            if p.startswith("md:"):
                module_id = p.replace("md:", "", 1)
                if module_id.startswith("hsa_"):
                    module_id = module_id.replace("hsa_", "", 1)
            elif p.startswith("hsa:"):
                hsa_id = p.replace("hsa:", "", 1)

        if module_id and hsa_id:
            pairs.append((module_id, hsa_id))

    time.sleep(REQUEST_DELAY)
    return pairs


def fetch_hsa_symbol_map() -> dict[str, str]:
    """GET /list/hsa. Returns {hsa_id: hgnc_symbol},
    symbol = token before first ';' in column 2."""
    response = requests.get(f"{KEGG_BASE}/list/hsa")

    if response.status_code != 200:
        raise RuntimeError(
            f"HTTP {response.status_code}: failed to fetch KEGG hsa list"
        )

    hsa_to_symbol = {}
    miss_count = 0

    for line in response.text.splitlines():
        fields = line.split("\t")

        if len(fields) < 4:
            miss_count += 1
            continue

        hsa_id = fields[0].replace("hsa:", "", 1)
        annotation = fields[3]

        if ";" not in annotation:
            miss_count += 1
            continue

        token = annotation.split(";", 1)[0].split(",", 1)[0].strip()

        if token:
            hsa_to_symbol[hsa_id] = token
        else:
            miss_count += 1

    print(f"Symbol-map misses: {miss_count}")
    time.sleep(REQUEST_DELAY)
    return hsa_to_symbol


def fetch_compound_module_links() -> list[tuple[str, str]]:
    """GET /link/compound/module. Returns list of (module_id, compound_id),
    prefixes stripped (md:, cpd:), column order resolved by prefix check."""
    response = requests.get(f"{KEGG_BASE}/link/compound/module")

    if response.status_code != 200:
        raise RuntimeError(f"HTTP {response.status_code}: failed compound-module fetch")

    pairs = []

    for line in response.text.splitlines():
        a, b = line.split("\t")

        module_id = None
        cpd_id = None

        for p in (a, b):
            if p.startswith("md:"):
                module_id = p.replace("md:", "", 1)
            elif p.startswith("cpd:"):
                cpd_id = p.replace("cpd:", "", 1)

        if module_id and cpd_id:
            pairs.append((module_id, cpd_id))

    time.sleep(REQUEST_DELAY)
    return pairs


def fetch_module_metadata() -> dict[str, str]:
    """GET /list/module. Returns {module_id: module_name}."""
    response = requests.get(f"{KEGG_BASE}/list/module")

    if response.status_code != 200:
        raise RuntimeError(
            f"HTTP {response.status_code}: failed to fetch KEGG module metadata"
        )

    modules = {}

    for line in response.text.splitlines():
        fields = line.split("\t")

        if len(fields) < 2:
            continue

        module_id = fields[0].replace("md:", "", 1).strip()
        module_name = fields[1].strip()

        if module_id and module_name:
            modules[module_id] = module_name

    time.sleep(REQUEST_DELAY)
    return modules


def main() -> None:
    """Orchestrate: verify -> fetch all -> join gene links to symbols ->
    serialize 3 TSVs -> SHA256 -> print summary."""
    release = verify_release()
    gene_links = fetch_gene_module_links()
    symbol_map = fetch_hsa_symbol_map()
    compound_links = fetch_compound_module_links()
    module_data = fetch_module_metadata()

    joined_gene = []
    join_miss = 0

    for module_id, hsa_id in gene_links:
        symbol = symbol_map.get(hsa_id)

        if symbol is None:
            join_miss += 1
            continue

        joined_gene.append((module_id, symbol))

    print("Join misses:", join_miss)

    os.makedirs(RAW_DIR, exist_ok=True)

    gene_df = pd.DataFrame(joined_gene, columns=["module_id", "hgnc_symbol"])
    compound_df = pd.DataFrame(
        compound_links, columns=["module_id", "kegg_compound_id"]
    )
    metadata_df = pd.DataFrame(
        list(module_data.items()), columns=["module_id", "module_name"]
    )

    print("Gene rows before:", len(gene_df))
    gene_df = gene_df.drop_duplicates()
    print("Gene rows after:", len(gene_df))

    print("Compound rows before:", len(compound_df))
    compound_df = compound_df.drop_duplicates()
    print("Compound rows after:", len(compound_df))

    gene_path = os.path.join(RAW_DIR, "kegg_module_gene_membership.tsv")
    compound_path = os.path.join(RAW_DIR, "kegg_module_compound_membership.tsv")
    metadata_path = os.path.join(RAW_DIR, "kegg_module_metadata.tsv")

    gene_df.to_csv(gene_path, sep="\t", index=False)
    compound_df.to_csv(compound_path, sep="\t", index=False)
    metadata_df.to_csv(metadata_path, sep="\t", index=False)

    gene_hash = get_sha256(gene_path)
    compound_hash = get_sha256(compound_path)
    metadata_hash = get_sha256(metadata_path)

    print("Gene SHA256:", gene_hash)
    print("Compound SHA256:", compound_hash)
    print("Metadata SHA256:", metadata_hash)

    print("\n=== SUMMARY ===")
    print("Release:", release.strip())

    print("Gene links:", len(gene_df))
    print("Unique modules (gene):", gene_df["module_id"].nunique())
    print("Unique HGNC symbols:", gene_df["hgnc_symbol"].nunique())

    print("Compound links:", len(compound_df))
    print("Unique modules (compound):", compound_df["module_id"].nunique())
    print("Unique compounds:", compound_df["kegg_compound_id"].nunique())

    print("Symbol-map size:", len(symbol_map))
    print("Join misses:", join_miss)

    print("Modules in metadata:", len(module_data))


if __name__ == "__main__":
    main()
