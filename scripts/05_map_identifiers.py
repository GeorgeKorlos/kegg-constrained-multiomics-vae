import os
import time
import json
import hashlib
import requests
import urllib.parse
import pandas as pd

KEGG_BASE = "https://rest.kegg.jp"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROCESSED_DIR = "data/processed"
METABOLITE_NAMES_PATH = "data/processed/metabolite_names.csv"
REVIEW_OUTPUT_PATH = "data/processed/metabolite_mapping_review.tsv"
KEGG_CACHE_PATH = "data/processed/.kegg_compound_cache.json"
PUBCHEM_CACHE_PATH = "data/processed/.pubchem_cid_cache.json"
MAPPING_OUTPUT_PATH = "data/processed/metabolite_kegg_mapping.csv"
REQUEST_DELAY = 0.4


def get_sha256(filepath):
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def normalize_name(raw_name: str) -> str:
    s = raw_name.strip()
    s = s.lower()
    s = s.replace("’", "'").replace("–", "-").replace("—", "-")
    s = " ".join(s.split())
    return s


def query_kegg_compound(
    normalized_name: str, cache: dict
) -> list[tuple[str, list[str]]]:
    if normalized_name in cache:
        return cache[normalized_name]

    encoded_name = urllib.parse.quote(normalized_name, safe="")
    url = f"{KEGG_BASE}/find/compound/{encoded_name}"

    response = requests.get(url)

    if response.status_code in (400, 404):
        results = []

    elif response.status_code == 200:

        results = []

        if response.text.strip():

            for line in response.text.splitlines():
                line = line.strip()

                if not line:
                    continue

                compound_col, names_col = line.split("\t")

                compound_id = compound_col.replace("cpd:", "", 1)

                names_list = [
                    name.strip() for name in names_col.split(";") if name.strip()
                ]

                results.append((compound_id, names_list))

    else:
        raise RuntimeError(
            f"HTTP {response.status_code}: KEGG/find/compound failed for {normalized_name}"
        )

    cache[normalized_name] = results
    save_cache(cache, KEGG_CACHE_PATH)

    time.sleep(REQUEST_DELAY)
    return results


def classify_kegg_result(
    ccle_name: str, normalized: str, kegg_results: list[tuple[str, list[str]]]
) -> tuple[str, str | list | None]:
    if not kegg_results:
        return ("unmapped", None)

    survivors = [
        (cid, names)
        for cid, names in kegg_results
        if any(syn.strip().lower() == normalized for syn in names)
    ]
    if len(survivors) == 1:
        return ("kegg_exact", survivors[0][0])

    return ("review", kegg_results)


def query_pubchem_cid(normalized_name, cache) -> str | None:
    if normalized_name in cache:
        return cache[normalized_name]

    encoded_name = urllib.parse.quote(normalized_name, safe="")
    url = f"{PUBCHEM_BASE}/compound/name/{encoded_name}/cids/JSON"

    response = requests.get(url)
    cid = None

    if response.status_code == 200:
        body = response.json()
        cid_list = body.get("IdentifierList", {}).get("CID", [])
        cid = str(cid_list[0]) if cid_list else None
    elif response.status_code == 404:
        cid = None

    else:
        raise RuntimeError(
            f"HTTP {response.status_code}: PubChem /compound/name failed for {normalized_name}"
        )
    cache[normalized_name] = cid
    save_cache(cache, PUBCHEM_CACHE_PATH)

    time.sleep(REQUEST_DELAY)
    return cid


def load_cache(path: str) -> dict:
    if not os.path.isfile(path):
        return {}
    else:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data


def save_cache(cache, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(cache, f)


def map_metabolites(metabolite_names_path: str = METABOLITE_NAMES_PATH) -> pd.DataFrame:

    df = pd.read_csv(metabolite_names_path)
    assert len(df) == 225

    kegg_cache = load_cache(KEGG_CACHE_PATH)
    pubchem_cache = load_cache(PUBCHEM_CACHE_PATH)

    rows = []
    review_queue = []

    kegg_exact_count = 0
    review_count = 0
    pubchem_only_count = 0
    unmapped_count = 0

    for i, ccle_name in enumerate(df.iloc[:, 0], start=1):

        normalized = normalize_name(ccle_name)
        kegg_results = query_kegg_compound(normalized, kegg_cache)
        verdict, payload = classify_kegg_result(ccle_name, normalized, kegg_results)

        row = {
            "ccle_name": ccle_name,
            "kegg_compound_id": "",
            "pubchem_cid": "",
            "chebi_id": "",
            "source": "",
            "notes": "",
        }

        if verdict == "kegg_exact":
            row["kegg_compound_id"] = payload
            row["source"] = "kegg_exact"
            kegg_exact_count += 1

        elif verdict == "review":
            row["source"] = "review"
            review_count += 1

            candidate_ids = [cid for cid, _ in kegg_results]
            candidate_names = [names for _, names in kegg_results]

            review_queue.append((ccle_name, normalized, candidate_ids, candidate_names))

        elif verdict == "unmapped":
            cid = query_pubchem_cid(normalized, pubchem_cache)

            if cid:
                row["pubchem_cid"] = cid
                row["source"] = "pubchem_only"
                pubchem_only_count += 1
            else:
                row["source"] = "unmapped"
                unmapped_count += 1

        rows.append(row)

        if i % 25 == 0:
            print(f"  {i}/225 ...")

    mapping_df = pd.DataFrame(
        rows,
        columns=[
            "ccle_name",
            "kegg_compound_id",
            "pubchem_cid",
            "chebi_id",
            "source",
            "notes",
        ],
    )

    mapping_df.to_csv(MAPPING_OUTPUT_PATH, index=False)

    if review_queue:
        review_df = pd.DataFrame(
            [
                {
                    "ccle_name": r[0],
                    "normalized_name": r[1],
                    "candidate_kegg_ids": "|".join(r[2]),
                    "candidate_kegg_names": "|".join(
                        [";".join(names) for names in r[3]]
                    ),
                }
                for r in review_queue
            ]
        )

        review_df.to_csv(REVIEW_OUTPUT_PATH, sep="\t", index=False)

    kegg_coverage = kegg_exact_count / 225 * 100
    cross_ref_coverage = (kegg_exact_count + pubchem_only_count) / 225 * 100

    print("\n=== METABOLITE MAPPING SUMMARY ===")
    print(f"Total CCLE metabolites: 225")
    print(f"Auto-mapped (kegg_exact): {kegg_exact_count}")
    print(f"Manual review queue:     {review_count}")
    print(f"PubChem-only cross-ref:  {pubchem_only_count}")
    print(f"Fully unmapped:          {unmapped_count}")
    print()
    print(f"KEGG coverage:           {kegg_exact_count}/225 ({kegg_coverage:.1f}%)")
    print(
        f"Cross-ref coverage:      {kegg_exact_count + pubchem_only_count}/225 ({cross_ref_coverage:.1f}%)"
    )

    return mapping_df


def main():
    df = map_metabolites()


if __name__ == "__main__":
    main()
