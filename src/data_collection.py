import json
import time
from pathlib import Path

import requests

BASE_URL = "https://dbaasp.org/peptides"
PAGE_SIZE = 100
OUTPUT_DIR = Path("data/raw/peptide_pages")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"Accept": "application/json"})


def fetch_page(offset):
    params = {
        "limit": PAGE_SIZE,
        "offset": offset
    }

    response = session.get(BASE_URL, params=params, timeout=60)
    response.raise_for_status()

    result = response.json()

    if "data" not in result or "totalCount" not in result:
        raise ValueError("Unexpected API response structure")

    return result


def main():
    first_page = fetch_page(0)
    total = first_page["totalCount"]

    print("Total records reported by API:", total)

    downloaded = 0

    for offset in range(0, total, PAGE_SIZE):
        output_file = OUTPUT_DIR / f"offset_{offset}.json"

        # Reuse previously downloaded pages.
        if output_file.exists():
            with output_file.open(encoding="utf-8") as f:
                cached = json.load(f)

            if cached.get("offset") == offset and cached.get("data"):
                downloaded += len(cached["data"])
                continue

        result = first_page if offset == 0 else fetch_page(offset)
        records = result["data"]

        if not records:
            print("No records returned at offset:", offset)
            break

        result["offset"] = offset

        with output_file.open("w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False)

        downloaded += len(records)

        if downloaded % 1000 == 0 or offset == 0:
            print(f"Processed approximately {downloaded}/{total} records")

        time.sleep(0.25)

    print("\nFinished.")
    print("Pages saved in:", OUTPUT_DIR)
    print("Records downloaded or reused:", downloaded)


if __name__ == "__main__":
    main()