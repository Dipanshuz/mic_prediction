import json
import time
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Input and output directories
INPUT_DIR = Path("data/raw/peptide_pages")
DETAILS_DIR = Path("data/raw/peptide_details")
OUTPUT_DIR = Path("data/processed")

DETAILS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://dbaasp.org/peptides"

# Retry temporary network/server errors
session = requests.Session()

retry_strategy = Retry(
    total=4,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
)

adapter = HTTPAdapter(max_retries=retry_strategy)
session.mount("https://", adapter)
session.headers.update({"Accept": "application/json"})


def load_peptide_ids():
    """Read unique peptide IDs from all downloaded metadata pages."""
    peptide_ids = set()

    for file_path in sorted(INPUT_DIR.glob("offset_*.json")):
        with file_path.open(encoding="utf-8") as f:
            page = json.load(f)

        for peptide in page.get("data", []):
            peptide_id = peptide.get("id")

            if peptide_id is not None:
                peptide_ids.add(peptide_id)

    return sorted(peptide_ids)


def fetch_details(peptide_id):
    """Fetch and cache one peptide's detailed API response."""
    output_file = DETAILS_DIR / f"{peptide_id}.json"

    # Reuse a previously saved response.
    if output_file.exists():
        try:
            with output_file.open(encoding="utf-8") as f:
                cached = json.load(f)

            if isinstance(cached, dict) and (
                "data" in cached or "targetActivities" in cached
            ):
                return cached
        except (json.JSONDecodeError, OSError):
            pass

    url = f"{BASE_URL}/{peptide_id}"

    response = session.get(url, timeout=45)
    response.raise_for_status()

    result = response.json()

    if not isinstance(result, dict):
        raise ValueError("Unexpected response format")

    # Save the original response without modifying its fields.
    with output_file.open("w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False)

    return result


def get_peptide_object(result):
    """Handle responses that wrap the peptide in a data field."""
    peptide = result.get("data", result)

    if isinstance(peptide, list):
        if not peptide:
            return {}
        peptide = peptide[0]

    return peptide if isinstance(peptide, dict) else {}


def extract_measurements(peptide):
    """Extract all available target activity measurements."""
    rows = []

    for item in peptide.get("targetActivities") or []:
        target = item.get("targetSpecies") or {}
        measure_group = item.get("activityMeasureGroup") or {}
        unit = item.get("unit") or {}
        medium = item.get("medium") or {}

        rows.append({
            "peptide_id": item.get("peptideId", peptide.get("id")),
            "peptide_name": peptide.get("name"),
            "sequence": peptide.get("sequence"),
            "measurement_id": item.get("id"),
            "target_species": target.get("name"),
            "activity_measure": measure_group.get("name"),
            "activity_measure_value": item.get("activityMeasureValue"),
            "activity": item.get("activity"),
            "concentration": item.get("concentration"),
            "unit": unit.get("name"),
            "medium": medium.get("name"),
            "cfu": item.get("cfu"),
            "cfu_group": (item.get("cfuGroup") or {}).get("name"),
            "ph": item.get("ph"),
            "ionic_strength": item.get("ionicStrength"),
            "salt_type": item.get("saltType"),
            "note": item.get("note"),
            "reference": item.get("reference"),
        })

    return rows


def main():
    if not INPUT_DIR.exists():
        raise FileNotFoundError(
            f"Metadata directory not found: {INPUT_DIR}"
        )

    peptide_ids = load_peptide_ids()
   # peptide_ids = load_peptide_ids()[:50]

    if not peptide_ids:
        raise ValueError("No peptide IDs found in metadata files.")

    print("Unique peptide IDs:", len(peptide_ids))
    print("Starting detailed activity extraction...")

    all_rows = []
    failed_ids = []
    processed = 0

    for peptide_id in peptide_ids:
        try:
            result = fetch_details(peptide_id)
            peptide = get_peptide_object(result)

            # A response without a peptide ID is treated as invalid.
            if not peptide or peptide.get("id") is None:
                raise ValueError("Peptide details missing from response")

            all_rows.extend(extract_measurements(peptide))

        except Exception as error:
            failed_ids.append({
                "peptide_id": peptide_id,
                "error": str(error),
            })

        processed += 1

        # Save progress every 100 peptide records.
        if processed % 100 == 0 or processed == len(peptide_ids):
            df = pd.DataFrame(all_rows)
            df.to_csv(
                OUTPUT_DIR / "activity_measurements_partial.csv",
                index=False,
            )

            pd.DataFrame(failed_ids).to_csv(
                OUTPUT_DIR / "failed_peptide_ids.csv",
                index=False,
            )

            print(
                f"Processed {processed}/{len(peptide_ids)} | "
                f"Measurements collected: {len(all_rows)} | "
                f"Failures: {len(failed_ids)}"
            )

        # Be polite to the API; this is a large collection job.
        time.sleep(0.3)

    # Save final output.
    df = pd.DataFrame(all_rows)

    final_file = OUTPUT_DIR / "activity_measurements.csv"
    df.to_csv(final_file, index=False)

    pd.DataFrame(failed_ids).to_csv(
        OUTPUT_DIR / "failed_peptide_ids.csv",
        index=False,
    )

    print("\nExtraction finished.")
    print("Peptides processed:", processed)
    print("Activity measurements:", len(df))

    if not df.empty:
        print(
            "MIC measurements:",
            df["activity_measure"].fillna("").str.upper().eq("MIC").sum(),
        )

    print("Failed peptide IDs:", len(failed_ids))
    print("Final dataset:", final_file)
    print("Failure log:", OUTPUT_DIR / "failed_peptide_ids.csv")


if __name__ == "__main__":
    main()