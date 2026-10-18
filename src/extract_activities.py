import json
import threading
import time
from concurrent.futures import ThreadPoolExecutor
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

# Number of parallel requests. The API is slow (~1.7 s per request),
# so a few concurrent workers help a lot. Keep this small to stay polite.
MAX_WORKERS = 5

# Delay after each real network request, per worker.
REQUEST_DELAY = 0.3

# How often to report progress and save the partial CSV.
SAVE_EVERY = 500

# Retry temporary network/server errors
retry_strategy = Retry(
    total=4,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
)

# One session per thread, since sessions are not guaranteed thread-safe.
thread_local = threading.local()


def get_session():
    session = getattr(thread_local, "session", None)

    if session is None:
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retry_strategy))
        session.headers.update({"Accept": "application/json"})
        thread_local.session = session

    return session


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

    try:
        response = get_session().get(url, timeout=45)
        response.raise_for_status()

        result = response.json()

        if not isinstance(result, dict):
            raise ValueError("Unexpected response format")

        # Save the original response without modifying its fields.
        # Write to a temp file first so an interrupted run never
        # leaves a half-written cache file behind.
        tmp_file = output_file.with_suffix(".json.tmp")
        with tmp_file.open("w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False)
        tmp_file.replace(output_file)

        return result
    finally:
        # Be polite to the API; only delay after real requests.
        time.sleep(REQUEST_DELAY)


def process_peptide(peptide_id):
    """Fetch one peptide and return (rows, error). Runs in a worker thread."""
    try:
        result = fetch_details(peptide_id)
        peptide = get_peptide_object(result)

        # A response without a peptide ID is treated as invalid.
        if not peptide or peptide.get("id") is None:
            raise ValueError("Peptide details missing from response")

        return extract_measurements(peptide), None

    except Exception as error:
        return [], {"peptide_id": peptide_id, "error": str(error)}


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

    partial_file = OUTPUT_DIR / "activity_measurements_partial.csv"
    partial_file.unlink(missing_ok=True)
    unsaved_rows = []
    start_time = time.time()

    # executor.map returns results in input order, so the output
    # stays sorted by peptide ID even though fetching is parallel.
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for rows, error in executor.map(process_peptide, peptide_ids):
            all_rows.extend(rows)
            unsaved_rows.extend(rows)

            if error:
                failed_ids.append(error)

            processed += 1

            # Save progress periodically by appending only the new rows,
            # instead of rewriting the whole file each time.
            if processed % SAVE_EVERY == 0 or processed == len(peptide_ids):
                if unsaved_rows:
                    pd.DataFrame(unsaved_rows).to_csv(
                        partial_file,
                        mode="a",
                        header=not partial_file.exists(),
                        index=False,
                    )
                    unsaved_rows = []

                pd.DataFrame(failed_ids).to_csv(
                    OUTPUT_DIR / "failed_peptide_ids.csv",
                    index=False,
                )

                elapsed = time.time() - start_time
                print(
                    f"Processed {processed}/{len(peptide_ids)} | "
                    f"Measurements collected: {len(all_rows)} | "
                    f"Failures: {len(failed_ids)} | "
                    f"Elapsed: {elapsed / 60:.1f} min"
                )

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