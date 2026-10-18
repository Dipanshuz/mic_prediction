import requests
import json
from pathlib import Path

BASE_URL = "https://dbaasp.org/peptides"

def main():
    response = requests.get(
        BASE_URL,
        headers={"Accept": "application/json"},
        timeout=30
    )

    print("Status:", response.status_code)
    response.raise_for_status()

    data = response.json()

    print("Response keys:", list(data.keys()))
    print("Total records:", data.get("totalCount"))
    print("Records in this response:", len(data.get("data", [])))

    output_dir = Path("data/raw")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "peptide_sample.json"

    with output_file.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    print(f"Saved sample to {output_file}")

    if data.get("data"):
        print("\nFirst record:")
        print(json.dumps(data["data"][0], indent=2))

if __name__ == "__main__":
    main()