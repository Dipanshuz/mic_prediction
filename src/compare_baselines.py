import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

data_dir = Path("data/processed/splits")
model_dir = Path("models")
report_dir = Path("reports")

sequence_model = joblib.load(model_dir / "ridge_baseline.joblib")

with open(report_dir / "target_only_metrics.json", "r") as f:
    target_metrics = json.load(f)

rows = []

for split_name in ["train", "validation", "test"]:
    df = pd.read_csv(data_dir / f"{split_name}.csv", low_memory=False)
    actual = df["log10_mic_uM"].to_numpy()

    predicted = sequence_model.predict(
        df[["sequence", "target_species"]]
    )

    rows.append({
        "model": "Sequence + target",
        "split": split_name,
        "rows": len(df),
        "MAE_log10": mean_absolute_error(actual, predicted),
        "RMSE_log10": np.sqrt(mean_squared_error(actual, predicted)),
        "R2": r2_score(actual, predicted),
    })

    values = target_metrics[split_name]
    rows.append({
        "model": "Target only",
        "split": split_name,
        "rows": values["rows"],
        "MAE_log10": values["MAE_log10"],
        "RMSE_log10": values["RMSE_log10"],
        "R2": values["R2"],
    })

comparison = pd.DataFrame(rows)
comparison.to_csv(report_dir / "baseline_comparison.csv", index=False)

test_comparison = comparison[comparison["split"] == "test"].copy()
test_comparison = test_comparison.sort_values("MAE_log10")
test_comparison.to_json(
    report_dir / "baseline_test_comparison.json",
    orient="records",
    indent=2,
)

print("BASELINE COMPARISON")
print(comparison.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

print("\nTEST-SET COMPARISON")
print(test_comparison.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

print("\nSaved:")
print("reports/baseline_comparison.csv")
print("reports/baseline_test_comparison.json")
