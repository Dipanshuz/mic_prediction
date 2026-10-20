import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

test_path = Path("data/processed/splits/test.csv")
model_path = Path("models/ridge_baseline.joblib")
output_dir = Path("reports")
output_dir.mkdir(parents=True, exist_ok=True)

test = pd.read_csv(test_path, low_memory=False)
model = joblib.load(model_path)

features = ["sequence", "target_species"]
actual_log = test["log10_mic_uM"].to_numpy()
predicted_log = model.predict(test[features])

results = pd.DataFrame({
    "peptide_id": test["peptide_id"],
    "sequence": test["sequence"],
    "target_species": test["target_species"],
    "actual_mic_uM": np.power(10, actual_log),
    "predicted_mic_uM": np.power(10, predicted_log),
    "actual_log10_mic": actual_log,
    "predicted_log10_mic": predicted_log,
    "absolute_log_error": np.abs(actual_log - predicted_log),
})

results.to_csv(output_dir / "baseline_test_predictions.csv", index=False)

metrics = {
    "test_rows": len(test),
    "MAE_log10": float(mean_absolute_error(actual_log, predicted_log)),
    "RMSE_log10": float(np.sqrt(mean_squared_error(actual_log, predicted_log))),
    "R2": float(r2_score(actual_log, predicted_log)),
}

with open(output_dir / "baseline_metrics.json", "w") as f:
    json.dump(metrics, f, indent=2)

print("Test predictions saved to reports/baseline_test_predictions.csv")
print("\nFirst 10 predictions:")
print(results[[
    "target_species",
    "actual_mic_uM",
    "predicted_mic_uM",
    "absolute_log_error",
]].head(10).to_string(index=False))

print("\nMetrics:")
print(json.dumps(metrics, indent=2))
