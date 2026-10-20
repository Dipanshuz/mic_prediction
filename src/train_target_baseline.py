
import json
import joblib
import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from src.config import (
    SPLIT_DIR,
    MODEL_DIR,
    REPORT_DIR,
    STRAIN_COLUMN,
    TARGET_COLUMN,
    SPLIT_NAMES,
    RIDGE_ALPHA,
)

MODEL_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

features = [STRAIN_COLUMN]
target_column = TARGET_COLUMN

preprocessor = ColumnTransformer([
    (
        "target",
        Pipeline([
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("encoder", OneHotEncoder(handle_unknown="ignore")),
        ]),
        features,
    ),
])

model = Pipeline([
    ("preprocessor", preprocessor),
    ("regressor", Ridge(alpha=RIDGE_ALPHA)),
])

metrics = {}

for split_name in SPLIT_NAMES:
    df = pd.read_csv(SPLIT_DIR / f"{split_name}.csv", low_memory=False)

    X = df[features]
    y = df[target_column].to_numpy()

    if split_name == "train":
        model.fit(X, y)

    predictions = model.predict(X)

    metrics[split_name] = {
        "rows": len(df),
        "MAE_log10": float(mean_absolute_error(y, predictions)),
        "RMSE_log10": float(np.sqrt(mean_squared_error(y, predictions))),
        "R2": float(r2_score(y, predictions)),
    }

joblib.dump(model, MODEL_DIR / "target_only_baseline.joblib")

with open(REPORT_DIR / "target_only_metrics.json", "w") as f:
    json.dump(metrics, f, indent=2)

print("Target-only baseline trained.")
print("\nMetrics:")
print(json.dumps(metrics, indent=2))
print("\nModel saved to:", MODEL_DIR / "target_only_baseline.joblib")
print("Metrics saved to:", REPORT_DIR / "target_only_metrics.json")
