
import json
import joblib
import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from src.config import (
    SPLIT_DIR,
    MODEL_DIR,
    REPORT_DIR,
    SEQUENCE_COLUMN,
    STRAIN_COLUMN,
    TARGET_COLUMN,
    RIDGE_ALPHA,
)

MODEL_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

# Load the existing splits
train = pd.read_csv(SPLIT_DIR / "train.csv", low_memory=False)
validation = pd.read_csv(SPLIT_DIR / "validation.csv", low_memory=False)
test = pd.read_csv(SPLIT_DIR / "test.csv", low_memory=False)

features = [SEQUENCE_COLUMN, STRAIN_COLUMN]
target = TARGET_COLUMN

# Convert sequence and target labels into machine-readable features
preprocessor = ColumnTransformer(
    transformers=[
        (
            "sequence",
            TfidfVectorizer(
                analyzer="char",
                ngram_range=(1, 3),
                min_df=2,
                max_features=50000,
                sublinear_tf=True,
            ),
            SEQUENCE_COLUMN,
        ),
        (
            "bacterial_target",
            OneHotEncoder(handle_unknown="ignore"),
            [STRAIN_COLUMN],
        ),
    ]
)

model = Pipeline(
    steps=[
        ("features", preprocessor),
        ("regression", Ridge(alpha=RIDGE_ALPHA)),
    ]
)

X_train = train[features]
y_train = train[target]

X_val = validation[features]
y_val = validation[target]

X_test = test[features]
y_test = test[target]

print("Training baseline model...")
model.fit(X_train, y_train)


def evaluate(name, X, y):
    predictions = model.predict(X)

    results = {
        "dataset": name,
        "rows": len(y),
        "MAE": float(mean_absolute_error(y, predictions)),
        "RMSE": float(np.sqrt(mean_squared_error(y, predictions))),
        "R2": float(r2_score(y, predictions)),
    }

    print(
        f"{name}: rows={results['rows']}, "
        f"MAE={results['MAE']:.4f}, "
        f"RMSE={results['RMSE']:.4f}, "
        f"R2={results['R2']:.4f}"
    )

    return results


results = [
    evaluate("train", X_train, y_train),
    evaluate("validation", X_val, y_val),
    evaluate("test", X_test, y_test),
]

# Save the trained model and evaluation results
joblib.dump(model, MODEL_DIR / "ridge_baseline.joblib")

with open(REPORT_DIR / "baseline_metrics.json", "w") as f:
    json.dump(results, f, indent=2)

print("\nModel saved to:", MODEL_DIR / "ridge_baseline.joblib")
print("Metrics saved to:", REPORT_DIR / "baseline_metrics.json")
