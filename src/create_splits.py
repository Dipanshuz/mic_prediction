import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from pathlib import Path

source = Path("data/processed/mic_benchmark.csv")
output_dir = Path("data/processed/splits")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source, low_memory=False)
df = df.dropna(subset=["sequence", "log10_mic_uM"]).copy()

# First: 80% train, 20% temporary data, grouping by peptide sequence
split1 = GroupShuffleSplit(
    n_splits=1, test_size=0.20, random_state=42
)
train_idx, temp_idx = next(
    split1.split(df, groups=df["sequence"])
)

train = df.iloc[train_idx].copy()
temp = df.iloc[temp_idx].copy()

# Second: divide the remaining 20% equally into validation and test
split2 = GroupShuffleSplit(
    n_splits=1, test_size=0.50, random_state=42
)
val_idx, test_idx = next(
    split2.split(temp, groups=temp["sequence"])
)

validation = temp.iloc[val_idx].copy()
test = temp.iloc[test_idx].copy()

datasets = {
    "train": train,
    "validation": validation,
    "test": test,
}

for name, data in datasets.items():
    data.to_csv(output_dir / f"{name}.csv", index=False)
    print(
        f"{name}: {len(data)} rows, "
        f"{data['sequence'].nunique()} unique sequences, "
        f"{data['target_species'].nunique()} targets"
    )

# Verify there is no sequence overlap
train_sequences = set(train["sequence"])
val_sequences = set(validation["sequence"])
test_sequences = set(test["sequence"])

print("\nSequence overlap checks:")
print("Train vs validation:", len(train_sequences & val_sequences))
print("Train vs test:", len(train_sequences & test_sequences))
print("Validation vs test:", len(val_sequences & test_sequences))
print("\nSaved split files in:", output_dir)
