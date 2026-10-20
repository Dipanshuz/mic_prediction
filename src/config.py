from pathlib import Path

# Project directories
ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "data" / "processed"
SPLIT_DIR = DATA_DIR / "splits"
MODEL_DIR = ROOT_DIR / "models"
REPORT_DIR = ROOT_DIR / "reports"

# Dataset configuration
TARGET_COLUMN = "log10_mic_uM"
SEQUENCE_COLUMN = "sequence"
STRAIN_COLUMN = "target_species"

# Existing reproducible splits
SPLIT_NAMES = ("train", "validation", "test")

# Initial baseline model configuration
RIDGE_ALPHA = 10.0

# Dataset used for the current baseline
MODEL_READY_DATA = DATA_DIR / "mic_model_ready.csv"
BENCHMARK_DATA = DATA_DIR / "mic_benchmark.csv"
