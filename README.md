# Antimicrobial Peptide (AMP) Prediction Pipeline

This project implements an end-to-end machine learning pipeline for predicting antimicrobial peptide (AMP) activity from peptide sequences. It downloads data from the DBAASP database and trains a gradient boosting model for classification.

## Project Structure

```
amp_prediction/
├── src/
│   ├── data_collection.py       # Downloads data from DBAASP API
│   ├── data_processing.py       # Cleans, transforms, and prepares data
│   ├── training_and_evaluation.py # Trains ML models and evaluates performance
│   └── utils.py                 # Utility functions
├── data/
│   ├── raw/                     # Raw downloaded data
│   └── processed/               # Processed datasets for modeling
├── models/                      # Trained model files
├── notebooks/                   # Jupyter notebooks for exploration
├── .gitignore                   # Specifies files to ignore from version control
└── README.md                    # Project overview
```

## Installation

1.  Clone the repository:
    ```bash
    git clone <repository-url>
    cd amp_prediction
    ```

2.  Create a Python virtual environment:
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate   # On Windows: .venv\Scripts\activate
    ```

3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### 1. Data Collection
Download AMP data from the DBAASP database:
```bash
python src/data_collection.py
```
This will download data to `data/raw/peptide_sample.json`.

### 2. Data Processing
Clean and prepare the dataset for modeling:
```bash
python src/data_processing.py
```
This will create processed datasets in `data/processed/`.

### 3. Model Training and Evaluation
Train machine learning models and evaluate their performance:
```bash
python src/training_and_evaluation.py
```
This will:
-   Train a Gradient Boosting classifier.
-   Save the trained model to `models/gradient_boosting_model.pkl`.
-   Generate evaluation metrics and save them to `models/evaluation_metrics.json`.
-   Create performance plots in `models/`.

## Model Evaluation

The pipeline evaluates models using the following metrics:
-   Accuracy
-   Precision
-   Recall (Sensitivity)
-   F1-Score
-   ROC-AUC Score

The best model's performance is saved to `models/evaluation_metrics.json`.

