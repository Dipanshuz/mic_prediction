# Strain-Specific MIC Prediction of Antimicrobial Peptides

## 1. Project Overview

The goal of this project is to predict the Minimum Inhibitory Concentration (MIC) of antimicrobial peptides against different bacterial species and strains using machine learning and protein language models.

We will collect experimental antimicrobial activity data from the DBAASP database, process and normalize the measurements, and build models that predict peptide activity against specific bacterial targets.

The main focus is strain-aware MIC prediction rather than predicting antimicrobial activity without considering the target.

## 2. What is MIC?

MIC stands for Minimum Inhibitory Concentration.

It represents the minimum concentration of an antimicrobial agent required to inhibit visible microbial growth under specified experimental conditions.

Lower MIC generally indicates stronger antimicrobial activity under comparable experimental conditions.

MIC values can be reported using different units, so unit normalization and data validation are important parts of this project.

## 3. Data Source

We use DBAASP (Database of Antimicrobial Activity and Structure of Peptides).

- Website: https://dbaasp.org/
- Peptide API: https://dbaasp.org/peptides
- Individual peptide details: `https://dbaasp.org/peptides/{peptide_id}`

The peptide listing API supports pagination using parameters such as `limit` and `offset`.

During our initial API test, the database reported 25,542 peptide records. This count represents peptide records, not the total number of antimicrobial activity measurements.

## 4. Project Workflow

The project is divided into the following stages:

1. Retrieve peptide records from DBAASP.
2. Fetch detailed peptide information and activity measurements.
3. Preserve relationships between peptides, bacterial targets, and individual experimental measurements.
4. Identify MIC measurements while retaining other activity data.
5. Clean and validate the collected dataset.
6. Normalize measurement units and handle concentration ranges carefully.
7. Extract peptide features and protein language model embeddings.
8. Train and evaluate MIC prediction models.
9. Compare baseline models with ESM-2-based approaches.

## 5. Current Implementation

The current data collection script is:

`src/data_collection.py`

The initial implementation has tested the DBAASP API and verified pagination with different offsets.

We have also tested extraction of detailed activity measurements for an individual peptide.

The detailed API response can contain fields such as:

- Peptide ID
- Peptide name
- Amino acid sequence
- Measurement ID
- Target species or strain
- Activity measure, such as MIC or MBC
- Activity value
- Concentration
- Measurement unit
- Experimental medium
- CFU information
- pH and ionic strength
- Notes and references

Not every record contains all these fields.

The current sample extraction saves its output to:

`data/processed/activity_sample.csv`

This sample currently covers one peptide record. It is not yet the complete project dataset.

## 6. Project Structure

```text
mic_prediction/
│
├── src/
│   └── data_collection.py
│
├── data/
│   ├── raw/
│   │   ├── peptide_sample.json
│   │   └── peptide_pages/
│   │
│   └── processed/
│       └── activity_sample.csv
│
├── .gitignore
├── requirements.txt
└── README.md
```

Some data files or directories may not exist until the relevant collection scripts have been run.

Large raw datasets should not normally be committed to GitHub.

## 7. Setup Instructions

### Prerequisites

- Python 3
- Git
- Internet connection

### Clone the Repository

```bash
git clone <repository-url>
cd mic_prediction
```

Replace `<repository-url>` with the actual GitHub repository URL.

### Create a Virtual Environment

On macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### Install Dependencies

```bash
python -m pip install --upgrade pip
python -m pip install requests pandas numpy scikit-learn matplotlib
```

## 8. Run the Current Sample Extraction

Run:

```bash
python src/data_collection.py
```

The current sample script retrieves the detailed record for peptide ID 1 and extracts its available target activity measurements.

The resulting CSV is saved to:

`data/processed/activity_sample.csv`

The script prints:

- Peptide name
- Sequence
- Total activity measurements
- Number of MIC measurements
- Number of other measurements

**Note:** The current script is a sample extraction test, not yet the complete database crawler.

## 9. Important Data Handling Rules

### Preserve Raw Data

Keep the original API responses unchanged.

Perform cleaning on separate processed copies so that the original information remains available for verification.

### Preserve Individual Measurements

One peptide can have multiple measurements against different targets, experimental conditions, or references.

Each measurement should remain a separate record.

Do not automatically average measurements from different experiments.

### MIC and Other Activity Measures

MIC is the primary target for our prediction task.

However, other measurements, such as MBC, should be retained in the extracted activity dataset.

A separate MIC-only dataset can be created for model training.

### Concentration and Activity Fields

The API provides both `activity` and `concentration` fields.

Their meanings and consistency must be verified before selecting the final prediction target.

Do not assume they are interchangeable.

Preserve the original values and units.

### Missing Sequences

Some records, including certain multimer peptide entries, may have an empty sequence at the parent-record level.

These records require additional inspection before they can be used for sequence-based modelling.

### Species and Strain

A bacterial species and a specific bacterial strain are not the same thing.

Preserve the original target name and distinguish explicit strain identifiers from species-level records.

## 10. Planned Machine Learning Approaches

We plan to compare the following approaches.

### Baseline Models

- Amino acid composition (AAC)
- Dipeptide composition (DPC)
- Physicochemical peptide descriptors
- Random Forest
- XGBoost, if appropriate

### Protein Language Model

- Frozen ESM-2 embeddings
- Regression model using peptide embeddings
- Strain-aware modelling using target information

The final modelling approach will depend on the quality and quantity of the collected data.

## 11. Evaluation

The models will be evaluated using appropriate regression metrics, including:

- Mean Absolute Error (MAE)
- Root Mean Squared Error (RMSE)
- R-squared (R²)
- Pearson correlation
- Spearman correlation

We also plan to use sequence-identity-controlled evaluation to reduce overly optimistic results caused by similar peptide sequences appearing across training and test sets.

## 12. Current Status and Next Steps

### Completed

- Tested the DBAASP peptide listing API.
- Verified that different pagination offsets return records.
- Inspected an individual peptide's detailed API response.
- Extracted activity measurements for one sample peptide.
- Saved the sample measurements as a CSV file.

### Pending

- Implement reliable collection of all peptide metadata.
- Retrieve detailed activity measurements systematically.
- Validate pagination and handle API errors and retries.
- Preserve raw API responses.
- Build the complete activity dataset.
- Audit MIC fields, units, ranges, and target information.
- Clean and validate the dataset.
- Train and evaluate the prediction models.


## 13. Dataset audit and cleaning pipeline

After collecting peptide metadata and activity measurements, run from the repository root:

```bash
python -m pip install -r requirement.txt
python src/clean_dataset.py
```

The script prefers `data/processed/activity_measurements.csv`; if the full extraction is absent, it uses `activity_sample.csv` and clearly labels its report as sample-only. It joins activity rows to cached peptide metadata to recover sequence information, preserves raw activity/concentration fields, parses exact values and ranges, standardizes molar units to µM, and estimates mass-to-molar conversions only for canonical linear sequences where a sequence-based molecular weight is available.

Outputs:

- `data/processed/activity_measurements_audited.csv`: all activity types with audit columns.
- `data/processed/mic_measurements_clean.csv`: MIC-only rows for downstream analysis.
- `data/processed/reports/unit_counts.csv`: activity-type/unit counts.
- `data/processed/reports/strain_coverage.csv`: MIC measurement and distinct-peptide counts by target label.
- `data/processed/reports/peptide_lengths.csv`: peptide IDs and sequence lengths.
- `data/processed/reports/cleaning_summary.json`: compact quality summary.

Important: a target label may be species-level rather than strain-level. The script preserves the database label and does not infer a strain when the source does not provide one. Ranges and censored values are retained as bounds/operators rather than replaced with a midpoint. Do not train on a single numeric MIC target until the treatment of ranges, censored values, experimental conditions, and mass-unit conversions has been decided and documented.
