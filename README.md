# MIC Prediction Using Antimicrobial Peptides

## 1. Project Overview

The goal of this project is to predict the Minimum Inhibitory Concentration (MIC) of antimicrobial peptides against different bacterial species and strains using machine learning and protein language models.

We collect experimental antimicrobial activity data from the DBAASP database, process and normalize measurements, and prepare the dataset for strain-aware MIC prediction.

The primary objective is to predict peptide activity against specific bacterial targets rather than treating antimicrobial activity as independent of the target organism.

## 2. What Is MIC?

MIC stands for Minimum Inhibitory Concentration. It is the minimum concentration of an antimicrobial agent required to inhibit visible microbial growth under specified experimental conditions.

- Lower MIC generally indicates stronger antimicrobial activity under comparable conditions.
- MIC values may be reported in different units.
- Experimental conditions, measurement ranges, and target organisms must be considered when comparing measurements.

## 3. Data Source

We use the Database of Antimicrobial Activity and Structure of Peptides (DBAASP).

- Website: https://dbaasp.org/
- Peptide listing API: https://dbaasp.org/peptides
- Individual peptide details: `https://dbaasp.org/peptides/{peptide_id}`

The listing API supports pagination using `limit` and `offset` parameters.

During the initial API test, DBAASP reported 25,542 peptide records. This was the count observed during that test, not the number of activity measurements or a verified current database total.

## 4. Project Workflow

The current data pipeline consists of three implemented stages:

1. **Data collection:** Download and cache paginated peptide metadata from DBAASP.
2. **Activity extraction:** Retrieve detailed peptide records and extract individual experimental activity measurements.
3. **Data cleaning and audit:** Recover peptide metadata, parse measurement values, normalize supported units, isolate MIC measurements, and generate data-quality reports.

The planned modelling stages are:

4. Extract peptide sequence features and protein language model embeddings.
5. Train baseline regression models.
6. Develop ESM-2-based prediction approaches.
7. Evaluate performance using appropriate regression metrics and sequence-aware data splits.

## 5. Current Implementation

### 5.1 Data Collection

Script: `src/data_collection.py`

The collection script:

- Retrieves the total record count reported by the API.
- Downloads peptide listing pages in batches of 100 records.
- Saves each page as a separate JSON file.
- Reuses previously downloaded pages when their cached content passes basic checks.
- Validates the API response structure and raises errors for unsuccessful requests.
- Introduces a delay between page requests.

Output directory:

`data/raw/peptide_pages/`

Example files:

- `offset_0.json`
- `offset_100.json`
- `offset_200.json`

### 5.2 Activity Extraction

Script: `src/extract_activities.py`

The extraction script:

- Reads peptide IDs from the downloaded metadata pages.
- Removes duplicate peptide IDs before fetching detailed records.
- Retrieves and caches individual peptide API responses.
- Uses up to five concurrent workers to fetch records.
- Retries specified temporary HTTP errors.
- Extracts available target activity measurements into a tabular dataset.
- Periodically saves progress and records failed peptide requests.

Extracted fields include:

- Peptide ID, name, and sequence
- Measurement ID
- Target species or strain label
- Activity measure and its reported value
- Activity and concentration fields
- Measurement unit
- Experimental medium
- CFU information
- pH, ionic strength, and salt type
- Notes and references

Field availability varies across records.

Output files:

- `data/processed/activity_measurements.csv`
- `data/processed/activity_measurements_partial.csv`
- `data/processed/failed_peptide_ids.csv`
- `data/raw/peptide_details/`

The final CSV contains individual activity measurements, not one averaged row per peptide.

### 5.3 Dataset Cleaning and Audit

Script: `src/clean_dataset.py`

The cleaning script prepares an audited dataset while preserving the original activity and concentration fields.

Its current functionality includes:

- Loading the full activity dataset, or falling back to the sample dataset if the full extraction is unavailable.
- Recovering peptide sequences and metadata from cached JSON records.
- Handling missing parent sequences using available monomer information.
- Calculating peptide sequence lengths and checking sequence alphabets.
- Parsing numeric values, ranges, and inequality operators.
- Classifying supported units into molar and mass-concentration families.
- Converting supported molar measurements to micromolar units (µM).
- Estimating mass-to-molar conversions for suitable canonical linear peptide sequences when molecular weight can be calculated.
- Separating MIC measurements from other activity types.
- Generating summaries of measurement units, target coverage, peptide lengths, and data quality.

The script retains non-MIC records in the audited dataset and preserves ranges and inequalities rather than automatically replacing them with midpoints.

**Important:** the generated MIC dataset still requires validation before it is used for model training. In particular, the appropriate measurement field, treatment of censored values, conversion assumptions, and duplicate records must be reviewed.

## 6. Project Structure

```text
mic_prediction/
│
├── src/
│   ├── data_collection.py
│   ├── extract_activities.py
│   └── clean_dataset.py
│
├── data/
│   ├── raw/
│   │   ├── peptide_pages/
│   │   └── peptide_details/
│   │
│   └── processed/
│       ├── activity_measurements.csv
│       ├── activity_measurements_partial.csv
│       ├── activity_measurements_audited.csv
│       ├── mic_measurements_clean.csv
│       ├── failed_peptide_ids.csv
│       └── reports/
│           ├── unit_counts.csv
│           ├── strain_coverage.csv
│           ├── peptide_lengths.csv
│           └── cleaning_summary.json
│
├── requirements.txt
├── .gitignore
└── README.md
```

This structure describes the expected outputs after running the pipeline. Files may be absent until the relevant scripts have executed successfully. The sample dataset may also be present, depending on the current repository state.

Large raw datasets and generated outputs should generally be excluded from Git commits unless there is a specific reason to version them.

## 7. Setup Instructions

### Prerequisites

- Python 3
- Git
- Internet connection

### Clone the repository

```bash
git clone https://github.com/dhriti046/mic_prediction.git
cd mic_prediction
```

### Create a virtual environment

On macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

On Windows:

```powershell
python -m venv .venv
.venv\Scripts\activate
```

### Install dependencies

From the repository root:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Ensure that `requirements.txt` contains the dependencies required by the scripts, including `requests`, `pandas`, and `numpy`. Biopython is needed for the molecular-weight estimation functionality.

## 8. Running the Data Pipeline

Run the scripts from the repository root in the following order.

### Step 1: Collect peptide metadata

```bash
python src/data_collection.py
```

Downloads and caches paginated peptide listing responses under `data/raw/peptide_pages/`.

### Step 2: Extract activity measurements

```bash
python src/extract_activities.py
```

Fetches individual peptide details and generates `data/processed/activity_measurements.csv`, along with a failure log.

### Step 3: Clean and audit the dataset

```bash
python src/clean_dataset.py
```

Uses `activity_measurements.csv` when available. Otherwise, it falls back to `activity_sample.csv`.

The main output for downstream MIC analysis is:

`data/processed/mic_measurements_clean.csv`

If the sample fallback is used, the cleaning statistics represent the sample rather than the full DBAASP dataset.

## 9. Important Data Handling Rules

### Preserve Raw Data

Original API responses are retained separately from processed datasets. Cleaning and normalization should be performed on derived copies so that measurements can be traced back to their source.

### Preserve Individual Measurements

One peptide may have multiple measurements across target organisms, conditions, and references. Each experimental measurement should remain a separate record.

Measurements from different experiments should not be averaged automatically.

### MIC and Other Activity Measures

MIC is the primary prediction target. Other activity measurements, such as MBC, are retained in the audited activity dataset but excluded from the MIC-only dataset.

### Activity and Concentration Fields

The API provides both `activity` and `concentration` fields. Their meanings must be checked against the source schema and individual records before choosing a final target.

The cleaning script currently parses `concentration` when that column exists, otherwise it falls back to `activity`. This selection is not sufficient by itself to establish that the chosen field is the correct MIC value for every record.

### Units and Measurement Ranges

Original values and units must be preserved. Unsupported units should be investigated rather than silently discarded.

Exact values, ranges, and inequality measurements represent different types of information. A value reported as `< 2` should not be treated as an exact measurement of 2 without an explicit modelling decision.

Mass-to-molar conversions are estimates and depend on the molecular-weight assumptions used.

### Missing Sequences and Multimers

Some peptide records may have incomplete parent-level sequences or multiple monomer sequences. These records require additional inspection before sequence-based modelling.

### Species and Strain

A species-level label is not equivalent to a specific bacterial strain. Preserve the original target label and distinguish explicit strain identifiers from species-level records when the source supports that distinction.

## 10. Current Status and Next Steps

### Implemented in Code

- [x] Tested the DBAASP peptide listing API and pagination.
- [x] Implemented paginated metadata collection and local caching.
- [x] Implemented detailed peptide-record retrieval.
- [x] Implemented extraction of target activity measurements.
- [x] Added concurrent requests, retry handling, and failure logging.
- [x] Implemented metadata recovery and sequence checks.
- [x] Implemented parsing of exact values, ranges, and inequalities.
- [x] Implemented supported unit conversions and conversion-status flags.
- [x] Implemented MIC-only dataset generation.
- [x] Implemented unit, target-coverage, peptide-length, and cleaning-summary reports.

These checkboxes indicate that the corresponding functionality exists in the code. They do not establish that the full database extraction has completed successfully or that all generated records have been validated.

### Remaining Work
- [ ] Review failed requests and missing records.
- [ ] Define how ranges and censored MIC measurements will be handled.
- [ ] Define target-species and strain representation.
- [ ] Finalize the modelling dataset and feature set.
- [ ] Train baseline regression models.
- [ ] Implement ESM-2-based prediction.
- [ ] Evaluate models using sequence-aware data splits.

## 11. Dataset Audit Outputs

After running the cleaning script, the following files are generated under `data/processed/`.

| File | Purpose |
|---|---|
| `activity_measurements_audited.csv` | All activity types with added metadata and audit columns |
| `mic_measurements_clean.csv` | MIC-only records with parsed values and supported unit conversions |
| `reports/unit_counts.csv` | Record counts by normalized activity type and reported unit |
| `reports/strain_coverage.csv` | MIC measurement and distinct-peptide counts by target label |
| `reports/peptide_lengths.csv` | Distinct peptide IDs and sequence lengths |
| `reports/cleaning_summary.json` | Summary of dataset size, MIC coverage, missing sequences, unit conversions, and parsing issues |


