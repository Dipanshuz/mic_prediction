# DBAASP MIC Dataset Analysis

- Input: `C:\Users\nehal\OneDrive\Desktop\Projects\mic_prediction\data\processed\activity_measurements.csv`
- Extraction complete: **False**
- Activity rows analyzed: **197,338**
- MIC measurements analyzed: **145,745**
- Metadata peptide count: **25542**
- Detail-file coverage: **13.86%**
- Failed peptide IDs: **0**

> If extraction is incomplete, this report describes only the current exported rows.

## MIC concentrations

- Missing concentration values: **19**
- Unique raw concentration/unit pairs: **7,344**
- Exact numeric measurements: **102,517**
- Ranges: **7,656**; censored: **31,509**; unparsed: **4,044**

Exact-value summaries are grouped by compatible canonical unit. Ranges and censored values are excluded from those numeric summaries.

| Unit | Exact n | Min | Q1 | Median | Q3 | Q90 | Q95 | Max | Mean | Missing | Range | Censored |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Unknown unit | 15 | 0.015 | 1.6 | 3.13 | 14.25 | 89.6 | 128 | 128 | 23.27 | 14 | 4 | 4 |
| uM | 51694 | 6.1e-05 | 3.1 | 8 | 25 | 64 | 128 | 3e+04 | 47.97 | 3 | 4452 | 16498 |
| ug/mL | 50808 | 2.2e-05 | 4 | 12.5 | 32 | 125 | 211 | 3.615e+04 | 57.34 | 2 | 3200 | 15007 |

Unit counts and percentages: see `mic_units.csv`. Per-measurement normalized values: `mic_measurements_normalized.csv`; raw concentration value frequencies: `mic_concentration_frequencies.csv`.

## Targets

- Unique non-empty target labels: **6,321**
- Most-tested target label(s) by distinct peptides: **Escherichia coli ATCC 25922** (7,464 peptides each)
- Most MIC measurements: **Escherichia coli ATCC 25922** (9,548 measurements each)
- Peptides tested against multiple target labels: **18,671**

Every target label, MIC measurement count, and unique peptide count is in `mic_targets.csv`. Labels are not automatically parsed into species and strain.

## Peptide lengths

- Unique peptides with MIC records: **20,608**
- Peptides with a length: **20,608**
- Min / max: **1 / 190 aa**
- Mean / median / mode: **17.32 / 14 / 13 aa**

Lengths are counted once per peptide, not once per measurement. Multimer lengths may be sums of component lengths. Per-peptide provenance is in `peptide_lengths.csv`; frequencies are in `peptide_length_frequencies.csv`.

## Data quality

- Exact duplicate rows (additional copies): **0**
- Duplicate measurement-ID groups / additional rows: **0 / 0**
- Repeated measurement-key groups / rows: **1,199 / 2,683**
- Peptide/target pairs with inconsistent unit dimensions: **315**

Molar units are standardized to uM; mass/volume units to ug/mL. No mass-to-molar conversion is attempted because molecular weights and modification states are not validated.

See `repeated_measurements.csv`, `inconsistent_unit_pairs.csv`, and `summary.json` for details.

## Plots

- `mic_concentration_histograms.png`: exact concentrations, separated by canonical unit.
- `peptide_length_histogram.png`: one length per peptide.

## Interpretation notes

- MIC numeric statistics use the raw concentration field because the unit column describes that field. The separate activity field is retained in the source data and is not assumed to be the concentration.
- Range and threshold measurements are counted separately and excluded from exact-value mean, median, quantile, and histogram calculations.
- Target labels are reported exactly as supplied by DBAASP; species-level labels are not treated as confirmed strains.
- For a multimer with no parent sequence, length uses the sum of monomer sequence lengths when available; this is a component-length sum, not a claim that the chains form one linear sequence.
