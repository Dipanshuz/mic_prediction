import argparse
import json
import math
import re
import unicodedata
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FINAL_INPUT = ROOT / "data" / "processed" / "activity_measurements.csv"
DEFAULT_PARTIAL_INPUT = ROOT / "data" / "processed" / "activity_measurements_partial.csv"
DEFAULT_OUTPUT_DIR = ROOT / "data" / "processed" / "analysis"
METADATA_DIR = ROOT / "data" / "raw" / "peptide_pages"
DETAILS_DIR = ROOT / "data" / "raw" / "peptide_details"
FAILURES_FILE = ROOT / "data" / "processed" / "failed_peptide_ids.csv"

NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
RANGE_PATTERN = re.compile(rf"^\s*({NUMBER})\s*[-–—]\s*({NUMBER})\s*$")
NUMBER_PATTERN = re.compile(rf"^\s*({NUMBER})\s*$")
THRESHOLD_PATTERN = re.compile(rf"^\s*(<=|>=|<|>)\s*({NUMBER})\s*$")

# Factors convert values to the canonical unit within each dimension.
UNIT_ALIASES = {
    "nm": ("uM", 0.001),
    "nanomolar": ("uM", 0.001),
    "um": ("uM", 1.0),
    "micromolar": ("uM", 1.0),
    "mm": ("uM", 1_000.0),
    "millimolar": ("uM", 1_000.0),
    "m": ("uM", 1_000_000.0),
    "molar": ("uM", 1_000_000.0),
    "ng/ml": ("ug/mL", 0.001),
    "ngml-1": ("ug/mL", 0.001),
    "ug/ml": ("ug/mL", 1.0),
    "mcg/ml": ("ug/mL", 1.0),
    "ugml-1": ("ug/mL", 1.0),
    "mg/l": ("ug/mL", 1.0),
    "mgl-1": ("ug/mL", 1.0),
    "mg/ml": ("ug/mL", 1_000.0),
    "g/l": ("ug/mL", 1_000.0),
    "gml-1": ("ug/mL", 1_000_000.0),
}


def normalize_unit(unit):
    if pd.isna(unit) or not str(unit).strip():
        return "", None, None

    raw = str(unit).strip()
    normalized = unicodedata.normalize("NFKC", raw).lower()
    normalized = normalized.replace("μ", "u").replace("µ", "u")
    normalized = normalized.replace("�", "u")
    normalized = re.sub(r"\s+", "", normalized)

    if normalized in UNIT_ALIASES:
        canonical_unit, factor = UNIT_ALIASES[normalized]
        return raw, canonical_unit, factor
    return raw, None, None


def parse_concentration(value):
    if pd.isna(value) or not str(value).strip():
        return "missing", None, None

    raw = str(value).strip().replace(",", "")
    match = NUMBER_PATTERN.fullmatch(raw)
    if match:
        number = float(match.group(1))
        return "exact", number, number

    match = THRESHOLD_PATTERN.fullmatch(raw)
    if match:
        operator, bound_text = match.groups()
        bound = float(bound_text)
        if operator in ("<", "<="):
            return "censored_upper", None, bound
        return "censored_lower", bound, None

    match = RANGE_PATTERN.fullmatch(raw)
    if match:
        lower, upper = (float(part) for part in match.groups())
        return "range", min(lower, upper), max(lower, upper)

    return "unparsed", None, None


def load_metadata():
    metadata = {}
    expected_ids = set()
    total_count = None

    if not METADATA_DIR.exists():
        return metadata, expected_ids, total_count

    for page_path in sorted(METADATA_DIR.glob("offset_*.json")):
        with page_path.open(encoding="utf-8") as file:
            page = json.load(file)
        total_count = page.get("totalCount", total_count)

        for peptide in page.get("data", []):
            peptide_id = peptide.get("id")
            if peptide_id is None:
                continue

            key = str(peptide_id)
            expected_ids.add(key)
            sequence = str(peptide.get("sequence") or "").strip()
            sequence_length = peptide.get("sequenceLength")
            length = None
            length_source = None

            if sequence_length:
                length = int(sequence_length)
                length_source = "api_sequence_length"
            elif sequence:
                length = len(re.sub(r"\s+", "", sequence))
                length_source = "sequence_character_count"
            else:
                monomers = peptide.get("monomers") or []
                monomer_lengths = []
                for monomer in monomers:
                    monomer_sequence = str(monomer.get("sequence") or "").strip()
                    monomer_length = monomer.get("sequenceLength")
                    if monomer_length:
                        monomer_lengths.append(int(monomer_length))
                    elif monomer_sequence:
                        monomer_lengths.append(
                            len(re.sub(r"\s+", "", monomer_sequence))
                        )
                if monomer_lengths:
                    length = sum(monomer_lengths)
                    length_source = "sum_of_monomer_lengths"

            metadata[key] = {
                "sequence": sequence,
                "length": length,
                "length_source": length_source,
                "complexity": peptide.get("complexity"),
            }

    return metadata, expected_ids, total_count


def load_detail_ids():
    if not DETAILS_DIR.exists():
        return set()
    return {path.stem for path in DETAILS_DIR.glob("*.json")}


def mic_mask(frame):
    measures = frame["activity_measure"].fillna("").astype(str).str.strip()
    return measures.str.upper().str.match(r"^MIC(?:\b|[^A-Z])", na=False)


def write_csv(frame, path):
    frame.to_csv(path, index=False, encoding="utf-8")


def summarize_distribution(values):
    numeric = pd.Series(values, dtype="float64").dropna()
    if numeric.empty:
        return {
            "count": 0,
            "minimum": None,
            "maximum": None,
            "mean": None,
            "median": None,
            "quantiles": {},
        }
    return {
        "count": int(numeric.count()),
        "minimum": float(numeric.min()),
        "maximum": float(numeric.max()),
        "mean": float(numeric.mean()),
        "median": float(numeric.median()),
        "quantiles": {
            str(level): float(numeric.quantile(level))
            for level in (0.25, 0.5, 0.75, 0.9, 0.95)
        },
    }


def make_histograms(unit_values, lengths, output_dir):
    groups = [
        (unit, values)
        for unit, values in sorted(unit_values.items())
        if values
    ]
    if groups:
        columns = min(3, len(groups))
        rows = math.ceil(len(groups) / columns)
        figure, axes = plt.subplots(
            rows,
            columns,
            figsize=(5 * columns, 3.5 * rows),
            squeeze=False,
        )
        for axis, (unit, values) in zip(axes.flat, groups):
            axis.hist(values, bins="auto", edgecolor="black")
            axis.set_title(f"MIC concentration ({unit})")
            axis.set_xlabel(unit)
            axis.set_ylabel("Measurement count")
        for axis in axes.flat[len(groups):]:
            axis.set_visible(False)
        figure.tight_layout()
        figure.savefig(output_dir / "mic_concentration_histograms.png", dpi=160)
        plt.close(figure)

    numeric_lengths = [value for value in lengths if value is not None]
    if numeric_lengths:
        figure, axis = plt.subplots(figsize=(8, 5))
        bins = min(40, max(5, int(math.sqrt(len(numeric_lengths)))))
        axis.hist(numeric_lengths, bins=bins, edgecolor="black")
        axis.set_title("Peptide length distribution")
        axis.set_xlabel("Length (amino acids; multimer lengths are component sums)")
        axis.set_ylabel("Unique peptides")
        figure.tight_layout()
        figure.savefig(output_dir / "peptide_length_histogram.png", dpi=160)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize DBAASP MIC measurements, targets, lengths, and quality."
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Activity CSV; defaults to the final file, or the partial file if final is absent.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for CSV tables, JSON summary, Markdown report, and plots.",
    )
    args = parser.parse_args()

    input_path = args.input
    if input_path is None:
        input_path = (
            DEFAULT_FINAL_INPUT
            if DEFAULT_FINAL_INPUT.exists()
            else DEFAULT_PARTIAL_INPUT
        )
    if not input_path.exists():
        raise FileNotFoundError(f"Activity CSV not found: {input_path}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    activities = pd.read_csv(input_path, dtype={"peptide_id": "string"}, low_memory=False)
    required_columns = {
        "peptide_id",
        "measurement_id",
        "target_species",
        "activity_measure",
        "concentration",
        "unit",
        "sequence",
    }
    missing_columns = required_columns.difference(activities.columns)
    if missing_columns:
        raise ValueError(f"Activity CSV is missing required columns: {sorted(missing_columns)}")

    metadata, expected_ids, api_total_count = load_metadata()
    detail_ids = load_detail_ids()
    if FAILURES_FILE.exists():
        try:
            failures = pd.read_csv(FAILURES_FILE)
        except pd.errors.EmptyDataError:
            failures = pd.DataFrame()
    else:
        failures = pd.DataFrame()
    mic = activities.loc[mic_mask(activities)].copy()
    mic["target_species"] = mic["target_species"].fillna("").astype(str).str.strip()
    mic["concentration_raw"] = mic["concentration"].fillna("").astype(str).str.strip()

    unit_info = mic["unit"].map(normalize_unit)
    mic["unit_raw"] = unit_info.map(lambda item: item[0])
    mic["canonical_unit"] = unit_info.map(lambda item: item[1])
    mic["unit_factor"] = unit_info.map(lambda item: item[2])
    parsed = mic["concentration"].map(parse_concentration)
    mic["concentration_kind"] = parsed.map(lambda item: item[0])
    mic["concentration_lower"] = parsed.map(lambda item: item[1])
    mic["concentration_upper"] = parsed.map(lambda item: item[2])
    mic["normalized_lower"] = mic["concentration_lower"] * mic["unit_factor"]
    mic["normalized_upper"] = mic["concentration_upper"] * mic["unit_factor"]

    # Only exact concentrations are used in numeric distribution summaries.
    mic["normalized_value"] = mic["normalized_lower"].where(
        mic["concentration_kind"].eq("exact")
    )
    mic["canonical_unit"] = mic["canonical_unit"].fillna(
        mic["unit_raw"].replace("", "Unknown unit")
    )
    mic["normalized_value"] = mic.apply(
        lambda row: row["normalized_value"]
        if pd.notna(row["unit_factor"])
        else float(row["concentration_lower"])
        if row["concentration_kind"] == "exact"
        else None,
        axis=1,
    )

    unit_counts = (
        mic.assign(unit_display=mic["unit_raw"].replace("", "(missing)"))
        .groupby("unit_display", dropna=False)
        .size()
        .rename("count")
        .reset_index()
        .sort_values(["count", "unit_display"], ascending=[False, True])
    )
    if len(mic):
        unit_counts["percentage"] = unit_counts["count"] / len(mic) * 100
    else:
        unit_counts["percentage"] = []
    write_csv(unit_counts, args.output_dir / "mic_units.csv")

    value_counts = (
        mic.groupby(
            [
                "concentration_kind",
                "concentration_raw",
                "unit_raw",
                "canonical_unit",
                "normalized_value",
            ],
            dropna=False,
        )
        .size()
        .rename("count")
        .reset_index()
        .sort_values(
            ["canonical_unit", "concentration_kind", "concentration_raw"],
            na_position="last",
        )
    )
    write_csv(value_counts, args.output_dir / "mic_concentration_frequencies.csv")
    normalized_measurements = mic[
        [
            "peptide_id",
            "measurement_id",
            "target_species",
            "concentration_raw",
            "unit_raw",
            "concentration_kind",
            "canonical_unit",
            "normalized_value",
            "normalized_lower",
            "normalized_upper",
        ]
    ]
    write_csv(
        normalized_measurements,
        args.output_dir / "mic_measurements_normalized.csv",
    )

    distribution_rows = []
    histogram_values = {}
    for unit, group in mic.groupby("canonical_unit", dropna=False):
        exact_values = group.loc[
            group["concentration_kind"].eq("exact"), "normalized_value"
        ].dropna()
        summary = summarize_distribution(exact_values)
        summary.update(
            {
                "unit": unit,
                "exact_count": int(group["concentration_kind"].eq("exact").sum()),
                "range_count": int(group["concentration_kind"].eq("range").sum()),
                "censored_count": int(
                    group["concentration_kind"].isin(
                        ["censored_upper", "censored_lower"]
                    ).sum()
                ),
                "unparsed_count": int(group["concentration_kind"].eq("unparsed").sum()),
                "missing_count": int(group["concentration_kind"].eq("missing").sum()),
            }
        )
        distribution_rows.append(summary)
        histogram_values[str(unit)] = exact_values.tolist()
    distribution = pd.DataFrame(distribution_rows)
    write_csv(distribution, args.output_dir / "mic_concentration_statistics.csv")

    target_counts = (
        mic.loc[mic["target_species"].ne("")]
        .groupby("target_species")
        .agg(
            mic_measurement_count=("measurement_id", "size"),
            unique_peptide_count=("peptide_id", "nunique"),
        )
        .reset_index()
        .sort_values(
            ["mic_measurement_count", "unique_peptide_count", "target_species"],
            ascending=[False, False, True],
        )
    )
    write_csv(target_counts, args.output_dir / "mic_targets.csv")

    target_per_peptide = (
        mic.loc[mic["target_species"].ne("")]
        .groupby("peptide_id")
        .agg(
            target_count=("target_species", "nunique"),
            targets=("target_species", lambda values: "; ".join(sorted(set(values)))),
        )
        .reset_index()
    )
    multiple_targets = target_per_peptide.loc[
        target_per_peptide["target_count"].gt(1)
    ].sort_values(["target_count", "peptide_id"], ascending=[False, True])
    write_csv(multiple_targets, args.output_dir / "peptides_tested_against_multiple_targets.csv")

    peptide_rows = []
    activity_sequence = (
        mic[["peptide_id", "sequence"]]
        .dropna(subset=["peptide_id"])
        .drop_duplicates("peptide_id")
        .set_index("peptide_id")["sequence"]
        .to_dict()
    )
    mic_peptide_ids = mic["peptide_id"].dropna().astype(str).unique()
    for peptide_id in mic_peptide_ids:
        info = metadata.get(peptide_id, {})
        activity_seq = activity_sequence.get(peptide_id)
        activity_seq = "" if pd.isna(activity_seq) else str(activity_seq).strip()
        sequence = activity_seq or str(info.get("sequence") or "").strip()
        length = info.get("length")
        length_source = info.get("length_source")
        if length is None and sequence:
            length = len(re.sub(r"\s+", "", sequence))
            length_source = "activity_sequence_character_count"
        peptide_rows.append(
            {
                "peptide_id": peptide_id,
                "sequence": sequence or None,
                "length_aa": length,
                "length_source": length_source or "missing",
                "complexity": info.get("complexity"),
                "has_sequence": bool(sequence),
            }
        )
    peptide_lengths = pd.DataFrame(peptide_rows)
    peptide_lengths = peptide_lengths.sort_values("peptide_id")
    write_csv(peptide_lengths, args.output_dir / "peptide_lengths.csv")

    lengths = pd.to_numeric(peptide_lengths["length_aa"], errors="coerce").dropna()
    length_counts = (
        lengths.value_counts()
        .rename_axis("length_aa")
        .rename("unique_peptide_count")
        .reset_index()
        .sort_values("length_aa")
    )
    write_csv(length_counts, args.output_dir / "peptide_length_frequencies.csv")

    duplicate_rows = int(activities.duplicated(keep="first").sum())
    nonmissing_measurement_ids = activities["measurement_id"].dropna()
    duplicate_measurement_rows = int(nonmissing_measurement_ids.duplicated(keep="first").sum())
    duplicate_measurement_groups = int(
        (nonmissing_measurement_ids.value_counts() > 1).sum()
    )
    repeated_key_columns = [
        "peptide_id",
        "target_species",
        "activity_measure",
        "concentration",
        "unit",
        "medium",
        "reference",
    ]
    existing_key_columns = [column for column in repeated_key_columns if column in activities]
    repeated = activities.loc[
        activities.duplicated(subset=existing_key_columns, keep=False)
    ].sort_values(existing_key_columns)
    write_csv(repeated, args.output_dir / "repeated_measurements.csv")

    mic["unit_dimension"] = mic["canonical_unit"].map(
        lambda unit: "molar" if unit == "uM" else "mass/volume" if unit == "ug/mL" else unit
    )
    inconsistent_pairs = (
        mic.loc[mic["target_species"].ne("")]
        .groupby(["peptide_id", "target_species"])
        .agg(
            unit_count=("unit_dimension", "nunique"),
            units=("unit_raw", lambda values: "; ".join(sorted(set(values)))),
            measurement_count=("measurement_id", "count"),
        )
        .reset_index()
    )
    inconsistent_pairs = inconsistent_pairs.loc[
        inconsistent_pairs["unit_count"].gt(1)
    ].sort_values(["unit_count", "measurement_count"], ascending=[False, False])
    write_csv(inconsistent_pairs, args.output_dir / "inconsistent_unit_pairs.csv")

    normalized_exact = mic.loc[
        mic["concentration_kind"].eq("exact"),
        ["canonical_unit", "normalized_value"],
    ].drop_duplicates()

    matched_detail_ids = detail_ids.intersection(expected_ids) if expected_ids else detail_ids
    coverage = {
        "input_file": str(input_path),
        "activity_rows": int(len(activities)),
        "mic_measurements": int(len(mic)),
        "metadata_records_expected": len(expected_ids) if expected_ids else None,
        "metadata_api_total_count": api_total_count,
        "detail_files_downloaded": len(matched_detail_ids),
        "failed_peptide_ids": int(len(failures)),
    }
    expected = len(expected_ids) or api_total_count
    coverage["detail_download_coverage_percent"] = (
        round(len(matched_detail_ids) / expected * 100, 2) if expected else None
    )
    coverage["extraction_complete"] = bool(
        expected
        and len(matched_detail_ids) >= expected
        and failures.empty
        and input_path.name == "activity_measurements.csv"
    )

    mic_missing = int(mic["concentration_kind"].eq("missing").sum())
    length_mode = (
        int(lengths.mode().iloc[0]) if not lengths.empty else None
    )
    summary = {
        "coverage": coverage,
        "mic": {
            "measurement_count": int(len(mic)),
            "missing_concentration_count": mic_missing,
            "raw_unique_concentration_values": int(
                mic[["concentration_raw", "unit_raw"]].drop_duplicates().shape[0]
            ),
            "exact_numeric_unique_values_by_unit": int(len(normalized_exact)),
            "exact_numeric_count": int(mic["concentration_kind"].eq("exact").sum()),
            "range_count": int(mic["concentration_kind"].eq("range").sum()),
            "censored_count": int(
                mic["concentration_kind"].isin(["censored_upper", "censored_lower"]).sum()
            ),
            "unparsed_count": int(mic["concentration_kind"].eq("unparsed").sum()),
            "exact_numeric_statistics_by_unit": distribution_rows,
        },
        "targets": {
            "unique_nonempty_target_labels": int(target_counts["target_species"].nunique()),
            "missing_target_count": int(mic["target_species"].eq("").sum()),
            "most_tested_target_labels": (
                target_counts.loc[
                    target_counts["unique_peptide_count"].eq(
                        target_counts["unique_peptide_count"].max()
                    ),
                    "target_species",
                ].tolist()
                if not target_counts.empty
                else []
            ),
            "most_tested_target_peptide_count": (
                int(target_counts["unique_peptide_count"].max())
                if not target_counts.empty
                else 0
            ),
            "most_measurement_target_labels": (
                target_counts.loc[
                    target_counts["mic_measurement_count"].eq(
                        target_counts["mic_measurement_count"].max()
                    ),
                    "target_species",
                ].tolist()
                if not target_counts.empty
                else []
            ),
            "most_measurement_target_count": (
                int(target_counts["mic_measurement_count"].max())
                if not target_counts.empty
                else 0
            ),
            "peptides_tested_against_multiple_target_labels": int(len(multiple_targets)),
        },
        "peptide_lengths": {
            "unique_peptides_with_mic": int(len(peptide_lengths)),
            "peptides_with_length": int(lengths.count()),
            "peptides_without_length": int(peptide_lengths["length_aa"].isna().sum()),
            "minimum_aa": int(lengths.min()) if not lengths.empty else None,
            "maximum_aa": int(lengths.max()) if not lengths.empty else None,
            "mean_aa": float(lengths.mean()) if not lengths.empty else None,
            "median_aa": float(lengths.median()) if not lengths.empty else None,
            "mode_aa": length_mode,
            "distribution": summarize_distribution(lengths),
        },
        "quality": {
            "exact_duplicate_rows_additional_copies": duplicate_rows,
            "duplicate_measurement_id_groups": duplicate_measurement_groups,
            "duplicate_measurement_id_rows_additional_copies": duplicate_measurement_rows,
            "repeated_measurement_key_rows": int(len(repeated)),
            "repeated_measurement_key_groups": int(
                len(repeated.drop_duplicates(subset=existing_key_columns))
            ),
            "peptide_target_pairs_with_incompatible_unit_dimensions": int(
                len(inconsistent_pairs)
            ),
            "unit_conversions": (
                "Molar units are standardized to uM; mass/volume units to ug/mL. "
                "No mass-to-molar conversion is attempted because molecular weights "
                "and modification states are not validated."
            ),
        },
        "interpretation": {
            "concentration_column": (
                "MIC numeric statistics use the raw concentration field because the "
                "unit column describes that field. The separate activity field is retained "
                "in the source data and is not assumed to be the concentration."
            ),
            "ranges_and_censoring": (
                "Range and threshold measurements are counted separately and excluded "
                "from exact-value mean, median, quantile, and histogram calculations."
            ),
            "targets": (
                "Target labels are reported exactly as supplied by DBAASP; species-level "
                "labels are not treated as confirmed strains."
            ),
            "multimers": (
                "For a multimer with no parent sequence, length uses the sum of monomer "
                "sequence lengths when available; this is a component-length sum, not a "
                "claim that the chains form one linear sequence."
            ),
        },
    }
    with (args.output_dir / "summary.json").open("w", encoding="utf-8") as file:
        json.dump(summary, file, indent=2, ensure_ascii=False)

    lines = [
        "# DBAASP MIC Dataset Analysis",
        "",
        f"- Input: `{input_path}`",
        f"- Extraction complete: **{coverage['extraction_complete']}**",
        f"- Activity rows analyzed: **{len(activities):,}**",
        f"- MIC measurements analyzed: **{len(mic):,}**",
        f"- Metadata peptide count: **{expected if expected else 'unavailable'}**",
        f"- Detail-file coverage: **{coverage['detail_download_coverage_percent']}%**",
        f"- Failed peptide IDs: **{coverage['failed_peptide_ids']:,}**",
        "",
        "> If extraction is incomplete, this report describes only the current exported rows.",
        "",
        "## MIC concentrations",
        "",
        f"- Missing concentration values: **{mic_missing:,}**",
        f"- Unique raw concentration/unit pairs: **{summary['mic']['raw_unique_concentration_values']:,}**",
        f"- Exact numeric measurements: **{summary['mic']['exact_numeric_count']:,}**",
        f"- Ranges: **{summary['mic']['range_count']:,}**; censored: **{summary['mic']['censored_count']:,}**; unparsed: **{summary['mic']['unparsed_count']:,}**",
        "",
        "Exact-value summaries are grouped by compatible canonical unit. Ranges and censored values are excluded from those numeric summaries.",
        "",
        "| Unit | Exact n | Min | Q1 | Median | Q3 | Q90 | Q95 | Max | Mean | Missing | Range | Censored |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in distribution_rows:
        quantiles = row["quantiles"]
        lines.append(
            "| {unit} | {count} | {minimum} | {q1} | {median} | {q3} | {q90} | {q95} | {maximum} | {mean} | {missing_count} | {range_count} | {censored_count} |".format(
                unit=row["unit"],
                count=row["count"],
                minimum=_format_number(row["minimum"]),
                q1=_format_number(quantiles.get("0.25")),
                median=_format_number(row["median"]),
                q3=_format_number(quantiles.get("0.75")),
                q90=_format_number(quantiles.get("0.9")),
                q95=_format_number(quantiles.get("0.95")),
                maximum=_format_number(row["maximum"]),
                mean=_format_number(row["mean"]),
                missing_count=row["missing_count"],
                range_count=row["range_count"],
                censored_count=row["censored_count"],
            )
        )
    lines.extend(
        [
            "",
            "Unit counts and percentages: see `mic_units.csv`. Per-measurement normalized values: `mic_measurements_normalized.csv`; raw concentration value frequencies: `mic_concentration_frequencies.csv`.",
            "",
            "## Targets",
            "",
            f"- Unique non-empty target labels: **{summary['targets']['unique_nonempty_target_labels']:,}**",
            f"- Most-tested target label(s) by distinct peptides: **{', '.join(summary['targets']['most_tested_target_labels']) or 'none'}** ({summary['targets']['most_tested_target_peptide_count']:,} peptides each)",
            f"- Most MIC measurements: **{', '.join(summary['targets']['most_measurement_target_labels']) or 'none'}** ({summary['targets']['most_measurement_target_count']:,} measurements each)",
            f"- Peptides tested against multiple target labels: **{len(multiple_targets):,}**",
            "",
            "Every target label, MIC measurement count, and unique peptide count is in `mic_targets.csv`. Labels are not automatically parsed into species and strain.",
            "",
            "## Peptide lengths",
            "",
            f"- Unique peptides with MIC records: **{len(peptide_lengths):,}**",
            f"- Peptides with a length: **{int(lengths.count()):,}**",
            f"- Min / max: **{summary['peptide_lengths']['minimum_aa']} / {summary['peptide_lengths']['maximum_aa']} aa**",
            f"- Mean / median / mode: **{_format_number(summary['peptide_lengths']['mean_aa'])} / {_format_number(summary['peptide_lengths']['median_aa'])} / {summary['peptide_lengths']['mode_aa']} aa**",
            "",
            "Lengths are counted once per peptide, not once per measurement. Multimer lengths may be sums of component lengths. Per-peptide provenance is in `peptide_lengths.csv`; frequencies are in `peptide_length_frequencies.csv`.",
            "",
            "## Data quality",
            "",
            f"- Exact duplicate rows (additional copies): **{duplicate_rows:,}**",
            f"- Duplicate measurement-ID groups / additional rows: **{duplicate_measurement_groups:,} / {duplicate_measurement_rows:,}**",
            f"- Repeated measurement-key groups / rows: **{summary['quality']['repeated_measurement_key_groups']:,} / {len(repeated):,}**",
            f"- Peptide/target pairs with inconsistent unit dimensions: **{len(inconsistent_pairs):,}**",
            "",
            summary["quality"]["unit_conversions"],
            "",
            "See `repeated_measurements.csv`, `inconsistent_unit_pairs.csv`, and `summary.json` for details.",
            "",
            "## Plots",
            "",
            "- `mic_concentration_histograms.png`: exact concentrations, separated by canonical unit.",
            "- `peptide_length_histogram.png`: one length per peptide.",
            "",
            "## Interpretation notes",
            "",
            f"- {summary['interpretation']['concentration_column']}",
            f"- {summary['interpretation']['ranges_and_censoring']}",
            f"- {summary['interpretation']['targets']}",
            f"- {summary['interpretation']['multimers']}",
        ]
    )
    (args.output_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    make_histograms(histogram_values, lengths.tolist(), args.output_dir)
    print(f"Analysis complete: {len(mic)} MIC measurements from {input_path}")
    print(f"Extraction complete: {coverage['extraction_complete']}")
    print(f"Report: {args.output_dir / 'report.md'}")


def _format_number(value):
    if value is None or pd.isna(value):
        return "n/a"
    return f"{value:.4g}"


if __name__ == "__main__":
    main()
