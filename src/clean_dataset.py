"""Audit and clean DBAASP activity measurements without discarding provenance.

Run from repository root:
    python src/clean_dataset.py

Expected input:
    data/processed/activity_measurements.csv (full extraction)
    or data/processed/activity_sample.csv (sample fallback)

Metadata is read from data/raw/peptide_pages/offset_*.json and/or
 data/raw/peptide_sample.json. Detailed records in data/raw/peptide_details/*.json
 are also used to recover sequences and monomer information.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
REPORTS = PROCESSED / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)


def read_json(path: Path):
    try:
        with path.open(encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def peptide_records():
    """Collect parent peptide metadata; prefer details over listing records."""
    records = {}
    for path in sorted((RAW / "peptide_pages").glob("offset_*.json")):
        obj = read_json(path) or {}
        for p in obj.get("data", []):
            if p.get("id") is not None:
                records[str(p["id"])] = p

    sample = read_json(RAW / "peptide_sample.json") or {}
    for p in sample.get("data", []):
        if p.get("id") is not None:
            records.setdefault(str(p["id"]), p)

    for path in sorted((RAW / "peptide_details").glob("*.json")):
        obj = read_json(path)
        if not isinstance(obj, dict):
            continue
        p = obj.get("data", obj)
        if isinstance(p, list):
            p = p[0] if p else {}
        if isinstance(p, dict) and p.get("id") is not None:
            records[str(p["id"])] = p
    return records


def sequence_info(p):
    seq = str(p.get("sequence") or "").strip().upper()
    monomers = p.get("monomers") or []
    monomer_seqs = [str(m.get("sequence") or "").strip().upper()
                    for m in monomers if str(m.get("sequence") or "").strip()]
    # Remove whitespace only. Do not silently delete non-canonical symbols.
    seq = re.sub(r"\s+", "", seq)
    monomer_seqs = [re.sub(r"\s+", "", s) for s in monomer_seqs]
    if not seq and len(monomer_seqs) == 1:
        seq = monomer_seqs[0]
    valid = bool(seq) and bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+", seq))
    if seq:
        length = len(seq)
    else:
        length = sum(map(len, monomer_seqs)) if monomer_seqs else np.nan
    return {
        "sequence": seq if seq else np.nan,
        "monomer_sequences": "|".join(monomer_seqs) if monomer_seqs else np.nan,
        "sequence_length": length,
        "sequence_valid_alphabet": valid if seq else False,
        "complexity": p.get("complexity"),
        "synthesis_type": p.get("synthesisType"),
        "peptide_name_metadata": p.get("name"),
    }


def parse_numeric_or_range(value):
    """Return low/high/operator from values like 8, '0.73-11.68', '<= 4'."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan, np.nan, "missing"
    s = str(value).strip().replace(",", ".").replace("−", "-").replace("–", "-")
    if not s:
        return np.nan, np.nan, "missing"
    s = re.sub(r"\s+", "", s)
    m = re.fullmatch(r"(<=|>=|<|>)?([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", s)
    if m:
        op, num = m.groups(); n = float(num)
        return n, n, op or "exact"
    m = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)-([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", s)
    if m:
        a, b = map(float, m.groups())
        return min(a, b), max(a, b), "range"
    return np.nan, np.nan, "unparsed"


def unit_kind(unit):
    u = str(unit or "").strip().lower().replace("μ", "µ").replace("µ", "u")
    u = u.replace(" ", "")
    if u in {"m", "mol/l", "mol·l-1", "molar"}: return "M", 1e6
    if u in {"mm", "mmol/l", "millimolar"}: return "M", 1e3
    if u in {"um", "µm", "micromolar", "µmol/l", "umol/l"}: return "M", 1.0
    if u in {"nm", "nmol/l", "nanomolar"}: return "M", 1e-3
    if u in {"pm", "pmol/l", "picomolar"}: return "M", 1e-6
    if u in {"mg/ml", "mgml-1"}: return "mass", 1000.0  # mg/mL -> ug/mL
    if u in {"ug/ml", "µg/ml", "mcg/ml", "ugml-1"}: return "mass", 1.0
    if u in {"ng/ml", "ngml-1"}: return "mass", 1e-3
    return "unknown", np.nan


def molecular_weight_g_mol(sequence):
    """Estimate unmodified linear peptide MW; returns NaN if sequence is unsuitable."""
    if not isinstance(sequence, str) or not sequence or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence):
        return np.nan
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        return float(ProteinAnalysis(sequence).molecular_weight())
    except ImportError:
        return np.nan


def main():
    full = PROCESSED / "activity_measurements.csv"
    sample = PROCESSED / "activity_sample.csv"
    source = full if full.exists() else sample
    if not source.exists():
        raise FileNotFoundError("No activity CSV found. Run src/data_collection.py and src/extract_activities.py first.")

    df = pd.read_csv(source, dtype={"peptide_id": "string", "measurement_id": "string"})
    df.columns = [str(c).strip().lower() for c in df.columns]
    if "activity_measure" not in df.columns:
        raise ValueError("Missing activity_measure column; inspect the DBAASP export schema.")

    # Metadata join restores sequences omitted in target-activity rows (notably multimers).
    metadata = peptide_records()
    meta_rows = []
    for pid, p in metadata.items():
        meta_rows.append({"peptide_id": str(pid), **sequence_info(p)})
    if meta_rows:
        meta = pd.DataFrame(meta_rows).drop_duplicates("peptide_id", keep="last")
        if "peptide_id" in df.columns:
            df["peptide_id"] = df["peptide_id"].astype("string")
            df = df.merge(meta, on="peptide_id", how="left", suffixes=("", "_metadata"))
            # The activity export may contain an empty parent sequence for multimers.
            if "sequence_metadata" in df.columns:
                df["sequence"] = df["sequence"].where(df["sequence"].notna() & df["sequence"].astype(str).str.strip().ne(""), df["sequence_metadata"]) if "sequence" in df.columns else df["sequence_metadata"]
            if "peptide_name_metadata" in df.columns and "peptide_name" in df.columns:
                df["peptide_name"] = df["peptide_name"].where(df["peptide_name"].notna() & df["peptide_name"].astype(str).str.strip().ne(""), df["peptide_name_metadata"])

    # Preserve original values. Parse activity/concentration separately; do not assume they mean the same thing.
    value_col = "concentration" if "concentration" in df.columns else "activity"
    parsed = df[value_col].apply(parse_numeric_or_range)
    df["mic_source_value"] = df[value_col]
    df["value_low_original_unit"] = parsed.apply(lambda x: x[0])
    df["value_high_original_unit"] = parsed.apply(lambda x: x[1])
    df["value_operator"] = parsed.apply(lambda x: x[2])

    if "unit" not in df.columns:
        df["unit"] = np.nan
    kinds = df["unit"].apply(unit_kind)
    df["unit_family"] = kinds.apply(lambda x: x[0])
    df["unit_to_uM_factor"] = kinds.apply(lambda x: x[1] if x[0] == "M" else np.nan)
    df["value_low_uM"] = np.where(df["unit_family"].eq("M"), df["value_low_original_unit"] * df["unit_to_uM_factor"], np.nan)
    df["value_high_uM"] = np.where(df["unit_family"].eq("M"), df["value_high_original_unit"] * df["unit_to_uM_factor"], np.nan)

    # Mass-to-molar conversion only for canonical single-chain sequences. This is an estimate and flagged.
    if "sequence" in df.columns:
        df["estimated_mw_g_mol"] = df["sequence"].apply(molecular_weight_g_mol)
    else:
        df["estimated_mw_g_mol"] = np.nan
    mass_mask = df["unit_family"].eq("mass") & df["estimated_mw_g_mol"].notna()
    # Values are represented as ug/mL after applying the mass-unit factor.
    mass_factor = kinds.apply(lambda x: x[1] if x[0] == "mass" else np.nan)
    mass_low_ug_ml = df["value_low_original_unit"] * mass_factor
    mass_high_ug_ml = df["value_high_original_unit"] * mass_factor
    df.loc[mass_mask, "value_low_uM"] = mass_low_ug_ml[mass_mask] * 1000 / df.loc[mass_mask, "estimated_mw_g_mol"]
    df.loc[mass_mask, "value_high_uM"] = mass_high_ug_ml[mass_mask] * 1000 / df.loc[mass_mask, "estimated_mw_g_mol"]
    df["unit_conversion_status"] = np.select(
        [df["unit_family"].eq("M"), mass_mask, df["unit_family"].eq("mass"), df["unit_family"].eq("unknown")],
        ["converted_from_molar_unit", "estimated_from_mass_unit_and_linear_sequence_MW", "mass_unit_needs_MW_or_manual_review", "unknown_unit_review"],
        default="unparsed_value_or_missing_unit",
    )

    # Activity class is normalized but original label is preserved.
    df["activity_measure_normalized"] = df["activity_measure"].astype("string").str.strip().str.upper()
    mic = df[df["activity_measure_normalized"].eq("MIC")].copy()
    mic.to_csv(PROCESSED / "mic_measurements_clean.csv", index=False)
    df.to_csv(PROCESSED / "activity_measurements_audited.csv", index=False)

    # Reproducible summaries.
    unit_summary = df.groupby(["activity_measure_normalized", "unit"], dropna=False).size().rename("records").reset_index()
    unit_summary.to_csv(REPORTS / "unit_counts.csv", index=False)
    strain_col = "target_species" if "target_species" in mic.columns else None
    if strain_col:
        strain_summary = mic.groupby(strain_col, dropna=False).agg(
            mic_measurements=("peptide_id", "size"),
            distinct_peptides=("peptide_id", "nunique"),
            distinct_measurement_ids=("measurement_id", "nunique") if "measurement_id" in mic.columns else ("peptide_id", "size"),
        ).reset_index().sort_values(["distinct_peptides", "mic_measurements"], ascending=False)
        strain_summary.to_csv(REPORTS / "strain_coverage.csv", index=False)
    if "sequence_length" in mic.columns:
        mic[["peptide_id", "sequence_length"]].drop_duplicates().to_csv(REPORTS / "peptide_lengths.csv", index=False)

    summary = {
        "input_file": str(source.relative_to(ROOT)),
        "input_activity_rows": int(len(df)),
        "unique_peptides_in_activity_rows": int(df["peptide_id"].nunique()) if "peptide_id" in df else None,
        "MIC_rows": int(len(mic)),
        "non_MIC_rows_retained_in_audit": int(len(df) - len(mic)),
        "unique_MIC_units": int(mic["unit"].nunique(dropna=True)) if "unit" in mic else None,
        "unknown_or_unparsed_units_MIC": int(mic["unit_family"].eq("unknown").sum()),
        "MIC_values_not_parsed": int(mic["value_operator"].isin(["unparsed", "missing"]).sum()),
        "MIC_rows_with_uM_conversion": int(mic["value_low_uM"].notna().sum()),
        "sequence_representation_missing_MIC_rows": int((mic["sequence"].isna() & mic.get("monomer_sequences", pd.Series(index=mic.index, dtype=object)).isna()).sum()) if "sequence" in mic else None,
        "sequence_length_min": float(mic["sequence_length"].min()) if "sequence_length" in mic and mic["sequence_length"].notna().any() else None,
        "sequence_length_max": float(mic["sequence_length"].max()) if "sequence_length" in mic and mic["sequence_length"].notna().any() else None,
        "note": "If the input is activity_sample.csv, these are sample-only statistics, not DBAASP-wide statistics. Mass-to-molar conversions are estimates for canonical linear sequences only.",
    }
    with (REPORTS / "cleaning_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2))
    print("\nWrote:")
    print(" -", PROCESSED / "activity_measurements_audited.csv")
    print(" -", PROCESSED / "mic_measurements_clean.csv")
    print(" -", REPORTS)

if __name__ == "__main__":
    main()
