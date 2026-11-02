#!/usr/bin/env python3
"""Prepare step 1 data for gallstone-related disease prediction.

This script deliberately stops before predictor selection, imputation, or modelling.
It creates one record per first study episode, keeps SEH and outpatient cohorts
separate, links only same-day investigations to the index visit, derives the
three-diagnosis outcome, and reports variable availability and data-quality issues.
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "OPERATE data" / "Buikpijn_na_bariatrische_chirurg_export_20260622.csv"
OUTPUT = ROOT / "results" / "step1_gallstone"

CASTOR_MISSING = {
    "-95": "__MEASUREMENT_FAILED__",
    "-96": "__NOT_APPLICABLE__",
    "-97": "__NOT_ASKED__",
    "-98": "__ASKED_UNKNOWN__",
    "-99": "__NOT_DONE__",
}

DATE_MISSING = {
    "01-01-2995": "__MEASUREMENT_FAILED__",
    "01-01-2996": "__NOT_APPLICABLE__",
    "01-01-2997": "__NOT_ASKED__",
    "01-01-2998": "__ASKED_UNKNOWN__",
    "01-01-2999": "__NOT_DONE__",
}

TIME_MISSING = {
    "00:00:35": "__MEASUREMENT_FAILED__",
    "00:00:36": "__NOT_APPLICABLE__",
    "00:00:37": "__NOT_ASKED__",
    "00:00:38": "__ASKED_UNKNOWN__",
    "00:00:39": "__NOT_DONE__",
}

TEXT_MISSING = {
    "niet omschreven": "__NOT_DESCRIBED__",
    "niet beschreven": "__NOT_DESCRIBED__",
    "-": "__NOT_DESCRIBED__",
    "--": "__NOT_DESCRIBED__",
    "nvt": "__NOT_APPLICABLE__",
    "n.v.t.": "__NOT_APPLICABLE__",
    "niet van toepassing": "__NOT_APPLICABLE__",
    "onbekend": "__ASKED_UNKNOWN__",
    "unknown": "__ASKED_UNKNOWN__",
    "niet verricht": "__NOT_DONE__",
    "niet gedaan": "__NOT_DONE__",
    "geen antwoord ingevuld": "__NOT_RECORDED__",
}

SPECIAL_TOKENS = set(CASTOR_MISSING.values()) | set(TEXT_MISSING.values())

DEF_DIAGNOSIS = {
    "1": "appendicitis",
    "2": "diverticulitis",
    "3": "ulcuslijden",
    "4": "interne_herniatie",
    "5": "cholecystolithiasis",
    "6": "cholecystitis",
    "7": "lekkage_perforatie",
    "8": "coprostase",
    "9": "ibs",
    "10": "buikwandneuralgie",
    "11": "stenose_gj",
    "12": "stenose_jj",
    "13": "ileus_streng",
    "14": "bezoar",
    "15": "buikwandbreuk",
    "16": "urolithiasis",
    "17": "pancreatitis",
    "18": "onbegrepen_klachten",
    "19": "overig",
}

OLD_DIAGNOSIS = {
    "1": "symptomatisch_galsteenlijden",
    "2": "interne_herniatie",
    "3": "cholecystitis",
    "4": "diverticulitis",
    "5": "appendicitis",
    "6": "ulcus",
    "7": "stenose",
    "8": "urolithiasis",
    "9": "lekkage",
    "10": "anders",
}

OUTCOME_RAW_COLUMNS = [
    "FU_einde_datum",
    "FU_einde_hoe",
    "FU_einde_diagnose",
    "FU_einde_diagnose_definitief",
    "FU_einde_diagnose_definitief_ordinaal",
    "FU_einde_definitieve_diagnose_1_DEF_1",
    "FU_einde_def_diagnose_overig_1",
    "definitieve_diagnose_datum",
    "FU_einde_diagnose_2",
    "FU_einde_diagnose_definitief_ordinaal_2",
    "FU_einde_def_diagnose_overig_2",
    "definitieve_diagnose_datum_2",
    "FU_einde_definitieve_diagnose_3",
    "FU_einde_definitieve_diagnose_3_ordinaal",
    "FU_einde_definitieve_diagnose_3_overig",
    "Definiteve_diagnose_datum_3",
    "FU_einde_persisterende_klachten",
]

DX_NEW_COLUMNS = [
    "FU_einde_definitieve_diagnose_1_DEF_1",
    "FU_einde_diagnose_definitief_ordinaal_2",
    "FU_einde_definitieve_diagnose_3_ordinaal",
]

DX_TEXT_COLUMNS = [
    "FU_einde_diagnose_definitief",
    "FU_einde_def_diagnose_overig_1",
    "FU_einde_def_diagnose_overig_2",
    "FU_einde_definitieve_diagnose_3_overig",
]

LEGACY_CLEAR_GALLSTONE_TEXT = re.compile(
    r"(?:cholodocholith|choledocholith|cholecystolith|cholelith|"
    r"cholecystitis|galsteen|biliaire pancreatitis|"
    r"cholangiosepsis bij cholecystitis|"
    r"cholangitis.{0,40}(?:chole|choledo|galsteen))",
    flags=re.IGNORECASE,
)

CHOLECYSTOLITHIASIS_FINAL_TEXT = re.compile(
    r"(?:cholodocholith|choledocholith|cholecystolith|cholelith|"
    r"symptomatisch(?:e)?\s+galsteenlijden|galwegste(?:en|nen)|"
    r"choledochusste(?:en|nen))",
    flags=re.IGNORECASE,
)

CHOLECYSTITIS_FINAL_TEXT = re.compile(r"\bcholecystitis\b", flags=re.IGNORECASE)

OTHER_POSSIBLE_GALLSTONE_RELATED_TEXT = re.compile(
    r"(?:biliaire\s+pancreatitis|pancreatitis.{0,50}galste|"
    r"\bcholangitis\b|\bcholangiosepsis\b)",
    flags=re.IGNORECASE,
)

ERCP_STONE_TEXT = re.compile(
    r"(?:cholodocholith|choledocholith|cholecystolith|cholelith|"
    r"galwegste(?:en|nen)|choledochusste(?:en|nen)|\bste(?:en|nen)\b|concrement)",
    flags=re.IGNORECASE,
)

CHOLE_OPERATION_TEXT = re.compile(
    r"(?:cholecystect|cholecystolithiot|cholelithiot|lap(?:aroscopisch|aroscopische|\.)?\s*chol)",
    flags=re.IGNORECASE,
)

ECHO_CHOLECYSTITIS_EVIDENCE = [
    (f"Echo{number}_datum", f"Echo{number}_cholecystitis") for number in range(1, 6)
]

CT_CHOLECYSTITIS_EVIDENCE = [
    ("CTabd1_datum", "CTabd1_CC"),
    ("CT_basis_date_1", "CT_basis_CC_1"),
    ("CT_basis_date_2", "CT_basis_CC_2"),
    ("CT_basis_date_3", "CT_basis_CC_3"),
    ("CT_basis_date_4", "CT_basis_CC_4"),
    ("CT_basis_date_4_1", "CT_basis_CC_4_1"),
]

SURGERY_EVIDENCE = [
    (
        "Her_op",
        "Her_op_datum",
        "Her_op_soort",
        "Her_Op_overig",
        "Diagnose_heroperatie_definitief_ordinaal",
        "Diagnose_heroperatie_definitief",
    ),
    (
        "Her_op_1",
        "Her_op_datum_1",
        "Her_op_soort_1",
        "Her_Op_overig_1",
        "Diagnose_heroperatie_2_definitief_ordinaal",
        "Diagnose_heroperatie_definitief_1",
    ),
    (
        "Her_op_2",
        "Her_op_datum_2",
        "Her_op_soort_2",
        "Her_Op_overig_2",
        "Diagnose_heroperatie_3_definitief_ordinaal",
        "Diagnose_heroperatie_definitief_2",
    ),
]

ERCP_EVIDENCE = [
    (
        f"AO_overig_ERCP_datum{suffix}",
        f"AO_overig_ERCP_afwijkend{suffix}",
        f"AO_overig_ERCP_verslag{suffix}",
        f"AO_overig_ERCP_diagnose{suffix}",
    )
    for suffix in ["", "_1", "_2", "_3"]
]

LAB_UNITS = {
    "HB": "mmol/L",
    "HT": "L/L",
    "ery": "10^12/L",
    "MCV": "fL",
    "trombos": "10^9/L",
    "leukos": "10^9/L",
    "natrium": "mmol/L",
    "kalium": "mmol/L",
    "ureum": "not_specified_in_crf",
    "kreat": "umol/L",
    "EGFR": "mL/min/1.73m2",
    "bilirubin_total": "umol/L",
    "AF": "U/L",
    "gGT": "U/L",
    "ASAT": "U/L",
    "ALAT": "U/L",
    "LD": "not_specified_in_crf",
    "lipase": "not_specified_in_crf",
    "calcium": "mmol/L",
    "CRP": "mg/L",
    "glucose": "not_specified_in_crf",
}

SANITY_RANGES = {
    "age_years_at_index": (18, 100),
    "years_since_bariatric_surgery": (0, 40),
    "SEH_LO_Temp_numeric": (30, 45),
    "SEH_LO_pols_numeric": (20, 250),
    "SEH_LO_SRR_numeric": (50, 300),
    "SEH_LO_DRR_numeric": (20, 200),
    "SEH_LO_Sat_numeric": (50, 100),
    "SEH_LO_AF_numeric": (4, 80),
    "SEH_LO_gewicht_numeric": (20, 400),
    "SEH_LO_TWL_numeric": (-100, 100),
    "SEH_LO_BMI_numeric": (10, 100),
    "POLI_LO_temp_numeric": (30, 45),
    "POLI_LO_pols_numeric": (20, 250),
    "POLI_LO_SRR_numeric": (50, 300),
    "POLI_LO_DRR_numeric": (20, 200),
    "POLI_LO_sat_numeric": (50, 100),
    "POLI_LO_AF_numeric": (4, 80),
    "POLI_LO_gewicht_numeric": (20, 400),
    "POLI_LO_TWL_numeric": (-100, 100),
    "POLI_LO_BMI_numeric": (10, 100),
    "index_lab_HB_numeric": (2, 15),
    "index_lab_HT_numeric": (0.1, 0.8),
    "index_lab_ery_numeric": (1, 10),
    "index_lab_MCV_numeric": (40, 150),
    "index_lab_trombos_numeric": (0, 1500),
    "index_lab_leukos_numeric": (0, 100),
    "index_lab_natrium_numeric": (100, 180),
    "index_lab_kalium_numeric": (1, 10),
    "index_lab_ureum_numeric": (0, 100),
    "index_lab_kreat_numeric": (0, 3000),
    "index_lab_EGFR_numeric": (0, 200),
    "index_lab_bilirubin_total_numeric": (0, 1000),
    "index_lab_AF_numeric": (0, 3000),
    "index_lab_gGT_numeric": (0, 5000),
    "index_lab_ASAT_numeric": (0, 10000),
    "index_lab_ALAT_numeric": (0, 10000),
    "index_lab_LD_numeric": (0, 10000),
    "index_lab_lipase_numeric": (0, 50000),
    "index_lab_calcium_numeric": (1, 4),
    "index_lab_CRP_numeric": (0, 700),
    "index_lab_glucose_numeric": (0, 50),
}


def clean_exact_text(value: object) -> object:
    if pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if text == "":
        return pd.NA
    if text in CASTOR_MISSING:
        return CASTOR_MISSING[text]
    if text in DATE_MISSING:
        return DATE_MISSING[text]
    if text in TIME_MISSING:
        return TIME_MISSING[text]
    simplified = text.strip("'\"").strip()
    # In this CRF, uppercase/lowercase "NO" was used as the abbreviation for
    # "Niet Omschreven" in narrative fields. Title-case "No" and Dutch
    # ja/nee are ordinary binary answers.
    if simplified in {"NO", "no"}:
        return "__NOT_DESCRIBED__"
    if simplified in {"Ja", "ja", "Yes", "yes"}:
        return "yes"
    if simplified in {"Nee", "nee", "No"}:
        return "no"
    normalized = re.sub(r"\s+", " ", simplified).casefold()
    return TEXT_MISSING.get(normalized, text)


def is_date_column(name: str) -> bool:
    low = name.casefold()
    return (
        name in {"Geboortedatum", "Participant Creation Date"}
        or "datum" in low
        or re.search(r"(?:^|_)date(?:_|$)", low) is not None
    )


def is_time_column(name: str) -> bool:
    low = name.casefold()
    return low.endswith("tijd") or "_tijd" in low or "_time" in low


def normalize_date(value: object) -> object:
    value = clean_exact_text(value)
    if pd.isna(value) or value in SPECIAL_TOKENS:
        return value
    text = str(value)
    parsed = (
        pd.to_datetime(text, format="%Y-%m-%d", errors="coerce")
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", text)
        else pd.to_datetime(text, dayfirst=True, errors="coerce")
    )
    if pd.isna(parsed):
        return value
    return parsed.strftime("%Y-%m-%d")


def normalize_time(value: object) -> object:
    value = clean_exact_text(value)
    if pd.isna(value) or value in SPECIAL_TOKENS:
        return value
    parsed = pd.to_datetime(str(value), errors="coerce")
    if pd.isna(parsed):
        return value
    return parsed.strftime("%H:%M:%S")


def clean_dataframe(frame: pd.DataFrame) -> pd.DataFrame:
    cleaned = frame.copy()
    for column in cleaned.columns:
        if is_date_column(column):
            cleaned[column] = cleaned[column].map(normalize_date)
        elif is_time_column(column):
            cleaned[column] = cleaned[column].map(normalize_time)
        else:
            cleaned[column] = cleaned[column].map(clean_exact_text)
    return cleaned


def parse_date_series(series: pd.Series) -> pd.Series:
    masked = series.where(~series.astype("string").isin(SPECIAL_TOKENS))
    return pd.to_datetime(masked, dayfirst=True, errors="coerce")


def date_equal(left: object, right: object) -> bool:
    left_date = pd.to_datetime(left, dayfirst=True, errors="coerce")
    right_date = pd.to_datetime(right, dayfirst=True, errors="coerce")
    return bool(pd.notna(left_date) and pd.notna(right_date) and left_date.date() == right_date.date())


def parse_numeric(value: object) -> tuple[float | None, str | None]:
    value = clean_exact_text(value)
    if pd.isna(value) or value in SPECIAL_TOKENS:
        return None, None
    text = str(value).strip().replace(",", ".")
    match = re.fullmatch(r"\s*(<=|>=|<|>)?\s*(-?\d+(?:\.\d+)?)\s*", text)
    if not match:
        return None, None
    return float(match.group(2)), match.group(1) or "="


def canonical_lab_field(name: str, selector: str) -> str:
    if name == selector:
        return "record_present"
    if re.fullmatch(r"Lab\d+_datum", name):
        return "date"
    normal = re.match(r"Lab\d+_(.+)$", name)
    if normal:
        field = normal.group(1)
    elif name.startswith("Lab_basis_date"):
        return "date"
    elif name.startswith("Lab_basis_"):
        field = name[len("Lab_basis_") :]
        suffix = selector[len("Lab_basis") :]
        if suffix and field.endswith(suffix):
            field = field[: -len(suffix)]
    elif name.startswith("Bevindingen_nav_Bloed"):
        return "bevindingen"
    else:
        field = name
    synonyms = {
        "bili_total": "bilirubin_total",
        "bili_totaal": "bilirubin_total",
        "bilirubin_total": "bilirubin_total",
        "af": "AF",
        "ggt": "gGT",
        "glucoseorien": "glucose",
        "glucose": "glucose",
    }
    return synonyms.get(field, field)


def canonical_first_group(name: str, prefix: str, selector: str) -> str:
    if name == selector:
        return "record_present"
    if name.startswith(prefix):
        field = name[len(prefix) :]
    else:
        field = name
    field = field.lstrip("_")
    if field.casefold() in {"datum", "date"}:
        return "date"
    if field.casefold() in {"tijd", "time"}:
        return "time"
    return field


def positional_groups(
    columns: list[str],
    starts: list[str],
    final_end: str,
    first_prefix: str,
    canonicalizer: Callable[[str, str], str] | None = None,
) -> list[dict[str, object]]:
    positions = {name: columns.index(name) for name in starts + [final_end]}
    first_start = positions[starts[0]]
    first_end = positions[starts[1]] if len(starts) > 1 else positions[final_end]
    first_names = columns[first_start:first_end]
    if canonicalizer is None:
        canonical = [canonical_first_group(name, first_prefix, starts[0]) for name in first_names]
    else:
        canonical = [canonicalizer(name, starts[0]) for name in first_names]
    groups: list[dict[str, object]] = []
    for idx, selector in enumerate(starts):
        start = positions[selector]
        end = positions[starts[idx + 1]] if idx + 1 < len(starts) else positions[final_end]
        source_names = columns[start:end]
        if canonicalizer is not None:
            names = [canonicalizer(name, selector) for name in source_names]
        else:
            names = canonical[: len(source_names)]
        date_candidates = [name for name, canon in zip(source_names, names) if canon == "date"]
        if not date_candidates:
            raise ValueError(f"No date column found for {selector}")
        groups.append(
            {
                "selector": selector,
                "date": date_candidates[0],
                "mapping": list(zip(names, source_names)),
            }
        )
    return groups


def add_index_records(
    raw: pd.DataFrame,
    cohort: pd.DataFrame,
    index_date_column: str,
    groups: list[dict[str, object]],
    output_prefix: str,
    performed_raw_column: str | None,
    performed_yes_value: str = "2",
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    canonical_fields = sorted({canon for group in groups for canon, _ in group["mapping"]})
    for row_id in cohort.index:
        source_row = raw.loc[row_id]
        index_date = source_row.get(index_date_column)
        matches = [group for group in groups if date_equal(source_row.get(group["date"]), index_date)]
        performed_value = source_row.get(performed_raw_column) if performed_raw_column else None
        performed = str(performed_value).strip() == performed_yes_value if pd.notna(performed_value) else False
        explicitly_not_performed = (
            pd.notna(performed_value) and str(performed_value).strip() == "1" if performed_raw_column else False
        )
        index_missing = pd.isna(pd.to_datetime(index_date, dayfirst=True, errors="coerce"))
        if index_missing:
            status = "index_date_missing"
        elif matches and performed:
            status = "matched" if len(matches) == 1 else "multiple_same_date_records"
        elif matches and explicitly_not_performed:
            status = "record_match_flag_no_conflict"
        elif matches:
            status = "record_match_performance_unknown"
        elif performed:
            status = "performed_or_requested_no_same_day_detail"
        elif explicitly_not_performed:
            status = "not_performed_or_not_requested"
        else:
            status = "performance_unknown"
        out: dict[str, object] = {
            f"{output_prefix}_link_status": status,
            f"{output_prefix}_same_day_record_count": len(matches),
        }
        if not matches and performed:
            for field in canonical_fields:
                if field not in {"record_present", "date", "time"}:
                    out[f"{output_prefix}_{field}"] = "__PERFORMED_DETAILS_UNAVAILABLE__"
        for number, group in enumerate(matches, start=1):
            repeat_prefix = output_prefix if number == 1 else f"{output_prefix}_repeat{number}"
            out[f"{repeat_prefix}_source_instance"] = group["selector"]
            for canonical, source_name in group["mapping"]:
                out[f"{repeat_prefix}_{canonical}"] = source_row.get(source_name)
        rows.append(out)
    detail = pd.DataFrame(rows, index=cohort.index)
    return pd.concat([cohort, detail], axis=1)


def add_other_diagnostic_blocks(raw: pd.DataFrame, cohort: pd.DataFrame, index_date_column: str) -> pd.DataFrame:
    suffixes = ["", "_1", "_2", "_3"]
    specs = {
        "colon": ("colon_datum", ["colon_afwijkend", "colon_verslag", "colon_diagnose"]),
        "mrcp": ("MRCP_datum", ["MRCP_afwijkend", "MRCP_verslag", "MRCP_diagnose"]),
        "ercp": ("ERCP_datum", ["ERCP_afwijkend", "ERCP_verslag", "ERCP_diagnose"]),
        "mri": ("MRI_datum", ["MRI_afwijkend", "MRI_verslag", "MRI_diagnose"]),
    }
    rows: list[dict[str, object]] = []
    for row_id in cohort.index:
        source_row = raw.loc[row_id]
        index_date = source_row.get(index_date_column)
        out: dict[str, object] = {}
        for modality, (date_stem, field_stems) in specs.items():
            matches: list[str] = []
            for suffix in suffixes:
                date_col = f"AO_overig_{date_stem}{suffix}"
                if date_col in raw.columns and date_equal(source_row.get(date_col), index_date):
                    matches.append(suffix)
            out[f"index_{modality}_same_day_record_count"] = len(matches)
            out[f"index_{modality}_link_status"] = "matched" if matches else "no_same_day_record"
            for number, suffix in enumerate(matches, start=1):
                prefix = f"index_{modality}" if number == 1 else f"index_{modality}_repeat{number}"
                out[f"{prefix}_date"] = source_row.get(f"AO_overig_{date_stem}{suffix}")
                for stem in field_stems:
                    out[f"{prefix}_{stem.split('_', 1)[-1]}"] = source_row.get(f"AO_overig_{stem}{suffix}")
        rows.append(out)
    return pd.concat([cohort, pd.DataFrame(rows, index=cohort.index)], axis=1)


def add_same_day_survey(raw: pd.DataFrame, cohort: pd.DataFrame, index_date_column: str) -> pd.DataFrame:
    survey_columns = list(raw.columns[raw.columns.get_loc("BasisVL_datum") : raw.columns.get_loc("vpo_2wk_datum")])
    rows: list[dict[str, object]] = []
    for row_id in cohort.index:
        row = raw.loc[row_id]
        matched = date_equal(row.get("BasisVL_datum"), row.get(index_date_column))
        values: dict[str, object] = {"index_survey_link_status": "matched" if matched else "no_same_day_survey"}
        if matched:
            values.update({f"index_survey_{column}": row.get(column) for column in survey_columns})
        rows.append(values)
    return pd.concat([cohort, pd.DataFrame(rows, index=cohort.index)], axis=1)


def evidence_on_or_after_index(
    raw_cohort: pd.DataFrame,
    index_dates: pd.Series,
    date_column: str,
    condition: pd.Series,
) -> pd.Series:
    """Return dated evidence occurring on or after the first presentation."""
    evidence_dates = parse_date_series(raw_cohort[date_column])
    return condition.fillna(False) & evidence_dates.notna() & index_dates.notna() & evidence_dates.ge(index_dates)


def derive_outcome(raw_cohort: pd.DataFrame, setting: str) -> pd.DataFrame:
    """Derive the outcome using the supplied definitive-diagnosis table.

    A final diagnosis label remains necessary. Cholecystitis additionally needs
    positive ultrasound, CT, or relevant surgical evidence. Cholecystolithiasis
    additionally needs laparoscopic cholecystectomy, laparoscopic
    cholecystolithotomy, or a stone-positive ERCP. Procedures without a matching
    final diagnosis label never create an event.
    """
    new_codes = raw_cohort[DX_NEW_COLUMNS].fillna("").astype(str).apply(lambda col: col.str.strip())
    old_code = raw_cohort["FU_einde_diagnose_definitief_ordinaal"].fillna("").astype(str).str.strip()
    diagnosis_text = raw_cohort[DX_TEXT_COLUMNS].fillna("").astype(str).agg(" | ".join, axis=1)

    fallback_allowed = new_codes[DX_NEW_COLUMNS[0]].eq("")
    cholecystolithiasis_label = (
        new_codes.eq("5").any(axis=1)
        | (fallback_allowed & old_code.eq("1"))
        | diagnosis_text.str.contains(CHOLECYSTOLITHIASIS_FINAL_TEXT)
    )
    cholecystitis_label = (
        new_codes.eq("6").any(axis=1)
        | (fallback_allowed & old_code.eq("3"))
        | diagnosis_text.str.contains(CHOLECYSTITIS_FINAL_TEXT)
    )
    other_possible_gallstone_related_label = diagnosis_text.str.contains(
        OTHER_POSSIBLE_GALLSTONE_RELATED_TEXT
    )

    index_date_column = "SEH_datum" if setting == "SEH" else "POLI_datum"
    index_dates = parse_date_series(raw_cohort[index_date_column])

    ultrasound_cholecystitis = pd.Series(False, index=raw_cohort.index)
    for date_column, finding_column in ECHO_CHOLECYSTITIS_EVIDENCE:
        ultrasound_cholecystitis |= evidence_on_or_after_index(
            raw_cohort,
            index_dates,
            date_column,
            raw_cohort[finding_column].eq("1"),
        )

    ct_cholecystitis = pd.Series(False, index=raw_cohort.index)
    for date_column, finding_column in CT_CHOLECYSTITIS_EVIDENCE:
        ct_cholecystitis |= evidence_on_or_after_index(
            raw_cohort,
            index_dates,
            date_column,
            raw_cohort[finding_column].eq("1"),
        )

    cholecystectomy = pd.Series(False, index=raw_cohort.index)
    cholecystolithotomy = pd.Series(False, index=raw_cohort.index)
    cholecystitis_surgery = pd.Series(False, index=raw_cohort.index)
    for performed_column, date_column, type_column, other_column, diagnosis_code_column, diagnosis_text_column in SURGERY_EVIDENCE:
        performed = raw_cohort[performed_column].eq("1")
        other_text = raw_cohort[other_column].fillna("").astype(str)
        diagnosis_text_surgery = raw_cohort[diagnosis_text_column].fillna("").astype(str)
        dated_performed = evidence_on_or_after_index(raw_cohort, index_dates, date_column, performed)
        type_cholecystectomy = raw_cohort[type_column].eq("2") | other_text.str.contains(CHOLE_OPERATION_TEXT)
        type_cholecystolithotomy = raw_cohort[type_column].eq("3")
        cholecystectomy |= dated_performed & type_cholecystectomy
        cholecystolithotomy |= dated_performed & type_cholecystolithotomy
        cholecystitis_surgery |= dated_performed & (
            type_cholecystectomy
            | type_cholecystolithotomy
            | raw_cohort[diagnosis_code_column].eq("6")
            | diagnosis_text_surgery.str.contains(CHOLECYSTITIS_FINAL_TEXT)
        )

    positive_ercp = pd.Series(False, index=raw_cohort.index)
    for date_column, abnormal_column, report_column, diagnosis_column in ERCP_EVIDENCE:
        report_text = raw_cohort[report_column].fillna("").astype(str)
        stone_documented = report_text.str.contains(ERCP_STONE_TEXT)
        ercp_positive = raw_cohort[abnormal_column].eq("1") & raw_cohort[diagnosis_column].eq("1") & stone_documented
        positive_ercp |= evidence_on_or_after_index(
            raw_cohort,
            index_dates,
            date_column,
            ercp_positive,
        )

    cholecystitis_evidence = ultrasound_cholecystitis | ct_cholecystitis | cholecystitis_surgery
    cholecystolithiasis_evidence = cholecystectomy | cholecystolithotomy | positive_ercp
    definition_met_cholecystitis = cholecystitis_label & cholecystitis_evidence
    definition_met_cholecystolithiasis = cholecystolithiasis_label & cholecystolithiasis_evidence
    positive = definition_met_cholecystitis | definition_met_cholecystolithiasis
    labelled_but_unverified = (cholecystitis_label | cholecystolithiasis_label) & ~positive
    other_possible_gallstone_related = (
        other_possible_gallstone_related_label & ~positive & ~labelled_but_unverified
    )

    # Preserve the earlier provisional rule only for a transparent change audit.
    legacy_strict = new_codes.isin(["5", "6"]).any(axis=1)
    legacy_strict |= fallback_allowed & old_code.isin(["1", "3"])
    legacy_explicit_text = diagnosis_text.str.contains(LEGACY_CLEAR_GALLSTONE_TEXT) & ~legacy_strict
    legacy_positive = legacy_strict | legacy_explicit_text

    outcome_known = raw_cohort["FU_einde_diagnose"].isin(["0", "1"])
    negative = outcome_known & ~positive & ~labelled_but_unverified & ~other_possible_gallstone_related
    previous_chole = raw_cohort["VG_cholecystectomie"].isin(["2", "3", "4"])
    previous_chole_no = raw_cohort["VG_cholecystectomie"].eq("1")
    previous_chole_unknown = ~previous_chole & ~previous_chole_no

    status = pd.Series("outcome_missing", index=raw_cohort.index, dtype="string")
    status.loc[positive] = "positive_table_definition_met"
    status.loc[labelled_but_unverified] = "adjudication_required_label_without_table_evidence"
    status.loc[other_possible_gallstone_related] = "adjudication_required_other_gallstone_related_diagnosis"
    status.loc[negative & previous_chole] = "negative_excluded_prior_cholecystectomy"
    status.loc[negative & previous_chole_unknown] = "negative_control_eligibility_unknown"
    status.loc[negative & previous_chole_no] = "negative_control_eligible"

    outcome = pd.Series(pd.NA, index=raw_cohort.index, dtype="Int64")
    outcome.loc[positive] = 1
    outcome.loc[negative & previous_chole_no] = 0
    include = outcome.notna().astype("Int64")

    previous_status = pd.Series("unknown", index=raw_cohort.index, dtype="string")
    previous_status.loc[previous_chole] = "yes_prior_cholecystectomy"
    previous_status.loc[previous_chole_no] = "no_prior_cholecystectomy"

    derived = pd.DataFrame(index=raw_cohort.index)
    derived["analysis_setting"] = setting
    derived["outcome_status"] = status
    derived["gallstone_related_disease"] = outcome
    derived["primary_analysis_include"] = include
    derived["outcome_observed"] = outcome_known.astype("Int64")
    derived["prior_cholecystectomy_status"] = previous_status
    derived["outcome_dx1_def_label"] = new_codes[DX_NEW_COLUMNS[0]].map(DEF_DIAGNOSIS)
    derived["outcome_dx2_def_label"] = new_codes[DX_NEW_COLUMNS[1]].map(DEF_DIAGNOSIS)
    derived["outcome_dx3_def_label"] = new_codes[DX_NEW_COLUMNS[2]].map(DEF_DIAGNOSIS)
    derived["outcome_dx1_old_label"] = old_code.map(OLD_DIAGNOSIS)
    derived["outcome_candidate_cholecystolithiasis"] = cholecystolithiasis_label.astype("Int64")
    derived["outcome_candidate_cholecystitis"] = cholecystitis_label.astype("Int64")
    derived["outcome_candidate_other_gallstone_related"] = other_possible_gallstone_related_label.astype("Int64")
    derived["outcome_evidence_ultrasound_cholecystitis"] = ultrasound_cholecystitis.astype("Int64")
    derived["outcome_evidence_ct_cholecystitis"] = ct_cholecystitis.astype("Int64")
    derived["outcome_evidence_surgery_cholecystitis"] = cholecystitis_surgery.astype("Int64")
    derived["outcome_evidence_laparoscopic_cholecystectomy"] = cholecystectomy.astype("Int64")
    derived["outcome_evidence_laparoscopic_cholecystolithotomy"] = cholecystolithotomy.astype("Int64")
    derived["outcome_evidence_positive_ercp"] = positive_ercp.astype("Int64")
    derived["outcome_definition_met_cholecystolithiasis"] = definition_met_cholecystolithiasis.astype("Int64")
    derived["outcome_definition_met_cholecystitis"] = definition_met_cholecystitis.astype("Int64")
    derived["outcome_previous_provisional_positive"] = legacy_positive.astype("Int64")
    derived["outcome_changed_from_previous_provisional"] = legacy_positive.ne(positive).astype("Int64")
    derived["outcome_positive_source"] = np.select(
        [
            definition_met_cholecystolithiasis & definition_met_cholecystitis,
            definition_met_cholecystolithiasis,
            definition_met_cholecystitis,
        ],
        [
            "both_table_definitions_met",
            "cholecystolithiasis_table_definition_met",
            "cholecystitis_table_definition_met",
        ],
        default="",
    )
    return derived


def add_date_derivations(raw: pd.DataFrame, cohort: pd.DataFrame, index_date_column: str) -> pd.DataFrame:
    index_date = parse_date_series(raw.loc[cohort.index, index_date_column])
    birth = parse_date_series(raw.loc[cohort.index, "Geboortedatum"])
    surgery = parse_date_series(raw.loc[cohort.index, "B_operatie_datum"])
    cohort["index_date_iso"] = index_date.dt.strftime("%Y-%m-%d")
    cohort["age_years_at_index"] = ((index_date - birth).dt.days / 365.2425).round(2)
    cohort["years_since_bariatric_surgery"] = ((index_date - surgery).dt.days / 365.2425).round(2)
    return cohort


def add_numeric_derivations(cohort: pd.DataFrame, setting: str) -> pd.DataFrame:
    vital_columns = {
        "SEH": ["SEH_LO_Temp", "SEH_LO_pols", "SEH_LO_SRR", "SEH_LO_DRR", "SEH_LO_Sat", "SEH_LO_AF", "SEH_LO_gewicht", "SEH_LO_TWL", "SEH_LO_BMI"],
        "POLI": ["POLI_LO_temp", "POLI_LO_pols", "POLI_LO_SRR", "POLI_LO_DRR", "POLI_LO_sat", "POLI_LO_AF", "POLI_LO_gewicht", "POLI_LO_TWL", "POLI_LO_BMI"],
    }[setting]
    additions: dict[str, pd.Series | str] = {}
    for column in vital_columns:
        if column in cohort:
            parsed = cohort[column].map(parse_numeric)
            additions[f"{column}_numeric"] = pd.Series(
                [item[0] for item in parsed], index=cohort.index, dtype="Float64"
            )
    for lab_name, unit in LAB_UNITS.items():
        source = f"index_lab_{lab_name}"
        if source in cohort:
            parsed = cohort[source].map(parse_numeric)
            additions[f"{source}_numeric"] = pd.Series(
                [item[0] for item in parsed], index=cohort.index, dtype="Float64"
            )
            additions[f"{source}_qualifier"] = pd.Series(
                [item[1] for item in parsed], index=cohort.index, dtype="string"
            )
            additions[f"{source}_unit"] = pd.Series(unit, index=cohort.index, dtype="string")
    return pd.concat([cohort, pd.DataFrame(additions, index=cohort.index)], axis=1)


def variable_role(column: str, setting: str) -> tuple[str, str, bool]:
    low = column.casefold()
    if column.startswith("outcome_raw_") or low.startswith("outcome_") or column in {
        "gallstone_related_disease",
        "primary_analysis_include",
        "outcome_observed",
    }:
        return "outcome_or_followup", "not_a_predictor", False
    if column in {
        "Participant Id",
        "Participant Status",
        "Site Abbreviation",
        "Participant Creation Date",
        "episodes",
        "Episode_hoeveel",
        "Episode_nummer_1ste episode_Nummer",
        "Episode_nummer_2de episode_Nummer",
        "Episode_nummer_3de episode_Nummer",
        "Episode_nummer_4de episode_Nummer",
        "analysis_setting",
        "index_date_iso",
        "prior_cholecystectomy_status",
    }:
        return "administrative_or_selection", "not_a_predictor", False
    if "link_status" in low or "source_instance" in low or "same_day_record_count" in low:
        return "audit_linkage", "not_a_predictor", False
    if low.startswith("index_lab_"):
        tier = "primary_SEH_lab" if setting == "SEH" else "secondary_exploratory_POLI_lab"
        return "index_laboratory", tier, True
    if low.startswith(("index_ct_", "index_echo_", "index_mri_", "index_mrcp_")):
        return "index_imaging", "extended_model_only", True
    if low.startswith(("index_urine_", "index_feces_", "index_xboz_", "index_xslik_", "index_gastro_", "index_ercp_", "index_colon_")):
        return "index_other_diagnostics", "review_before_candidate_use", True
    if low.startswith("index_survey_"):
        return "same_day_patient_questionnaire", "review_before_candidate_use", True
    if any(term in low for term in ["_wd", "_dd", "conclusie", "diagnose", "beleid", "startmed", "_fu"]):
        return "initial_assessment_or_care_process", "review_for_leakage_or_subjectivity", True
    if column in {"age_years_at_index", "years_since_bariatric_surgery"}:
        return "derived_index_predictor", "candidate_pool_no_selection", True
    return "first_presentation_information", "candidate_pool_no_selection", True


def availability_state(value: object) -> str:
    if pd.isna(value) or str(value).strip() == "":
        return "blank_or_dependency_not_triggered"
    text = str(value).strip()
    if text.startswith("__") and text.endswith("__"):
        return text.strip("_").casefold()
    return "substantive"


def build_availability(cohort: pd.DataFrame, setting: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    denominator_sets = {
        "all_first_presentations": pd.Series(True, index=cohort.index),
        "primary_analysis_ready": cohort["primary_analysis_include"].eq(1),
    }
    for denominator_name, mask in denominator_sets.items():
        subset = cohort.loc[mask]
        n = len(subset)
        rows: list[dict[str, object]] = []
        for column in cohort.columns:
            states = subset[column].map(availability_state)
            counts = Counter(states)
            source_role, model_tier, candidate = variable_role(column, setting)
            row: dict[str, object] = {
                "setting": setting,
                "denominator_set": denominator_name,
                "variable": column,
                "source_role": source_role,
                "model_tier_preselection": model_tier,
                "candidate_pool_eligible": candidate,
                "denominator_n": n,
                "substantive_n": counts.get("substantive", 0),
                "substantive_pct": round(100 * counts.get("substantive", 0) / n, 1) if n else np.nan,
                "blank_or_dependency_n": counts.get("blank_or_dependency_not_triggered", 0),
                "distinct_substantive_n": int(subset.loc[states.eq("substantive"), column].nunique(dropna=True)),
            }
            for state in sorted(set(states) - {"substantive", "blank_or_dependency_not_triggered"}):
                row[f"state_{state}_n"] = counts.get(state, 0)
            rows.append(row)
        frames.append(pd.DataFrame(rows))
    return pd.concat(frames, ignore_index=True)


def build_categorical_inventory(cohort: pd.DataFrame, setting: str, max_levels: int = 50) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    denominator = len(cohort)
    for column in cohort.columns:
        states = cohort[column].map(availability_state)
        substantive = cohort.loc[states.eq("substantive"), column]
        distinct = int(substantive.nunique(dropna=True))
        if distinct == 0 or distinct > max_levels:
            continue
        source_role, model_tier, candidate = variable_role(column, setting)
        for value, number in substantive.value_counts(dropna=False).items():
            rows.append(
                {
                    "setting": setting,
                    "variable": column,
                    "value": value,
                    "n": int(number),
                    "pct_all_first_presentations": round(100 * int(number) / denominator, 2),
                    "distinct_substantive_n": distinct,
                    "source_role": source_role,
                    "model_tier_preselection": model_tier,
                    "candidate_pool_eligible": candidate,
                }
            )
    return pd.DataFrame(rows)


def build_quality_flags(cohort: pd.DataFrame, setting: str) -> pd.DataFrame:
    issues: list[dict[str, object]] = []
    for column, (low, high) in SANITY_RANGES.items():
        if column not in cohort:
            continue
        numeric = pd.to_numeric(cohort[column], errors="coerce")
        bad = numeric.notna() & ((numeric < low) | (numeric > high))
        for idx, value in numeric.loc[bad].items():
            issues.append(
                {
                    "Participant Id": cohort.at[idx, "Participant Id"],
                    "setting": setting,
                    "variable": column,
                    "value": value,
                    "issue": f"outside_broad_sanity_range_{low}_to_{high}",
                    "action": "review_source; value not automatically changed",
                }
            )
    for column in [name for name in cohort if name.endswith("_numeric")]:
        raw_column = column[: -len("_numeric")]
        if raw_column not in cohort:
            continue
        raw_state = cohort[raw_column].map(availability_state)
        bad = raw_state.eq("substantive") & cohort[column].isna()
        for idx in cohort.index[bad]:
            issues.append(
                {
                    "Participant Id": cohort.at[idx, "Participant Id"],
                    "setting": setting,
                    "variable": raw_column,
                    "value": cohort.at[idx, raw_column],
                    "issue": "substantive_value_not_numeric_parseable",
                    "action": "manual unit/value review; raw value retained",
                }
            )
    return pd.DataFrame(issues)


def flow_rows(cohort: pd.DataFrame, setting: str) -> list[dict[str, object]]:
    status = cohort["outcome_status"]
    return [
        {"setting": setting, "stage": "all_first_presentations", "n": len(cohort)},
        {"setting": setting, "stage": "outcome_observed_flag", "n": int(cohort["outcome_observed"].sum())},
        {"setting": setting, "stage": "gallstone_positive_total", "n": int(cohort["gallstone_related_disease"].eq(1).sum())},
        {"setting": setting, "stage": "positive_table_definition_met", "n": int(status.eq("positive_table_definition_met").sum())},
        {"setting": setting, "stage": "positive_cholecystolithiasis_definition_met", "n": int(cohort["outcome_definition_met_cholecystolithiasis"].eq(1).sum())},
        {"setting": setting, "stage": "positive_cholecystitis_definition_met", "n": int(cohort["outcome_definition_met_cholecystitis"].eq(1).sum())},
        {"setting": setting, "stage": "adjudication_required", "n": int(status.str.startswith("adjudication_required", na=False).sum())},
        {"setting": setting, "stage": "adjudication_label_without_table_evidence", "n": int(status.eq("adjudication_required_label_without_table_evidence").sum())},
        {"setting": setting, "stage": "adjudication_other_gallstone_related_diagnosis", "n": int(status.eq("adjudication_required_other_gallstone_related_diagnosis").sum())},
        {"setting": setting, "stage": "negative_before_control_exclusions", "n": int(status.str.startswith("negative_", na=False).sum())},
        {"setting": setting, "stage": "negative_excluded_prior_cholecystectomy", "n": int(status.eq("negative_excluded_prior_cholecystectomy").sum())},
        {"setting": setting, "stage": "negative_control_eligibility_unknown", "n": int(status.eq("negative_control_eligibility_unknown").sum())},
        {"setting": setting, "stage": "negative_control_eligible", "n": int(status.eq("negative_control_eligible").sum())},
        {"setting": setting, "stage": "outcome_missing_not_imputed", "n": int(status.eq("outcome_missing").sum())},
        {"setting": setting, "stage": "primary_analysis_ready", "n": int(cohort["primary_analysis_include"].eq(1).sum())},
    ]


def build_cohort(
    raw: pd.DataFrame,
    setting: str,
    flag_column: str,
    index_date_column: str,
    setting_start: str,
    setting_end: str,
    modality_flags: dict[str, str],
    groups: dict[str, list[dict[str, object]]],
) -> pd.DataFrame:
    mask = raw["episodes"].eq("1") & raw[flag_column].eq("1")
    common_columns = list(raw.columns[: raw.columns.get_loc("SEH_1stecontact")])
    section_columns = list(raw.columns[raw.columns.get_loc(setting_start) : raw.columns.get_loc(setting_end)])
    raw_cohort = raw.loc[mask].copy()
    cohort = raw_cohort[common_columns + section_columns].copy()
    outcome = derive_outcome(raw_cohort, setting)
    raw_outcome = raw_cohort[OUTCOME_RAW_COLUMNS].rename(columns=lambda c: f"outcome_raw_{c}")
    cohort = pd.concat([outcome, cohort, raw_outcome], axis=1)
    cohort = add_date_derivations(raw, cohort, index_date_column)
    for modality, group in groups.items():
        cohort = add_index_records(
            raw,
            cohort,
            index_date_column,
            group,
            f"index_{modality}",
            modality_flags.get(modality),
        )
    cohort = add_other_diagnostic_blocks(raw, cohort, index_date_column)
    cohort = add_same_day_survey(raw, cohort, index_date_column)
    cohort = clean_dataframe(cohort)
    outcome_integer_columns = [
        "gallstone_related_disease",
        "primary_analysis_include",
        "outcome_observed",
        "outcome_candidate_cholecystolithiasis",
        "outcome_candidate_cholecystitis",
        "outcome_candidate_other_gallstone_related",
        "outcome_evidence_ultrasound_cholecystitis",
        "outcome_evidence_ct_cholecystitis",
        "outcome_evidence_surgery_cholecystitis",
        "outcome_evidence_laparoscopic_cholecystectomy",
        "outcome_evidence_laparoscopic_cholecystolithotomy",
        "outcome_evidence_positive_ercp",
        "outcome_definition_met_cholecystolithiasis",
        "outcome_definition_met_cholecystitis",
        "outcome_previous_provisional_positive",
        "outcome_changed_from_previous_provisional",
    ]
    for column in outcome_integer_columns:
        cohort[column] = pd.to_numeric(cohort[column], errors="coerce").astype("Int64")
    cohort = add_numeric_derivations(cohort, setting)
    front = [
        "Participant Id",
        "analysis_setting",
        "index_date_iso",
        "outcome_status",
        "gallstone_related_disease",
        "primary_analysis_include",
        "outcome_observed",
        "prior_cholecystectomy_status",
        "age_years_at_index",
        "years_since_bariatric_surgery",
        "outcome_dx1_def_label",
        "outcome_dx2_def_label",
        "outcome_dx3_def_label",
        "outcome_dx1_old_label",
        "outcome_positive_source",
    ]
    return cohort[front + [column for column in cohort if column not in front]]


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(SOURCE, sep=";", dtype=str, encoding="utf-8-sig", low_memory=False)
    columns = list(raw.columns)

    lab_starts = [f"Lab{i}" for i in range(1, 12)] + [
        "Lab_basis_6",
        "Lab_basis_7",
        "Lab_basis_7_1",
        "Lab_basis_7_1_1",
        "Lab_basis_7_1_1_1",
        "Lab_basis_7_1_1_1_1",
    ]
    urine_starts = [f"Urine{i}" for i in range(1, 6)]
    xboz_starts = [f"XBOZ{i}" for i in range(1, 6)]
    ct_starts = ["CTabd1", "CT_basis_1", "CT_basis_2", "CT_basis_3", "CT_basis_4", "CT_basis_4_1"]
    echo_starts = [f"Echo{i}" for i in range(1, 6)]
    xslik_starts = [f"Xslik{i}" for i in range(1, 3)]
    gastro_starts = [f"gastro{i}" for i in range(1, 10)]

    groups = {
        "lab": positional_groups(columns, lab_starts, "Urine1", "Lab1_", canonicalizer=canonical_lab_field),
        "urine": positional_groups(columns, urine_starts, "Feces", "Urine1_"),
        "feces": positional_groups(columns, ["Feces"], "XBOZ1", "Feces_"),
        "xboz": positional_groups(columns, xboz_starts, "CTabd1", "XBOZ1_"),
        "ct": positional_groups(columns, ct_starts, "Echo1", "CTabd1_"),
        "echo": positional_groups(columns, echo_starts, "Xslik1", "Echo1_"),
        "xslik": positional_groups(columns, xslik_starts, "gastro1", "Xslik1_"),
        "gastro": positional_groups(columns, gastro_starts, "AO_overig", "gastro1_"),
    }

    seh_flags = {
        "lab": "SEH_AO_lab",
        "urine": "SEH_AO_urine",
        "feces": "SEH_AO_feces",
        "xboz": "SEH_AO_XBOZ",
        "ct": "SEH_AO_CTabd",
        "echo": "SEH_AO_echo",
        "xslik": "SEH_AO_XSLIK",
        "gastro": "SEH_AO_gastro",
    }
    poli_flags = {
        "lab": "POLI_AO_lab",
        "urine": "POLI_AO_urine",
        "feces": "POLI_AO_feces",
        "xboz": "POLI_AO_XBOZ",
        "ct": "POLI_AO_CTabd",
        "echo": "POLI_AO_echo",
        "xslik": "POLI_AO_XSLIK",
        "gastro": "POLI_AO_gastro",
    }

    seh = build_cohort(raw, "SEH", "SEH_1stecontact", "SEH_datum", "SEH_1stecontact", "POLI_1stecontact", seh_flags, groups)
    poli = build_cohort(raw, "POLI", "POLI_1stecontact", "POLI_datum", "POLI_1stecontact", "Lab1", poli_flags, groups)

    seh.to_csv(OUTPUT / "cohort_SEH_first_presentations.csv", index=False, encoding="utf-8-sig")
    poli.to_csv(OUTPUT / "cohort_POLI_first_presentations.csv", index=False, encoding="utf-8-sig")

    availability_seh = build_availability(seh, "SEH")
    availability_poli = build_availability(poli, "POLI")
    availability_seh.to_csv(OUTPUT / "variable_availability_SEH.csv", index=False, encoding="utf-8-sig")
    availability_poli.to_csv(OUTPUT / "variable_availability_POLI.csv", index=False, encoding="utf-8-sig")
    pd.concat(
        [build_categorical_inventory(seh, "SEH"), build_categorical_inventory(poli, "POLI")],
        ignore_index=True,
    ).to_csv(OUTPUT / "categorical_value_inventory.csv", index=False, encoding="utf-8-sig")

    flow = pd.DataFrame(flow_rows(seh, "SEH") + flow_rows(poli, "POLI"))
    flow.to_csv(OUTPUT / "cohort_flow.csv", index=False, encoding="utf-8-sig")

    linkage_rows: list[dict[str, object]] = []
    for setting, frame in [("SEH", seh), ("POLI", poli)]:
        for modality in ["lab", "urine", "feces", "xboz", "ct", "echo", "xslik", "gastro", "mrcp", "ercp", "mri", "colon"]:
            column = f"index_{modality}_link_status"
            if column not in frame:
                continue
            for status, number in frame[column].value_counts(dropna=False).items():
                linkage_rows.append(
                    {"setting": setting, "modality": modality, "link_status": status, "n": int(number)}
                )
    pd.DataFrame(linkage_rows).to_csv(
        OUTPUT / "index_investigation_linkage_summary.csv", index=False, encoding="utf-8-sig"
    )

    quality = pd.concat([build_quality_flags(seh, "SEH"), build_quality_flags(poli, "POLI")], ignore_index=True)
    quality.to_csv(OUTPUT / "value_quality_flags.csv", index=False, encoding="utf-8-sig")

    adjudication_columns = [
        "Participant Id",
        "analysis_setting",
        "index_date_iso",
        "outcome_status",
        "gallstone_related_disease",
        "outcome_candidate_cholecystolithiasis",
        "outcome_candidate_cholecystitis",
        "outcome_candidate_other_gallstone_related",
        "outcome_evidence_ultrasound_cholecystitis",
        "outcome_evidence_ct_cholecystitis",
        "outcome_evidence_surgery_cholecystitis",
        "outcome_evidence_laparoscopic_cholecystectomy",
        "outcome_evidence_laparoscopic_cholecystolithotomy",
        "outcome_evidence_positive_ercp",
        "outcome_previous_provisional_positive",
        "outcome_changed_from_previous_provisional",
        "outcome_raw_FU_einde_diagnose_definitief",
        "outcome_raw_FU_einde_definitieve_diagnose_1_DEF_1",
        "outcome_raw_FU_einde_def_diagnose_overig_1",
        "outcome_raw_FU_einde_diagnose_definitief_ordinaal_2",
        "outcome_raw_FU_einde_def_diagnose_overig_2",
        "outcome_raw_FU_einde_definitieve_diagnose_3_ordinaal",
        "outcome_raw_FU_einde_definitieve_diagnose_3_overig",
    ]
    review = pd.concat([seh, poli], ignore_index=True)
    review = review.loc[
        review["outcome_status"].str.startswith("adjudication_required", na=False),
        adjudication_columns,
    ].copy()
    review["adjudication_note"] = np.where(
        review["outcome_status"].eq("adjudication_required_label_without_table_evidence"),
        "A final cholecystolithiasis/cholecystitis label is present, but the required evidence from the supplied table is not verifiable in the export.",
        "Possible other gallstone-related diagnosis (for example biliary pancreatitis or cholangitis); manual clinical adjudication is required.",
    )
    review.to_csv(OUTPUT / "outcome_adjudication_review.csv", index=False, encoding="utf-8-sig")

    change_review = pd.concat([seh, poli], ignore_index=True)
    change_review = change_review.loc[
        change_review["outcome_changed_from_previous_provisional"].eq(1),
        adjudication_columns,
    ].copy()
    change_review.to_csv(
        OUTPUT / "outcome_definition_change_review.csv", index=False, encoding="utf-8-sig"
    )

    first_episode = raw[raw["episodes"].eq("1")]
    setting_overlap = int((first_episode["SEH_1stecontact"].eq("1") & first_episode["POLI_1stecontact"].eq("1")).sum())
    no_setting = int((~first_episode["SEH_1stecontact"].eq("1") & ~first_episode["POLI_1stecontact"].eq("1")).sum())
    metadata = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_rows": int(len(raw)),
        "source_columns": int(raw.shape[1]),
        "participant_id_unique": bool(raw["Participant Id"].is_unique),
        "first_episode_rows": int(len(first_episode)),
        "first_episode_no_setting": no_setting,
        "first_episode_both_settings": setting_overlap,
        "seh_rows": int(len(seh)),
        "poli_rows": int(len(poli)),
        "castor_missing_mapping": CASTOR_MISSING,
        "outcome_definition": {
            "reference_image": "Presumed and definitive table picture diagnosis.jpeg",
            "matching_final_diagnosis_required": True,
            "cholecystolithiasis": "final diagnosis label plus laparoscopic cholecystectomy, laparoscopic cholecystolithotomy, or stone-positive ERCP on/after index",
            "cholecystitis": "final diagnosis label plus positive ultrasound, CT, or relevant surgical evidence on/after index",
            "ordinal_DEF_candidate_codes": {"5": "cholecystolithiasis", "6": "cholecystitis"},
            "fallback_old_first_diagnosis_codes_if_DEF_missing": {"1": "symptomatic gallstone disease", "3": "cholecystitis"},
            "procedure_without_matching_final_diagnosis": "not_positive",
            "other_possible_gallstone_related_diagnosis": "adjudication_required",
            "label_without_required_table_evidence": "adjudication_required",
            "outcome_imputed": False,
        },
        "control_rule": "Negative records with pre-index cholecystectomy codes 2/3/4 are excluded; unknown history remains ineligible pending review.",
        "index_investigation_rule": "Only records with an investigation date equal to the first-visit date are linked as index information.",
    }
    (OUTPUT / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    def count(frame: pd.DataFrame, status: str) -> int:
        return int(frame["outcome_status"].eq(status).sum())

    def count_adjudication(frame: pd.DataFrame) -> int:
        return int(frame["outcome_status"].str.startswith("adjudication_required", na=False).sum())

    def same_day(frame: pd.DataFrame, modality: str) -> int:
        return int(pd.to_numeric(frame[f"index_{modality}_same_day_record_count"], errors="coerce").fillna(0).gt(0).sum())

    def concordant(frame: pd.DataFrame, modality: str) -> int:
        return int(
            frame[f"index_{modality}_link_status"].isin(["matched", "multiple_same_date_records"]).sum()
        )

    summary = f"""# Stap 1 - galsteen-gerelateerde uitkomst na bariatrische chirurgie

## Cohorten en outcome

| Cohort | Eerste presentaties | Positieve uitkomsten | Adjudicatie nodig | Geschikte negatieve controles | Negatieven uitgesloten: eerdere cholecystectomie | Controle-eligibiliteit onbekend | Outcome ontbreekt | Primair analyseerbaar |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SEH | {len(seh)} | {int(seh['gallstone_related_disease'].eq(1).sum())} | {count_adjudication(seh)} | {count(seh, 'negative_control_eligible')} | {count(seh, 'negative_excluded_prior_cholecystectomy')} | {count(seh, 'negative_control_eligibility_unknown')} | {count(seh, 'outcome_missing')} | {int(seh['primary_analysis_include'].eq(1).sum())} |
| Polikliniek | {len(poli)} | {int(poli['gallstone_related_disease'].eq(1).sum())} | {count_adjudication(poli)} | {count(poli, 'negative_control_eligible')} | {count(poli, 'negative_excluded_prior_cholecystectomy')} | {count(poli, 'negative_control_eligibility_unknown')} | {count(poli, 'outcome_missing')} | {int(poli['primary_analysis_include'].eq(1).sum())} |

De uitkomst volgt nu de definities uit `Presumed and definitive table picture diagnosis.jpeg`. Alle drie definitieve diagnoses aan het einde van follow-up worden beoordeeld, maar een diagnoselabel alleen is niet voldoende:

- **Cholecystolithiasis:** passend definitief diagnoselabel én op of na het indexmoment een laparoscopische cholecystectomie, laparoscopische cholecystolithiotomie of steen-positieve ERCP.
- **Cholecystitis:** passend definitief diagnoselabel én op of na het indexmoment een positieve echo, CT of relevante chirurgische bevestiging.

Een procedure zonder passend definitief diagnoselabel geldt niet als event. Mogelijk andere galsteen-gerelateerde diagnoses, zoals biliaire pancreatitis of cholangitis, en labels zonder het vereiste bewijs blijven outcome-missing tot klinische adjudicatie; ze worden niet als controles gebruikt.

In het SEH-cohort voldoen {int(seh['outcome_definition_met_cholecystolithiasis'].eq(1).sum())} records aan de cholecystolithiasisdefinitie en {int(seh['outcome_definition_met_cholecystitis'].eq(1).sum())} aan de cholecystitisdefinitie; {int((seh['outcome_definition_met_cholecystolithiasis'].eq(1) & seh['outcome_definition_met_cholecystitis'].eq(1)).sum())} voldoen aan beide. In de polikliniek zijn dit respectievelijk {int(poli['outcome_definition_met_cholecystolithiasis'].eq(1).sum())}, {int(poli['outcome_definition_met_cholecystitis'].eq(1).sum())} en {int((poli['outcome_definition_met_cholecystolithiasis'].eq(1) & poli['outcome_definition_met_cholecystitis'].eq(1)).sum())}. Ten opzichte van de eerdere voorlopige codering zijn {int(seh['outcome_changed_from_previous_provisional'].eq(1).sum())} SEH-records en {int(poli['outcome_changed_from_previous_provisional'].eq(1).sum())} poli-records gewijzigd; deze staan in `outcome_definition_change_review.csv`.

## Beschikbaarheid van indexonderzoek

| Cohort | Lab zelfde dag (CRF-concordant) | CT zelfde dag (CRF-concordant) | Echo zelfde dag (CRF-concordant) | Baselinevragenlijst zelfde dag |
|---|---:|---:|---:|---:|
| SEH | {same_day(seh, 'lab')} ({concordant(seh, 'lab')}) / {len(seh)} | {same_day(seh, 'ct')} ({concordant(seh, 'ct')}) / {len(seh)} | {same_day(seh, 'echo')} ({concordant(seh, 'echo')}) / {len(seh)} | {int(seh['index_survey_link_status'].eq('matched').sum())} / {len(seh)} ({100*seh['index_survey_link_status'].eq('matched').mean():.1f}%) |
| Polikliniek | {same_day(poli, 'lab')} ({concordant(poli, 'lab')}) / {len(poli)} | {same_day(poli, 'ct')} ({concordant(poli, 'ct')}) / {len(poli)} | {same_day(poli, 'echo')} ({concordant(poli, 'echo')}) / {len(poli)} | {int(poli['index_survey_link_status'].eq('matched').sum())} / {len(poli)} ({100*poli['index_survey_link_status'].eq('matched').mean():.1f}%) |

Poliklinisch onderzoek dat tijdens het indexbezoek alleen was aangevraagd en op een latere datum is uitgevoerd, is niet als indexpredictor gekoppeld. De waarde vóór haakjes telt ieder onderzoeksrecord op de indexdatum; de waarde tussen haakjes vereist daarnaast overeenstemming met de CRF-vlag voor verricht/aangevraagd onderzoek. De volledige beschikbaarheid per variabele staat in de twee `variable_availability`-bestanden, zowel voor alle eerste presentaties als voor de primair analyseerbare subset.

## Belangrijkste kwaliteitsbevindingen

- Bronbestand: {len(raw)} records en {raw.shape[1]} variabelen; `Participant Id` is uniek.
- Eerste studie-episode: {len(first_episode)} records. Daarvan hebben {no_setting} geen eenduidige eerste setting; overlap tussen SEH en poli: {setting_overlap}.
- Missende outcome wordt niet geïmputeerd.
- `niet verricht`, `niet van toepassing`, `niet gevraagd`, `gevraagd maar onbekend`, meetfout, niet omschreven en lege/afhankelijke velden blijven afzonderlijke toestanden.
- De aangeleverde zip bevat geen afzonderlijke Castor-variabelenlijst met waarde-labels. Exacte tekstlabels en officiële missingcodes zijn geharmoniseerd; overige numerieke categorische codes blijven bewust ongewijzigd en staan met frequenties in `categorical_value_inventory.csv` totdat de officiële variabelenlijst beschikbaar is.
- De export bevat datum/CRF-vlagconflicten voor sommige onderzoeken. Bovendien is binnen-dagvolgorde niet beschikbaar; “zelfde dag” is daarom een bovengrens voor informatie die werkelijk vóór de klinische beoordeling bekend was.
- Waarden buiten brede plausibiliteitsgrenzen zijn alleen gemarkeerd en niet automatisch gewijzigd.

## Nog vereist vóór modellering

1. Klinische adjudicatie van {count_adjudication(seh)} SEH-records en {count_adjudication(poli)} poli-records die niet automatisch aan de tabeldefinitie konden worden toegewezen.
2. Broncontrole van de zes negatieve controles met onbekende cholecystectomie-anamnese (drie per cohort).
3. Gezamenlijke kandidaatpredictorselectie en bepaling van het toegestane aantal parameters op basis van de definitieve eventtelling.
"""
    (OUTPUT / "step1_summary.md").write_text(summary, encoding="utf-8")

    print(flow.to_string(index=False))
    print(f"\nWrote outputs to {OUTPUT}")


if __name__ == "__main__":
    main()
