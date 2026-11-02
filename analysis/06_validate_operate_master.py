#!/usr/bin/env python3
"""Independently validate the conservative master against the untouched export."""
import hashlib
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/cleaned_2026-09-30"
SOURCE = ROOT / "OPERATE data/Buikpijn_na_bariatrische_chirurg_export_20260622.csv"
raw = pd.read_csv(SOURCE, sep=";", dtype=str, keep_default_na=False, encoding="utf-8-sig")
clean = pd.read_csv(OUT / "operate_cleaned_master.csv.gz", dtype=str, keep_default_na=False)
dictionary = pd.read_csv(OUT / "variable_dictionary.csv", dtype=str, keep_default_na=False)
ledger = pd.read_csv(OUT / "missing_reasons_explicit.csv.gz", dtype=str, keep_default_na=False)
summary = json.loads((OUT / "validation_summary.json").read_text())
mapping = dict(zip(dictionary.source_variable, dictionary.variable))
assert clean.shape == raw.shape
assert raw["Participant Id"].tolist() == clean[mapping["Participant Id"]].tolist()
assert len(mapping) == raw.shape[1] and len(set(mapping.values())) == raw.shape[1]
for source_col in ["GIQLI_control_operatie_date", "POLI_tijd"]:
    a, b = raw[source_col].str.strip(), clean[mapping[source_col]]
    valid = ~a.isin(["-95", "-96", "-97", "-98", "-99"])
    assert a[valid].equals(b[valid]), source_col
for token in ["NA", "NO", "no", "No"]:
    before = sum(int(c.str.strip().eq(token).sum()) for _, c in raw.items())
    after = int(clean.eq(token).sum().sum())
    assert before == after, token
for record in ledger.itertuples(index=False):
    row = int(record.source_data_row_1based) - 1
    assert raw.at[row, record.source_variable] == record.raw_value
    assert clean.at[row, record.variable] == ""
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == summary["source_sha256_before"]
result = {
    "independent_master_validation": "passed",
    "same_dimensions_and_record_order": True,
    "bijection_of_variable_names": True,
    "categorical_date_time_name_exceptions_preserved": True,
    "literal_NA_and_NO_variants_preserved": True,
    "every_explicit_missing_ledger_cell_traced_to_raw_and_clean": True,
    "ledger_cells_checked": len(ledger),
    "raw_source_hash_unchanged": True
}
(OUT / "independent_validation.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
