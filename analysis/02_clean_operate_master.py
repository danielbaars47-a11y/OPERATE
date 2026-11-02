#!/usr/bin/env python3
"""Conservatively clean the complete OPERATE Castor export without imputation.

Run from any directory. All source records and columns remain in original order.
Read every cell as a string with pandas NA inference disabled. The only value
changes are edge-whitespace removal, explicit Castor missing sentinels, and exact
calendar-validated date/datetime normalization in fields supported by both their
name and all substantive values. No clinical type, category, outlier correction,
patient linkage, exclusion, or imputation is inferred.

Empty cleaned cells represent original blanks or explicit sentinels. Their
provenance is available in the variable dictionary and explicit missing ledger.
Numeric codes are not decoded without a source data dictionary. Coded fields
with names containing 'date' or 'tijd' are not automatically dates or times.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import unicodedata
from collections import Counter
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / 'OPERATE data' / 'Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
DEFAULT_OUTPUT = ROOT / 'results' / 'cleaned_2026-09-30'
MISSING = {
    '-95': 'measurement_failed', '-96': 'not_applicable',
    '-97': 'not_asked', '-98': 'asked_unknown', '-99': 'not_done',
}
DATE_MISSING = {f'01-01-29{abs(int(code))}': reason for code, reason in MISSING.items()}
DATE_MISSING.update({f'29{abs(int(code))}-01-01': reason for code, reason in MISSING.items()})
TIME_MISSING = {f'00:00:{abs(int(code)) - 60}': reason for code, reason in MISSING.items()}
DATE_SHAPE = re.compile(r'(?:\d{2}-\d{2}-\d{4}|\d{4}-\d{2}-\d{2})\Z')
DATETIME_SHAPE = re.compile(r'(?:\d{2}-\d{2}-\d{4}|\d{4}-\d{2}-\d{2})[ T]\d{2}:\d{2}:\d{2}\Z')
TIME_SHAPE = re.compile(r'\d{2}:\d{2}(?::\d{2})?\Z')
NUMERIC_SHAPE = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z')
DATE_NAME = re.compile(r'datum|(?:^|[^a-z])date(?:$|[^a-z])|birth', re.I)
TIME_NAME = re.compile(r'(?:^|[^a-z])(?:tijd|time)(?:$|[^a-z])', re.I)
# This name is a coded selection, established from inspection, not a date.
KNOWN_NONDATE_FIELDS = {'GIQLI_control_operatie_date'}
LEDGER_FIELDS = ['source_data_row_1based', 'participant_id', 'variable', 'source_variable',
                 'raw_value', 'trimmed_value', 'sentinel_kind', 'missing_reason']
COUNT_FIELDS = ['variable', 'source_variable', 'reason', 'n', 'fraction_all_rows']


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def snake_names(columns: list[str]) -> list[str]:
    """ASCII names in source order; resolve collisions deterministically."""
    used, result = set(), []
    for original in columns:
        base = unicodedata.normalize('NFKD', original).encode('ascii', 'ignore').decode()
        base = re.sub(r'[^a-z0-9]+', '_', base.lower()).strip('_') or 'unnamed_variable'
        if base[0].isdigit():
            base = 'v_' + base
        candidate, suffix = base, 2
        while candidate in used:
            candidate = f'{base}_{suffix}'
            suffix += 1
        used.add(candidate)
        result.append(candidate)
    return result


def evidence_format(name: str, values: list[str]) -> str:
    """Infer syntactic format, never semantic type from cardinality."""
    substantive = [value for value in values if value and value not in MISSING]
    if not substantive:
        return 'no_substantive_values'
    if name not in KNOWN_NONDATE_FIELDS and DATE_NAME.search(name):
        if all(DATE_SHAPE.fullmatch(value) for value in substantive):
            return 'date'
        if all(DATETIME_SHAPE.fullmatch(value) for value in substantive):
            return 'datetime_without_timezone'
    if TIME_NAME.search(name) and all(TIME_SHAPE.fullmatch(value) for value in substantive):
        return 'time'
    return 'unclassified'


def normalize_date(value: str, field_format: str) -> tuple[str, bool]:
    """Return an ISO value only after exact calendar parsing; preserve failures."""
    if not value or field_format not in {'date', 'datetime_without_timezone'}:
        return value, False
    if field_format == 'date':
        formats = ('%d-%m-%Y', '%Y-%m-%d')
    else:
        formats = ('%d-%m-%Y %H:%M:%S', '%Y-%m-%d %H:%M:%S',
                   '%d-%m-%YT%H:%M:%S', '%Y-%m-%dT%H:%M:%S')
    for fmt in formats:
        try:
            parsed = datetime.strptime(value, fmt)
        except ValueError:
            continue
        if parsed.strftime(fmt) != value:
            continue
        return (parsed.date().isoformat() if field_format == 'date'
                else parsed.isoformat(timespec='seconds')), False
    return value, True


def syntactic_shape(value: str) -> str:
    if NUMERIC_SHAPE.fullmatch(value):
        return 'numeric_looking'
    if DATE_SHAPE.fullmatch(value):
        return 'date_looking'
    if DATETIME_SHAPE.fullmatch(value):
        return 'datetime_looking'
    if TIME_SHAPE.fullmatch(value):
        return 'time_looking'
    return 'other_text'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    if source == output or source.is_relative_to(output):
        raise ValueError('Source must not be inside the output directory.')
    output.mkdir(parents=True, exist_ok=True)
    source_hash_before = sha256(source)
    with source.open(encoding='utf-8-sig', newline='') as handle:
        header = next(csv.reader(handle, delimiter=';'))
    if len(set(header)) != len(header):
        raise ValueError('Duplicate source headers require explicit adjudication before import.')
    raw = pd.read_csv(source, sep=';', encoding='utf-8-sig', dtype=str,
                      keep_default_na=False, na_filter=False, low_memory=False)
    assert raw.columns.tolist() == header, 'Source header changed during import.'
    assert not raw.isna().any().any(), 'Unexpected parser NA values.'
    assert 'Participant Id' in raw.columns, 'Expected export record key missing.'
    names = snake_names(header)
    id_raw = raw['Participant Id'].tolist()
    id_trimmed = raw['Participant Id'].str.strip().tolist()
    assert all(id_trimmed), 'Blank record IDs require review.'
    assert len(set(id_trimmed)) == len(id_trimmed), 'Duplicate record IDs require review.'
    cleaned = raw.copy()
    cleaned.columns = names
    dictionary, ledger, missing_counts = [], [], []
    totals, format_counts = Counter(), Counter()
    protected_tokens = ['NA', 'NO', 'no', 'No', 'unknown', 'onbekend', 'nvt', 'n.v.t.']
    protected_before, protected_after = Counter(), Counter()
    all_changes_authorized = True
    excluded_date_field_preserved = True
    for position, (source_name, clean_name) in enumerate(zip(header, names), start=1):
        original = raw[source_name]
        trimmed = original.str.strip()
        trimmed_list = trimmed.tolist()
        unique_values = trimmed.unique().tolist()
        field_format = evidence_format(source_name, unique_values)
        format_counts[field_format] += 1
        sentinel_map = {value: ('numeric_code', reason) for value, reason in MISSING.items()}
        if field_format == 'date':
            sentinel_map.update({value: ('date_code', reason) for value, reason in DATE_MISSING.items()})
        if field_format == 'time':
            sentinel_map.update({value: ('time_code', reason) for value, reason in TIME_MISSING.items()})
        value_map, invalid_dates, date_changes = {}, set(), set()
        for value in unique_values:
            if value in sentinel_map:
                value_map[value] = ''
            else:
                replacement, invalid = normalize_date(value, field_format)
                value_map[value] = replacement
                if invalid:
                    invalid_dates.add(value)
                if replacement != value:
                    date_changes.add(value)
        result = trimmed.map(value_map)
        cleaned[clean_name] = result
        sentinel_mask = trimmed.isin(sentinel_map)
        raw_blanks = original.eq('')
        whitespace_only = original.ne('') & trimmed.eq('')
        whitespace_changed = original.ne(trimmed)
        conversion_mask = trimmed.isin(date_changes)
        invalid_mask = trimmed.isin(invalid_dates)
        reason_counts = Counter(sentinel_map[value][1] for value in trimmed[sentinel_mask])
        for row in trimmed.index[sentinel_mask]:
            kind, reason = sentinel_map[trimmed.iat[row]]
            ledger.append({
                'source_data_row_1based': int(row) + 1,
                'participant_id': id_trimmed[row], 'variable': clean_name,
                'source_variable': source_name, 'raw_value': original.iat[row],
                'trimmed_value': trimmed.iat[row], 'sentinel_kind': kind,
                'missing_reason': reason,
            })
        column_counts = {'original_blank': int(raw_blanks.sum()),
                         'whitespace_only': int(whitespace_only.sum()), **reason_counts}
        for reason, count in column_counts.items():
            if count:
                missing_counts.append({'variable': clean_name, 'source_variable': source_name,
                                       'reason': reason, 'n': count,
                                       'fraction_all_rows': count / len(raw)})
        nonmissing = result[result.ne('')]
        shape_counts = Counter(syntactic_shape(value) for value in nonmissing)
        nunique = int(nonmissing.nunique())
        record = {
            'source_position_1based': position, 'source_variable': source_name,
            'variable': clean_name, 'storage_type': 'string',
            'syntactic_format': field_format,
            'semantic_type': 'not_inferred_from_export_values',
            'n_rows': len(raw), 'n_original_blank': int(raw_blanks.sum()),
            'n_whitespace_only': int(whitespace_only.sum()),
            'n_whitespace_trimmed': int(whitespace_changed.sum()),
            'n_explicit_sentinel': int(sentinel_mask.sum()),
            'n_clean_missing': int(result.eq('').sum()),
            'missing_fraction': float(result.eq('').mean()),
            'n_clean_observed': len(nonmissing), 'n_distinct_observed': nunique,
            'all_missing_retained': len(nonmissing) == 0,
            'constant_observed_retained': nunique == 1,
            'constant_all_rows_retained': int(result.nunique()) == 1,
            'n_date_values_normalized': int(conversion_mask.sum()),
            'n_invalid_calendar_values_preserved': int(invalid_mask.sum()),
            **{f'n_{shape}': shape_counts[shape] for shape in
               ['numeric_looking', 'date_looking', 'datetime_looking', 'time_looking', 'other_text']},
            **{f'n_{reason}': reason_counts[reason] for reason in MISSING.values()},
            'minimum_observed_character_length': int(nonmissing.str.len().min()) if len(nonmissing) else '',
            'maximum_observed_character_length': int(nonmissing.str.len().max()) if len(nonmissing) else '',
            'interpretation_note': ('Coded selection; date-like name is not a date.'
                                    if source_name in KNOWN_NONDATE_FIELDS else
                                    'Syntactic profiles are not clinical variable definitions.'),
        }
        dictionary.append(record)
        totals.update({key: record[key] for key in ['n_original_blank', 'n_whitespace_only',
                       'n_whitespace_trimmed', 'n_explicit_sentinel', 'n_clean_missing',
                       'n_clean_observed', 'n_date_values_normalized',
                       'n_invalid_calendar_values_preserved']})
        for token in protected_tokens:
            protected_before[token] += int(trimmed.eq(token).sum())
            protected_after[token] += int(result.eq(token).sum())
        unchanged_mask = ~sentinel_mask & ~conversion_mask
        all_changes_authorized &= bool(result[unchanged_mask].equals(trimmed[unchanged_mask]))
        assert int(result.eq('').sum()) == sum(column_counts.values()), source_name
        if source_name in KNOWN_NONDATE_FIELDS:
            excluded_date_field_preserved &= result[~sentinel_mask].equals(trimmed[~sentinel_mask])
    assert all_changes_authorized
    assert excluded_date_field_preserved
    assert protected_before == protected_after, 'Protected literal text was changed.'
    assert cleaned.shape == raw.shape
    assert len(set(names)) == len(header)
    assert cleaned['participant_id'].tolist() == id_trimmed
    master = output / 'operate_cleaned_master.csv.gz'
    cleaned.to_csv(master, index=False, encoding='utf-8', compression={'method': 'gzip', 'mtime': 0})
    pd.DataFrame(dictionary).to_csv(output / 'variable_dictionary.csv', index=False)
    pd.DataFrame(missing_counts, columns=COUNT_FIELDS).to_csv(output / 'missing_reason_counts.csv', index=False)
    pd.DataFrame(ledger, columns=LEDGER_FIELDS).to_csv(
        output / 'missing_reasons_explicit.csv.gz', index=False,
        compression={'method': 'gzip', 'mtime': 0})
    reread = pd.read_csv(master, dtype=str, keep_default_na=False, na_filter=False)
    roundtrip_equal = reread.equals(cleaned)
    assert roundtrip_equal, 'Saved CSV failed exact string roundtrip.'
    source_hash_after = sha256(source)
    assert source_hash_before == source_hash_after, 'Raw source changed during processing.'
    summary = {
        'purpose': 'Conservative full-export master; no imputation, exclusions, recoding, or inferred patient linkage.',
        'source_file': str(source.relative_to(ROOT.resolve())),
        'source_sha256_before': source_hash_before, 'source_sha256_after': source_hash_after,
        'script_file': 'analysis/02_clean_operate_master.py', 'script_sha256': sha256(Path(__file__)),
        'software': {'python': sys.version.split()[0], 'pandas': pd.__version__},
        'n_source_rows': len(raw), 'n_source_columns': len(header),
        'n_clean_rows': len(cleaned), 'n_clean_columns': len(cleaned.columns),
        'n_unique_record_ids': len(set(id_trimmed)),
        'record_unit': 'Export record/episode; Participant Id is not assumed to be a unique person.',
        'n_id_edge_whitespace_changes': sum(a != b for a, b in zip(id_raw, id_trimmed)),
        'n_all_missing_columns_retained': sum(item['all_missing_retained'] for item in dictionary),
        'n_constant_observed_columns_retained': sum(item['constant_observed_retained'] for item in dictionary),
        'n_constant_all_rows_columns_retained': sum(item['constant_all_rows_retained'] for item in dictionary),
        'format_counts': dict(format_counts), 'cell_counts': dict(totals),
        'explicit_sentinel_counts_by_kind': dict(Counter(item['sentinel_kind'] for item in ledger)),
        'explicit_sentinel_counts_by_reason': dict(Counter(item['missing_reason'] for item in ledger)),
        'protected_literal_counts_before': dict(protected_before),
        'protected_literal_counts_after': dict(protected_after),
        'validation': {
            'raw_file_unchanged': source_hash_before == source_hash_after,
            'row_and_column_counts_preserved': cleaned.shape == raw.shape,
            'source_order_preserved': True,
            'record_ids_unique_and_order_preserved': cleaned['participant_id'].tolist() == id_trimmed,
            'cleaned_column_names_unique': len(set(names)) == len(names),
            'unchanged_nondate_nonsentinel_values_after_trim': all_changes_authorized,
            'nondate_giqli_control_field_preserved_after_sentinel_handling': excluded_date_field_preserved,
            'protected_literal_tokens_preserved': protected_before == protected_after,
            'csv_exact_string_roundtrip': roundtrip_equal,
            'missing_reasons_account_for_every_clean_missing_cell':
                sum(item['n'] for item in missing_counts) == totals['n_clean_missing'],
            'explicit_ledger_matches_sentinel_count': len(ledger) == totals['n_explicit_sentinel'],
        },
        'rules': [
            'Read as strings; pandas default missing-token inference disabled.',
            'Trim leading/trailing whitespace only; do not normalize internal spaces or case.',
            'Exact -95/-96/-97/-98/-99 tokens use Castor export missing-code conventions.',
            'Date/time sentinels apply only to fields with date/time name and unanimous compatible value-shape evidence.',
            'Date formats DD-MM-YYYY and YYYY-MM-DD; datetime adds HH:MM:SS, no timezone assumption.',
            'Calendar-invalid dates remain unchanged and are counted for review.',
            'Date ranges, clinical plausibility, code labels and structural applicability are not guessed.',
            'Literal NA, unknown text, and NO/no/No remain observed strings.',
            'No blank/constant column removal, row exclusion, imputation, or imputed-value fabrication.',
            'Read output using keep_default_na=False to preserve legitimate literal NA values.',
        ],
        'outputs': {path.name: {'sha256': sha256(path), 'bytes': path.stat().st_size}
                    for path in [master, output / 'variable_dictionary.csv',
                                 output / 'missing_reason_counts.csv', output / 'missing_reasons_explicit.csv.gz']},
    }
    assert all(summary['validation'].values()), 'One or more validation checks failed.'
    (output / 'validation_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'output': str(output), 'shape': list(cleaned.shape),
                      'cell_counts': dict(totals), 'validation': summary['validation']}, indent=2))


if __name__ == '__main__':
    main()
