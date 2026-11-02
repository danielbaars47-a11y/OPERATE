#!/usr/bin/env python3
"""Build complete, source-retaining OPERATE first-visit cohorts before imputation.

Reuses the historical cohort/outcome definitions in 01 without editing that file.
Replaces unsafe scalar cleaning/date inference and adds explicit timing gates,
source/analysis missing-status matrices, typed numeric companions, and audits.
No predictor selection, missingness-threshold exclusion, imputation or modelling.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from datetime import datetime
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'OPERATE data' / 'Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
OUTPUT = ROOT / 'results' / 'first_visit_2026-09-30'

def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'analysis' / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

old = load_module('operate_existing_step1', '01_prepare_gallstone_step1.py')
master = load_module('operate_master_cleaner', '02_clean_operate_master.py')
TOKENS = {f'__{reason.upper()}__': reason for reason in master.MISSING.values()}
TOKENS['__PERFORMED_DETAILS_UNAVAILABLE__'] = 'performed_details_unavailable'
TEXT_MISSING_STATES = {
    'unknown': 'asked_unknown', 'onbekend': 'asked_unknown',
    'niet omschreven': 'not_described', 'niet beschreven': 'not_described',
    '-': 'not_described', '--': 'not_described',
    'nvt': 'not_applicable', 'n.v.t.': 'not_applicable',
    'niet van toepassing': 'not_applicable',
    'niet verricht': 'not_done', 'niet gedaan': 'not_done',
    'niet gevraagd': 'not_asked', 'geen antwoord ingevuld': 'not_recorded',
}
NUM = re.compile(r'(<=|>=|<|>)?\s*([+-]?(?:\d+(?:[.,]\d*)?|[.,]\d+)(?:[eE][+-]?\d+)?)\Z')
BASE_NUM = {'Startgewicht': ('kg', (20, 400)), 'Lengte': ('cm', (100, 250)),
            'Start_BMI': ('kg/m2', (10, 100)), 'B_gewicht_OKdag': ('not_specified_in_crf_field', (20, 400))}
VITAL_NUM = {'temp': ('degree_C', (30, 45)), 'pols': ('beats/min', (20, 250)),
             'srr': ('mmHg', (50, 300)), 'drr': ('mmHg', (20, 200)),
             'sat': ('percent', (50, 100)), 'af': ('breaths/min', (4, 80)),
             'gewicht': ('not_specified_in_crf_field', (20, 400)),
             'twl': ('percent', (-100, 100)), 'bmi': ('kg/m2', (10, 100))}
LAB_UNITS = dict(old.LAB_UNITS)
LAB_UNITS['EGFR'] = 'unresolved_CRF_prints_kl/1.73m2'
LAB_RANGE = {name: old.SANITY_RANGES[f'index_lab_{name}_numeric'] for name in LAB_UNITS}
EXPLICIT_OUTCOME = {'gallstone_related_disease', 'primary_analysis_include', 'outcome_observed'}


def empty(value):
    return pd.isna(value) or str(value).strip() == ''


def safe_text(value):
    if empty(value):
        return pd.NA
    text = str(value).strip()
    return old.CASTOR_MISSING.get(text, text)


@lru_cache(maxsize=40000)
def exact_date(text):
    if not text or text in master.DATE_MISSING or text in TOKENS or text in master.MISSING:
        return None
    if not (re.fullmatch(r'\d{2}-\d{2}-\d{4}', text) or re.fullmatch(r'\d{4}-\d{2}-\d{2}', text)):
        return None
    for fmt in ('%d-%m-%Y', '%Y-%m-%d'):
        try:
            date = datetime.strptime(text, fmt)
        except ValueError:
            continue
        if date.strftime(fmt) == text:
            return date
    return None


def date_series(series):
    values = [None if empty(v) else exact_date(str(v).strip()) for v in series]
    return pd.Series(pd.to_datetime(values, errors='coerce'), index=series.index)


def dates_equal(left, right):
    if empty(left) or empty(right):
        return False
    a, b = exact_date(str(left).strip()), exact_date(str(right).strip())
    return a is not None and b is not None and a == b


def safe_clean(frame):
    result = {}
    for column in frame:
        s = frame[column].map(safe_text)
        vals = [str(v) for v in s.dropna().unique()]
        fmt = master.evidence_format(column, [v for v in vals if v not in TOKENS])
        lookup = {}
        for v in vals:
            if fmt == 'date' and v in master.DATE_MISSING:
                lookup[v] = f'__{master.DATE_MISSING[v].upper()}__'
            elif fmt == 'time' and v in master.TIME_MISSING:
                lookup[v] = f'__{master.TIME_MISSING[v].upper()}__'
            else:
                lookup[v] = master.normalize_date(v, fmt)[0]
        result[column] = s.map(lambda v: pd.NA if pd.isna(v) else lookup[str(v)])
    return pd.DataFrame(result, index=frame.index)


# Patch only this imported in-memory module; the original file is never changed.
old.clean_exact_text = safe_text
old.clean_dataframe = safe_clean
old.parse_date_series = date_series
old.date_equal = dates_equal
old.add_numeric_derivations = lambda frame, setting: frame


def strict_index_records(raw, cohort, index_date_column, groups, output_prefix,
                         performed_raw_column, performed_yes_value='2'):
    """Same retained linkage rule as 01; all date validation is now exact."""
    rows = []
    canonical_fields = sorted({canon for group in groups for canon, _ in group['mapping']})
    for row_id in cohort.index:
        source_row = raw.loc[row_id]
        index_date = source_row.get(index_date_column)
        matches = [group for group in groups if dates_equal(source_row.get(group['date']), index_date)]
        performed_value = source_row.get(performed_raw_column) if performed_raw_column else None
        performed = not empty(performed_value) and str(performed_value).strip() == performed_yes_value
        no = bool(performed_raw_column and not empty(performed_value) and str(performed_value).strip() == '1')
        index_missing = empty(index_date) or exact_date(str(index_date).strip()) is None
        if index_missing:
            status = 'index_date_missing'
        elif matches and performed:
            status = 'matched' if len(matches) == 1 else 'multiple_same_date_records'
        elif matches and no:
            status = 'record_match_flag_no_conflict'
        elif matches:
            status = 'record_match_performance_unknown'
        elif performed:
            status = 'performed_or_requested_no_same_day_detail'
        elif no:
            status = 'not_performed_or_not_requested'
        else:
            status = 'performance_unknown'
        out = {f'{output_prefix}_link_status': status, f'{output_prefix}_same_day_record_count': len(matches)}
        if not matches and performed:
            for field in canonical_fields:
                if field not in {'record_present', 'date', 'time'}:
                    out[f'{output_prefix}_{field}'] = '__PERFORMED_DETAILS_UNAVAILABLE__'
        for number, group in enumerate(matches, start=1):
            prefix = output_prefix if number == 1 else f'{output_prefix}_repeat{number}'
            out[f'{prefix}_source_instance'] = group['selector']
            for canonical, source_name in group['mapping']:
                out[f'{prefix}_{canonical}'] = source_row.get(source_name)
        rows.append(out)
    return pd.concat([cohort, pd.DataFrame(rows, index=cohort.index)], axis=1)


old.add_index_records = strict_index_records


def state(value):
    if empty(value):
        return 'blank_or_dependency_not_triggered'
    return TOKENS.get(str(value), TEXT_MISSING_STATES.get(str(value).strip().casefold(), 'present_code_or_text'))


def make_states(frame):
    return frame.apply(lambda s: s.map(state))


def missing_to_empty(frame):
    return frame.apply(lambda s: s.map(lambda v: '' if empty(v) or str(v) in TOKENS else v))


def outcome_column(column):
    return column.startswith('outcome_') or column in EXPLICIT_OUTCOME


def role(column, setting):
    low = column.lower()
    if outcome_column(column):
        return 'outcome_or_reference_standard_separate_file'
    if column in {'Participant Id', 'Participant Status', 'Participant Creation Date', 'Site Abbreviation',
                  'episodes', 'Episode_hoeveel', 'analysis_setting', 'index_date_iso',
                  'prior_cholecystectomy_status', 'duplicate_review_required', 'duplicate_review_group'} or column.startswith('Episode_nummer_'):
        return 'administrative_or_selection'
    if any(token in low for token in ['_source_instance', '_link_status', '_same_day_record_count']):
        return 'linkage_audit'
    if final_diagnosis(column) or retrospective_diagnostic(column):
        return 'reference_standard_timing_review'
    if column.startswith(('GIQLI', 'Andere_studies', 'Reguliere_FU')):
        return 'undated_context_timing_unconfirmed'
    if column.startswith(('RFU_', 'Date_RFU_', 'verleden_')) or history30(column):
        return 'historical_information_with_timing_gate'
    if low.startswith('index_lab_'):
        if any(term in low for term in ['diagnose', 'bevinding', '_wd', '_dd', 'conclusie', 'beleid']):
            return 'assessment_or_care_process_review'
        return 'same_day_laboratory'
    if low.startswith('index_survey_'):
        return 'same_day_questionnaire'
    if low.startswith('index_'):
        if any(term in low for term in ['diagnose', 'bevinding', '_wd', '_dd', 'conclusie', 'beleid']):
            return 'assessment_or_care_process_review'
        return 'same_day_diagnostic_information'
    if any(term in low for term in ['ingevuld_door', 'jaren_ervaring']):
        return 'administrative_assessor_metadata'
    if re.search(r'^(seh|poli)_(?:m_|fu|dd|wd|injectie|verwijzing)', low) or any(
            term in low for term in ['_wd', '_dd', 'conclusie', 'beleid', 'startmed', 'overlegd', '_fu']):
        return 'assessment_or_care_process_review'
    if column in {'age_years_at_index', 'years_since_bariatric_surgery'}:
        return 'derived_index_measurement'
    return 'baseline_or_first_visit_source_field'


def final_diagnosis(column):
    return bool(re.match(r'^(?:SEH_+diagnose|Diagnose_+SEH_+definitief|Poli_+diagnose|POLI_+diagnose_+DEF)', column, re.I))


def retrospective_diagnostic(column):
    # CRF diagnosis-contribution questions and CT IH truth are retrospective.
    return bool(re.fullmatch(r'index_(?:lab|urine|feces|xboz|ct|echo|xslik|gastro|colon|mrcp|ercp|mri)(?:_repeat\d+)?_diagnose', column, re.I)
                or re.fullmatch(r'index_ct(?:_repeat\d+)?_IH_JN', column, re.I))


def history30(column):
    return column.startswith(('B_operatie_complicaties30dgn', 'B_compl_30dgn', 'B_AO_30dgn',
                              'B_SEH_30dgn', 'B_heropname_30dgn', 'B_heroperatie_30dgn'))


def build_groups(columns):
    labs = [f'Lab{i}' for i in range(1, 12)] + ['Lab_basis_6', 'Lab_basis_7', 'Lab_basis_7_1',
           'Lab_basis_7_1_1', 'Lab_basis_7_1_1_1', 'Lab_basis_7_1_1_1_1']
    specs = {
        'lab': (labs, 'Urine1', 'Lab1_', old.canonical_lab_field),
        'urine': ([f'Urine{i}' for i in range(1, 6)], 'Feces', 'Urine1_', None),
        'feces': (['Feces'], 'XBOZ1', 'Feces_', None),
        'xboz': ([f'XBOZ{i}' for i in range(1, 6)], 'CTabd1', 'XBOZ1_', None),
        'ct': (['CTabd1', 'CT_basis_1', 'CT_basis_2', 'CT_basis_3', 'CT_basis_4', 'CT_basis_4_1'], 'Echo1', 'CTabd1_', None),
        'echo': ([f'Echo{i}' for i in range(1, 6)], 'Xslik1', 'Echo1_', None),
        'xslik': (['Xslik1', 'Xslik2'], 'gastro1', 'Xslik1_', None),
        'gastro': ([f'gastro{i}' for i in range(1, 10)], 'AO_overig', 'gastro1_', None),
    }
    return {key: old.positional_groups(columns, starts, end, prefix, canonicalizer=canon)
            for key, (starts, end, prefix, canon) in specs.items()}


def timing_blocks(columns):
    blocks = []
    dates = [c for c in columns if c.startswith('Date_RFU_')]
    for i, date_col in enumerate(dates):
        start = columns.index(date_col)
        end = columns.index(dates[i + 1]) if i + 1 < len(dates) else columns.index('RFU_FEtest')
        blocks.append((date_col, date_col, columns[start:end]))
    blocks.append(('RFU_FEtest', 'RFU_FEtest_date', [c for c in columns if c.startswith('RFU_FEtest')]))
    history_groups = {}
    for c in columns:
        match = re.match(r'^(verleden_(?:SEH|heropname|heroperatie)_\d+)_', c)
        if match:
            history_groups.setdefault(match.group(1), []).append(c)
    for group, cols in history_groups.items():
        date_candidates = [c for c in cols if 'datum' in c.lower()]
        assert len(date_candidates) == 1, (group, date_candidates)
        blocks.append((group, date_candidates[0], cols))
    return blocks


def timing_gate(snapshot, source_states, setting):
    clean, states = snapshot.copy(), source_states.copy()
    # Exact unambiguous textual missing values are harmonized only in the analysis
    # copy; the before-timing snapshot retains the original cleaned source text.
    clean = clean.mask(source_states.isin(set(TEXT_MISSING_STATES.values())))
    index = date_series(snapshot['index_date_iso'])
    audit, issues = [], []
    for block, date_col, fields in timing_blocks(list(snapshot.columns)):
        dates = date_series(snapshot[date_col])
        status = pd.Series('timing_date_missing_or_invalid', index=snapshot.index, dtype='string')
        status.loc[index.isna()] = 'index_date_missing_or_invalid'
        status.loc[index.notna() & dates.notna() & dates.lt(index)] = 'before_index_retained'
        status.loc[index.notna() & dates.notna() & dates.eq(index)] = 'same_day_timing_unconfirmed'
        status.loc[index.notna() & dates.notna() & dates.gt(index)] = 'after_index_withheld'
        hide = status.ne('before_index_retained')
        clean.loc[hide, fields] = pd.NA
        for col in fields:
            states.loc[hide, col] = status.loc[hide]
        for ridx in snapshot.index:
            audit.append({'Participant Id': snapshot.at[ridx, 'Participant Id'], 'setting': setting,
                          'block': block, 'date_variable': date_col, 'status': status.at[ridx],
                          'n_fields': len(fields), 'n_source_values_present': int(
                              source_states.loc[ridx, fields].eq('present_code_or_text').sum())})
    surgery = date_series(snapshot['B_operatie_datum'])
    days = (index - surgery).dt.days
    status = pd.Series('surgery_or_index_date_missing_or_invalid', index=snapshot.index, dtype='string')
    status.loc[days.notna() & days.lt(0)] = 'surgery_after_index_withheld'
    status.loc[days.notna() & days.ge(0) & days.lt(30)] = 'postoperative_30day_window_incomplete'
    status.loc[days.notna() & days.ge(30)] = 'postoperative_30day_window_complete'
    fields = [c for c in snapshot if history30(c)]
    hide = status.ne('postoperative_30day_window_complete')
    clean.loc[hide, fields] = pd.NA
    for col in fields:
        states.loc[hide, col] = status.loc[hide]
    for ridx in snapshot.index:
        audit.append({'Participant Id': snapshot.at[ridx, 'Participant Id'], 'setting': setting,
                      'block': 'postoperative_30day_history', 'date_variable': 'B_operatie_datum',
                      'status': status.at[ridx], 'n_fields': len(fields),
                      'n_source_values_present': int(source_states.loc[ridx, fields].eq('present_code_or_text').sum())})
        if pd.notna(days.at[ridx]) and days.at[ridx] < 0:
            issues.append({'Participant Id': snapshot.at[ridx, 'Participant Id'], 'setting': setting,
                           'variable': 'B_operatie_datum', 'issue': 'surgery_after_index_review_required',
                           'value': str(snapshot.at[ridx, 'B_operatie_datum'])})
    for predicate, reason in [(lambda c: c.startswith(('GIQLI', 'Andere_studies', 'Reguliere_FU')), 'timing_unconfirmed'),
                               (lambda c: final_diagnosis(c) or retrospective_diagnostic(c), 'reference_standard_timing_review')]:
        fields = [c for c in snapshot if predicate(c)]
        clean.loc[:, fields] = pd.NA
        states.loc[:, fields] = reason
        for ridx in snapshot.index:
            audit.append({'Participant Id': snapshot.at[ridx, 'Participant Id'], 'setting': setting,
                          'block': reason, 'date_variable': '', 'status': reason, 'n_fields': len(fields),
                          'n_source_values_present': int(source_states.loc[ridx, fields].eq('present_code_or_text').sum())})
    return clean, states, audit, issues


def add_linkage_missing_context(frame, states):
    """Explain unavailable detail cells without overwriting recorded findings."""
    contextual = {'not_performed_or_not_requested', 'index_date_missing', 'performance_unknown',
                  'no_same_day_record', 'no_same_day_survey'}
    replaceable = {'blank_or_dependency_not_triggered', 'performed_details_unavailable'}
    for link_col in [c for c in frame if c.endswith('_link_status')]:
        prefix = link_col[:-len('_link_status')] + '_'
        detail_cols = [c for c in frame if c.startswith(prefix) and not any(
            c.endswith(suffix) for suffix in ['_link_status', '_source_instance', '_same_day_record_count'])]
        link_states = frame[link_col].astype('string')
        for col in detail_cols:
            use = link_states.isin(contextual) & states[col].isin(replaceable)
            states.loc[use, col] = link_states.loc[use]
    return states


def numeric_specs(columns):
    specs = {}
    for col in columns:
        if col in BASE_NUM:
            specs[col] = BASE_NUM[col]
        match = re.fullmatch(r'(SEH|POLI)_LO_(Temp|temp|pols|SRR|DRR|Sat|sat|AF|gewicht|TWL|BMI)', col)
        if match:
            specs[col] = VITAL_NUM[match.group(2).lower()]
        match = re.fullmatch(r'index_lab(?:_repeat\d+)?_(.+)', col)
        if match and match.group(1) in LAB_UNITS:
            name = match.group(1)
            specs[col] = (LAB_UNITS[name], LAB_RANGE[name])
    return specs


def add_numeric(frame, states, setting):
    added, added_states, metadata, flags = {}, {}, {}, []
    for col, (unit, bounds) in numeric_specs(frame.columns).items():
        nums, thresholds, qualifiers, reasons = [], [], [], []
        for ridx, value in frame[col].items():
            original_reason = states.at[ridx, col]
            if original_reason != 'present_code_or_text' or empty(value) or str(value) in TOKENS:
                nums.append(pd.NA); thresholds.append(pd.NA); qualifiers.append(pd.NA)
                reasons.append(original_reason)
                continue
            match = NUM.fullmatch(str(value).strip())
            if not match:
                nums.append(pd.NA); thresholds.append(pd.NA); qualifiers.append(pd.NA)
                reasons.append('numeric_parse_unresolved')
                flags.append({'Participant Id': frame.at[ridx, 'Participant Id'], 'setting': setting,
                              'variable': col, 'issue': 'numeric_parse_unresolved', 'value': str(value)})
                continue
            number, qualifier = float(match.group(2).replace(',', '.')), match.group(1) or '='
            thresholds.append(number); qualifiers.append(qualifier)
            if not bounds[0] <= number <= bounds[1]:
                nums.append(pd.NA); reasons.append('numeric_plausibility_review')
                flags.append({'Participant Id': frame.at[ridx, 'Participant Id'], 'setting': setting,
                              'variable': col, 'issue': 'numeric_plausibility_review', 'value': str(value),
                              'screen_min': bounds[0], 'screen_max': bounds[1]})
            elif qualifier != '=':
                nums.append(pd.NA); reasons.append('censored_measurement_bound_retained')
            else:
                nums.append(number); reasons.append('numeric_observed')
        for suffix, values, dtype in [('numeric', nums, 'Float64'), ('bound', thresholds, 'Float64'),
                                      ('qualifier', qualifiers, 'string')]:
            name = f'{col}_{suffix}'
            added[name] = pd.Series(values, index=frame.index, dtype=dtype)
            if suffix == 'numeric':
                cell_states = reasons
            else:
                cell_states = [('numeric_bound_or_qualifier_recorded' if not empty(v) else r)
                               for v, r in zip(values, reasons)]
            added_states[name] = pd.Series(cell_states, index=frame.index, dtype='string')
            metadata[name] = {'source_column': col, 'numeric_companion': suffix, 'unit': unit,
                              'screen_min': bounds[0], 'screen_max': bounds[1]}
    for col, bounds in [('age_years_at_index', (18, 100)), ('years_since_bariatric_surgery', (0, 40))]:
        numbers = pd.to_numeric(frame[col], errors='coerce')
        bad = numbers.notna() & ~numbers.between(*bounds)
        for ridx in frame.index[bad]:
            flags.append({'Participant Id': frame.at[ridx, 'Participant Id'], 'setting': setting,
                          'variable': col, 'issue': 'numeric_plausibility_review', 'value': str(frame.at[ridx, col]),
                          'screen_min': bounds[0], 'screen_max': bounds[1]})
        frame.loc[bad, col] = pd.NA
        states.loc[bad, col] = 'numeric_plausibility_review'
    return pd.concat([frame, pd.DataFrame(added, index=frame.index)], axis=1), pd.concat(
        [states, pd.DataFrame(added_states, index=frame.index)], axis=1), metadata, flags


def write_csv(frame, path):
    compression = {'method': 'gzip', 'mtime': 0} if path.name.endswith('.gz') else None
    missing_to_empty(frame).to_csv(path, index=False, compression=compression)


def availability(frame, states, source_states, subset_mask, subset_name, setting, names):
    rows = []
    for col in frame:
        counts = states.loc[subset_mask, col].value_counts()
        observed = int(counts.get('present_code_or_text', 0) + counts.get('numeric_observed', 0)
                       + counts.get('numeric_bound_or_qualifier_recorded', 0))
        n = int(subset_mask.sum())
        row = {'setting': setting, 'subset': subset_name, 'variable': names[col],
               'original_or_derived_name': col, 'role': role(col, setting), 'n_total': n,
               'n_value_present': observed, 'n_missing_or_withheld': n - observed,
               'missing_or_withheld_percent': round(100 * (n - observed) / n, 3) if n else None,
               'n_unique_nonmissing': int(frame.loc[subset_mask & states[col].isin(['present_code_or_text', 'numeric_observed', 'numeric_bound_or_qualifier_recorded']), col].dropna().astype(str).nunique()),
               'interpretation': 'Present code/text does not establish clinical meaning or applicability.'}
        row.update({f'n_{key}': int(count) for key, count in counts.items()})
        if col in source_states:
            row['n_source_present_before_timing_gate'] = int(source_states.loc[subset_mask, col].eq('present_code_or_text').sum())
        rows.append(row)
    return pd.DataFrame(rows)


def source_options(column, groups):
    if column.startswith('outcome_raw_'):
        return [column[len('outcome_raw_'):]]
    if column.startswith('index_survey_'):
        return [column[len('index_survey_'):]]
    match = re.fullmatch(r'index_([^_]+)(?:_repeat\d+)?_(.+)', column)
    if match and match.group(1) in groups:
        return sorted({source for group in groups[match.group(1)] for canonical, source in group['mapping']
                       if canonical == match.group(2)})
    if match and match.group(1) in {'colon', 'mrcp', 'ercp', 'mri'}:
        modality, field = match.groups()
        label = {'colon':'colon', 'mrcp':'MRCP', 'ercp':'ERCP', 'mri':'MRI'}[modality]
        source_field = 'datum' if field == 'date' else field
        if source_field in {'datum','afwijkend','verslag','diagnose'}:
            return [f'AO_overig_{label}_{source_field}{suffix}' for suffix in ['', '_1', '_2', '_3']]
    return []


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    source_hash = master.sha256(SOURCE)
    old_hash = master.sha256(ROOT / 'analysis' / '01_prepare_gallstone_step1.py')
    raw = pd.read_csv(SOURCE, sep=';', encoding='utf-8-sig', dtype=str, keep_default_na=False, na_filter=False)
    raw = safe_clean(raw)
    assert raw.shape == (2536, 6221)
    assert raw['Participant Id'].notna().all() and raw['Participant Id'].is_unique
    columns = raw.columns.tolist()
    groups = build_groups(columns)
    first = raw['episodes'].eq('1')
    seh_first = first & raw['SEH_1stecontact'].eq('1')
    poli_first = first & raw['POLI_1stecontact'].eq('1')
    unassigned = first & ~raw['SEH_1stecontact'].eq('1') & ~raw['POLI_1stecontact'].eq('1')
    both = seh_first & poli_first
    assert int(first.sum()) == 1593
    assert int(unassigned.sum()) == 28
    assert not both.any(), 'Setting overlap requires adjudication.'
    unassigned_export = raw.loc[unassigned].copy()
    unassigned_export.columns = master.snake_names(columns)
    write_csv(unassigned_export, OUTPUT / 'first_episodes_unassigned_setting.csv.gz')
    flow = [{'stage': 'all_source_records', 'n': len(raw)}, {'stage': 'first_episode_records', 'n': int(first.sum())},
            {'stage': 'not_marked_first_episode_records', 'n': int((~first).sum())},
            {'stage': 'first_episode_SEH', 'n': int(seh_first.sum())}, {'stage': 'first_episode_POLI', 'n': int(poli_first.sum())},
            {'stage': 'first_episode_unassigned_setting', 'n': int(unassigned.sum())},
            {'stage': 'first_episode_both_settings', 'n': int(both.sum())}]
    pd.DataFrame(flow).to_csv(OUTPUT / 'cohort_flow.csv', index=False)
    duplicate_groups = {}
    duplicate_path = ROOT / 'results' / 'cleaned_2026-09-30' / 'candidate_duplicate_first_episodes.csv'
    if duplicate_path.exists():
        pairs = pd.read_csv(duplicate_path, dtype=str, keep_default_na=False)
        for number, (_, pair) in enumerate(pairs.iterrows(), start=1):
            group = f'possible_duplicate_{number}'
            for key in ['participant_id_a', 'participant_id_b']:
                duplicate_groups[str(pair[key])] = group
    unassigned_export['duplicate_review_required'] = unassigned_export['participant_id'].isin(duplicate_groups).astype(int)
    unassigned_export['duplicate_review_group'] = unassigned_export['participant_id'].map(duplicate_groups).fillna('')
    write_csv(unassigned_export, OUTPUT / 'first_episodes_unassigned_setting.csv.gz')
    summaries, all_audits, all_flags, inventories = {}, [], [], []
    for setting, stop, mask in [('SEH', 'POLI_1stecontact', seh_first), ('POLI', 'Lab1', poli_first)]:
        flags = {key: f'{setting}_AO_{suffix}' for key, suffix in {'lab':'lab', 'urine':'urine', 'feces':'feces',
                 'xboz':'XBOZ', 'ct':'CTabd', 'echo':'echo', 'xslik':'XSLIK', 'gastro':'gastro'}.items()}
        cohort = old.build_cohort(raw, setting, f'{setting}_1stecontact', f'{setting}_datum',
                                 f'{setting}_1stecontact', stop, flags, groups)
        # Retain every defined diagnostic field even when this cohort has no linked value.
        absent = {}
        for modality, these_groups in groups.items():
            for field in sorted({field for group in these_groups for field, _ in group['mapping']}):
                name = f'index_{modality}_{field}'
                if name not in cohort:
                    absent[name] = pd.NA
        cohort = pd.concat([cohort, pd.DataFrame(absent, index=cohort.index)], axis=1)
        cohort['duplicate_review_required'] = cohort['Participant Id'].isin(duplicate_groups).astype(int)
        cohort['duplicate_review_group'] = cohort['Participant Id'].map(duplicate_groups).fillna('')
        cohort['outcome_no_definitive_diagnosis'] = raw.loc[cohort.index, 'FU_einde_diagnose'].eq('0').astype(int)
        cohort['outcome_definition_status'] = 'provisional_existing_rule_not_analysis_freeze'
        cohort['outcome_observed_interpretation'] = 'definitive_diagnosis_yes_no_flag_answered_not_verified_outcome'
        snapshot = cohort.copy()
        source_states = make_states(snapshot)
        clean, analysis_states, audits, issues = timing_gate(snapshot, source_states, setting)
        analysis_states = add_linkage_missing_context(clean, analysis_states)
        clean, analysis_states, numeric_meta, numeric_flags = add_numeric(clean, analysis_states, setting)
        all_audits.extend(audits); all_flags.extend(issues + numeric_flags)
        all_columns = list(snapshot.columns) + [c for c in clean if c not in snapshot]
        renamed = master.snake_names(all_columns)
        names = dict(zip(all_columns, renamed))
        outcols = [c for c in snapshot if outcome_column(c)]
        outcomes = snapshot[['Participant Id', 'analysis_setting'] + outcols].copy()
        primary = pd.to_numeric(snapshot['primary_analysis_include'], errors='coerce').eq(1)
        events = pd.to_numeric(snapshot['gallstone_related_disease'], errors='coerce').eq(1)
        analysis_cols = [c for c in clean if not outcome_column(c)]
        cleaned_analysis = clean[analysis_cols]
        cleaned_states = analysis_states[analysis_cols]
        # Status matrices have ID/setting values as keys; other cells contain state strings.
        source_state_export = source_states.copy()
        analysis_state_export = cleaned_states.copy()
        for key in ['Participant Id', 'analysis_setting']:
            source_state_export[key] = snapshot[key]
            analysis_state_export[key] = clean[key]
        write_csv(snapshot.rename(columns=names), OUTPUT / f'cohort_{setting}_before_timing_filter.csv.gz')
        write_csv(cleaned_analysis.rename(columns=names), OUTPUT / f'cohort_{setting}_cleaned.csv.gz')
        write_csv(source_state_export.rename(columns=names), OUTPUT / f'cohort_{setting}_source_missing_status.csv.gz')
        write_csv(analysis_state_export.rename(columns=names), OUTPUT / f'cohort_{setting}_analysis_missing_status.csv.gz')
        write_csv(outcomes.rename(columns=names), OUTPUT / f'outcomes_{setting}_existing_definition.csv')
        dictionary = []
        for col in all_columns:
            meta = numeric_meta.get(col, {})
            parent = meta.get('source_column', col)
            sources = [parent] if parent in raw else source_options(parent, groups)
            dictionary.append({'setting': setting, 'variable': names[col], 'original_or_derived_name': col,
                               'source_columns': json.dumps(sources), 'role': role(parent, setting),
                               'output_location': 'outcomes_separate_file' if outcome_column(col) else 'cleaned_cohort',
                               'storage': 'numeric' if meta.get('numeric_companion') in {'numeric','bound'} else 'string_or_original_derived_value',
                               'value_labels': 'unresolved_without_machine_readable_codebook' if not meta else 'not_applicable_numeric_companion',
                               'numeric_companion': meta.get('numeric_companion', ''),
                               'companion_of': names.get(parent, '') if meta else '', 'unit': meta.get('unit', ''),
                               'screen_min': meta.get('screen_min', ''), 'screen_max': meta.get('screen_max', ''),
                               'note': ('Screening limits are review rules, not clinical reference intervals; raw values retained.' if meta
                                        else 'Source code/text retained; no code labels or binary semantics inferred.')})
        pd.DataFrame(dictionary).to_csv(OUTPUT / f'variable_dictionary_{setting}.csv', index=False)
        for subset_name, subset in [('all_first_visits', pd.Series(True, index=clean.index)), ('existing_primary_subset', primary)]:
            report = availability(cleaned_analysis, cleaned_states, source_states, subset, subset_name, setting, names)
            report.fillna(0).to_csv(OUTPUT / f'missingness_{setting}_{subset_name}.csv', index=False)
        # Controlled categorical code inventory; excludes narrative/free text.
        for col in snapshot:
            values = snapshot[col].dropna().astype(str)
            observed = values[~values.isin(TOKENS)]
            unique = sorted(observed.unique())
            if not outcome_column(col) and 0 < len(unique) <= 20 and all(re.fullmatch(r'-?\d+', v) for v in unique):
                for value, n in observed.value_counts().items():
                    inventories.append({'setting': setting, 'variable': names[col], 'original_or_derived_name': col,
                                        'code': value, 'n': int(n), 'label_status': 'unresolved_numeric_code_or_measurement_no_inferred_mapping'})
        assert len(clean) == int(mask.sum())
        assert clean['Participant Id'].tolist() == raw.loc[mask, 'Participant Id'].tolist()
        assert clean['Participant Id'].is_unique and clean.columns.is_unique
        assert cleaned_analysis.columns.equals(cleaned_states.columns)
        assert not any(outcome_column(c) for c in cleaned_analysis)
        path = OUTPUT / f'cohort_{setting}_cleaned.csv.gz'
        reread = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
        expected = missing_to_empty(cleaned_analysis.rename(columns=names)).astype(str)
        assert reread.shape == expected.shape
        assert reread.reset_index(drop=True).equals(expected.reset_index(drop=True)), 'String CSV roundtrip differs.'
        summaries[setting] = {'n_first_visits': len(clean), 'n_snapshot_columns': len(snapshot.columns),
                             'n_analysis_columns': len(cleaned_analysis.columns), 'n_numeric_measurement_sources': len(numeric_specs(clean.columns)),
                             'n_existing_primary_subset': int(primary.sum()), 'n_existing_positive_events': int(events.sum()),
                             'n_existing_primary_nonevents': int((primary & ~events).sum()),
                             'outcome_status_counts': {str(k): int(v) for k, v in snapshot['outcome_status'].value_counts().items()},
                             'index_date_missing': int(date_series(snapshot['index_date_iso']).isna().sum()),
                             'n_numeric_quality_flags': len(numeric_flags),
                             'n_duplicate_review_records': int(snapshot['duplicate_review_required'].sum()),
                             'n_existing_primary_controls_without_definitive_diagnosis': int((primary & ~events & snapshot['outcome_no_definitive_diagnosis'].eq(1)).sum()),
                             'same_day_linkage': {c: {str(k): int(v) for k,v in snapshot[c].value_counts(dropna=False).items()}
                                                  for c in snapshot if c.endswith('_link_status')},
                             'validation': {'ids_unique_and_order_preserved': True, 'snapshot_row_count_preserved': True,
                                            'analysis_status_columns_match': True, 'outcomes_separate_from_analysis': True,
                                            'exact_saved_string_roundtrip': True}}
    pd.DataFrame(all_audits).to_csv(OUTPUT / 'timing_audit.csv.gz', index=False, compression={'method':'gzip', 'mtime':0})
    audits = pd.DataFrame(all_audits)
    audits.groupby(['setting','block','status'], dropna=False).agg(n_records=('Participant Id','size'),
        source_values_present=('n_source_values_present','sum')).reset_index().to_csv(OUTPUT / 'timing_summary.csv', index=False)
    pd.DataFrame(all_flags).to_csv(OUTPUT / 'quality_flags.csv', index=False)
    pd.DataFrame(inventories).to_csv(OUTPUT / 'unresolved_code_inventory.csv', index=False)
    assert master.sha256(SOURCE) == source_hash
    assert master.sha256(ROOT / 'analysis' / '01_prepare_gallstone_step1.py') == old_hash
    owned_outputs = {'cohort_flow.csv', 'first_episodes_unassigned_setting.csv.gz',
                     'timing_audit.csv.gz', 'timing_summary.csv', 'quality_flags.csv', 'unresolved_code_inventory.csv'}
    for setting in ['SEH', 'POLI']:
        owned_outputs.update({
            f'cohort_{setting}_before_timing_filter.csv.gz', f'cohort_{setting}_cleaned.csv.gz',
            f'cohort_{setting}_source_missing_status.csv.gz', f'cohort_{setting}_analysis_missing_status.csv.gz',
            f'outcomes_{setting}_existing_definition.csv', f'variable_dictionary_{setting}.csv',
            f'missingness_{setting}_all_first_visits.csv', f'missingness_{setting}_existing_primary_subset.csv'})
    summary = {'stage': 'complete_first_visit_cleaning_before_imputation', 'imputation_performed': False,
               'predictor_selection_performed': False, 'source_sha256_before_and_after': source_hash,
               'existing_definition_script_sha256_unchanged': old_hash, 'script_sha256': master.sha256(Path(__file__)),
               'software': {'python': sys.version.split()[0], 'pandas': pd.__version__}, 'cohort_flow': flow, 'cohorts': summaries,
               'rules': ['All defined baseline and first-visit section fields retained; all linked same-day diagnostic records retained.',
                         'Same-day linkage proves calendar-date agreement, not availability before a clinical decision.',
                         'Historical RFU and prior-episode blocks require valid dates strictly before index; other timings withheld with distinct reasons.',
                         'Thirty-day complication histories require surgery-to-index interval of at least 30 days.',
                         'Undated GIQLI, study context and regular follow-up checkboxes withheld from cleaned analysis; source snapshot retained.',
                         'Visit final-diagnosis, diagnostic contribution-to-final-diagnosis questions and CT IH follow-up/surgery truth withheld for reference-standard timing review; assessment/process fields retained and labelled.',
                         'Source codes and ambiguous NO/no/No text preserved. Exact unambiguous textual missing terms harmonized only in analysis copy with source states retained. No inferred categorical labels or eligibility selection.',
                         'Numeric source strings preserved; numeric companions screened using documented review ranges, bounds/qualifiers retained.',
                         'Missing sentinel/text provenance and timing withholding preserved in separate aligned state matrices; empty diagnostic cells receive linkage context without erasing observed conflicting findings.',
                         'Existing gallstone outcome rule retained without reinterpretation; outcome validity/code mappings still require scientific review.',
                         'No imputation performed. Primary subset missingness uses existing outcome eligibility only, not a final analysis lock.',
                         'Read CSVs with keep_default_na=False to preserve legitimate literal NA text.'],
               'exact_text_missing_mapping': TEXT_MISSING_STATES,
               'ambiguous_tokens_preserved': ['NA', 'NO', 'no', 'No'],
               'duplicate_review_source': str(duplicate_path.relative_to(ROOT)) if duplicate_path.exists() else None,
               'outputs': {p.name: {'bytes': p.stat().st_size, 'sha256': master.sha256(p)} for p in sorted(OUTPUT.iterdir()) if p.is_file() and p.name in owned_outputs}}
    (OUTPUT / 'validation_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (OUTPUT / 'first_visit_metadata.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'output': str(OUTPUT), 'cohorts': {k: {key:value for key,value in v.items() if key not in ['same_day_linkage']} for k,v in summaries.items()}, 'flow': flow}, indent=2))

if __name__ == '__main__':
    main()
