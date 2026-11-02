#!/usr/bin/env python3
"""Audit provisional gallstone outcomes without changing any outcome or source value.

Run with the bundled Python runtime (pandas and numpy are required). The existing
step-1 definition is imported but its main() is never called. Numeric CRF category
labels remain unverified until the official value dictionary is supplied.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'OPERATE data' / 'Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
OUTPUT = ROOT / 'results' / 'first_visit_2026-09-30'
STEP1 = ROOT / 'analysis' / '01_prepare_gallstone_step1.py'
CODE_CAVEAT = (
    'The CRF prints option labels and their order, but not the stored numeric codes. '
    'FU_einde_hoe codes 2/3 and FU_einde_persisterende_klachten code 7 are therefore '
    'review triggers with unverified label mappings; they are not automatic exclusions.'
)
INTERPRETATION = (
    'FU_einde_diagnose records whether a definitive diagnosis was established. '
    'An answered No does not itself verify absence of gallstone disease. Existing '
    'controls remain provisional; no outcome is recoded or imputed by this audit.'
)


def strict_dates(values: pd.Series) -> pd.Series:
    values = values.astype('string').str.strip()
    parsed = pd.Series(pd.NaT, index=values.index, dtype='datetime64[ns]')
    for pattern, fmt in ((r'\d{4}-\d{2}-\d{2}', '%Y-%m-%d'), (r'\d{2}-\d{2}-\d{4}', '%d-%m-%Y')):
        mask = values.str.fullmatch(pattern, na=False)
        parsed.loc[mask] = pd.to_datetime(values.loc[mask], format=fmt, errors='coerce')
    return parsed


def support_dates(raw: pd.DataFrame, derived: pd.DataFrame, setting: str, step1) -> pd.DataFrame:
    """Replicate the existing qualifying evidence conditions, retaining dates.

    This does not introduce a new evidence criterion or follow-up cutoff.
    """
    index = strict_dates(raw[f'{setting}_datum'])
    evidence: dict[str, pd.Series] = {}
    for date_col, finding_col in step1.ECHO_CHOLECYSTITIS_EVIDENCE + step1.CT_CHOLECYSTITIS_EVIDENCE:
        dates = strict_dates(raw[date_col])
        condition = raw[finding_col].eq('1') & derived.outcome_candidate_cholecystitis.eq(1)
        evidence[finding_col] = dates.where(condition & dates.ge(index))
    for performed, date_col, type_col, other_col, diagnosis_col, diagnosis_text_col in step1.SURGERY_EVIDENCE:
        dates = strict_dates(raw[date_col])
        operation = raw[type_col].isin(['2', '3']) | raw[other_col].str.contains(step1.CHOLE_OPERATION_TEXT)
        gallstones = operation & derived.outcome_candidate_cholecystolithiasis.eq(1)
        inflammation = (
            operation | raw[diagnosis_col].eq('6')
            | raw[diagnosis_text_col].str.contains(step1.CHOLECYSTITIS_FINAL_TEXT)
        ) & derived.outcome_candidate_cholecystitis.eq(1)
        evidence[performed] = dates.where(raw[performed].eq('1') & (gallstones | inflammation) & dates.ge(index))
    for date_col, abnormal_col, report_col, diagnosis_col in step1.ERCP_EVIDENCE:
        dates = strict_dates(raw[date_col])
        condition = (
            raw[abnormal_col].eq('1') & raw[diagnosis_col].eq('1')
            & raw[report_col].str.contains(step1.ERCP_STONE_TEXT)
            & derived.outcome_candidate_cholecystolithiasis.eq(1)
        )
        evidence[date_col] = dates.where(condition & dates.ge(index))
    return pd.DataFrame(evidence, index=raw.index)


def main() -> None:
    spec = importlib.util.spec_from_file_location('operate_step1_outcomes', STEP1)
    step1 = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(step1)
    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    raw = pd.read_csv(SOURCE, sep=';', dtype=str, encoding='utf-8-sig', keep_default_na=False)
    assert raw['Participant Id'].is_unique
    output_rows: list[dict] = []
    summaries = {}
    validation = {}
    diagnosis_dates = ['definitieve_diagnose_datum', 'definitieve_diagnose_datum_2', 'Definiteve_diagnose_datum_3']
    for setting in ('SEH', 'POLI'):
        cohort = raw.loc[raw.episodes.eq('1') & raw[f'{setting}_1stecontact'].eq('1')].copy()
        outcome = step1.derive_outcome(cohort, setting)
        reference_path = ROOT / 'results' / 'step1_gallstone' / f'cohort_{setting}_first_presentations.csv'
        reference = pd.read_csv(reference_path, dtype=str, encoding='utf-8-sig', keep_default_na=False).set_index('Participant Id')
        current = pd.DataFrame({
            'outcome_status': outcome.outcome_status.to_numpy(),
            'gallstone_related_disease': outcome.gallstone_related_disease.astype('string').fillna('').to_numpy(),
        }, index=cohort['Participant Id'])
        assert set(current.index) == set(reference.index), 'First-presentation ID sets changed versus existing step1 outputs'
        for name in current.columns:
            assert current[name].equals(reference.loc[current.index, name]), f'Outcome changed versus step1: {setting} {name}'
        validation[setting] = {'reference_cohort': str(reference_path.relative_to(ROOT)), 'id_sets_match': True, 'all_existing_outcomes_and_statuses_match': True}
        control = outcome.gallstone_related_disease.eq(0).fillna(False)
        positive = outcome.gallstone_related_disease.eq(1).fillna(False)
        primary = control | positive
        index = strict_dates(cohort[f'{setting}_datum'])
        end = strict_dates(cohort.FU_einde_datum)
        supports = support_dates(cohort, outcome, setting, step1)
        first_support = supports.min(axis=1)
        assert first_support.loc[positive].notna().all(), 'Positive evidence dates could not be reproduced'
        new_codes = cohort[step1.DX_NEW_COLUMNS].astype(str).apply(lambda values: values.str.strip())
        new_first = new_codes.iloc[:, 0]
        old_first = cohort['FU_einde_diagnose_definitief_ordinaal'].str.strip()
        old_gallstone = old_first.isin(['1', '3'])
        new_first_special = new_first.isin(['-95', '-96', '-97', '-98', '-99'])
        legacy_conflict = old_gallstone & new_first.ne('') & ~new_first.isin(['5', '6']) & ~new_first_special
        new_later_gallstone = new_codes.iloc[:, 1:].isin(['5', '6']).any(axis=1)
        flags = {
            'legacy_gallstone_vs_new_first_other_code': legacy_conflict & ~new_later_gallstone,
            'legacy_gallstone_with_new_first_special_code_fallback_blocked': old_gallstone & new_first_special,
            'control_no_definitive_diagnosis_flag_0': control & cohort.FU_einde_diagnose.eq('0'),
            'control_followup_route_code_2_or_3_mapping_unverified': control & cohort.FU_einde_hoe.isin(['2', '3']),
            'control_persistent_symptoms_code_7_mapping_unverified': control & cohort.FU_einde_persisterende_klachten.eq('7'),
            'control_followup_route_missing_or_special': control & cohort.FU_einde_hoe.isin(['', '-95', '-96', '-97', '-98', '-99']),
            'index_date_missing_or_invalid': primary & index.isna(),
            'followup_end_date_missing_or_invalid': primary & end.isna(),
            'followup_end_before_index': primary & end.notna() & index.notna() & end.lt(index),
            'positive_first_support_after_followup_end': positive & end.notna() & first_support.gt(end),
            'positive_first_diagnosis_flag_no': positive & cohort.FU_einde_diagnose.eq('0'),
        }
        for number, name in enumerate(diagnosis_dates, start=1):
            dates = strict_dates(cohort[name])
            flags[f'diagnosis_{number}_date_before_index'] = primary & dates.notna() & index.notna() & dates.lt(index)
            flags[f'diagnosis_{number}_date_after_followup_end'] = primary & dates.notna() & end.notna() & dates.gt(end)
        counts = Counter()
        for row_id in cohort.index:
            reasons = [name for name, mask in flags.items() if bool(mask.loc[row_id])]
            if not reasons:
                continue
            counts.update(reasons)
            source_row = cohort.loc[row_id]
            record = {
                'Participant Id': source_row['Participant Id'],
                'setting': setting,
                'gallstone_related_disease_unchanged': int(outcome.loc[row_id, 'gallstone_related_disease']) if pd.notna(outcome.loc[row_id, 'gallstone_related_disease']) else '',
                'outcome_status_unchanged': outcome.loc[row_id, 'outcome_status'],
                'index_date_raw': source_row[f'{setting}_datum'],
                'followup_end_date_raw': source_row['FU_einde_datum'],
                'followup_route_raw_code': source_row['FU_einde_hoe'],
                'first_definitive_diagnosis_present_raw_code': source_row['FU_einde_diagnose'],
                'persistent_symptoms_raw_code': source_row['FU_einde_persisterende_klachten'],
                'first_qualifying_support_date': first_support.loc[row_id].strftime('%Y-%m-%d') if pd.notna(first_support.loc[row_id]) else '',
                'qualifying_evidence_dates': '; '.join(f'{name}={value:%Y-%m-%d}' for name, value in supports.loc[row_id].items() if pd.notna(value)),
                'review_reasons': '; '.join(reasons),
                'review_reason_count': len(reasons),
                'action': 'clinical_or_source_review; no automatic outcome change',
            }
            record['legacy_first_diagnosis_raw_code'] = old_first.loc[row_id]
            for number, name in enumerate(step1.DX_NEW_COLUMNS, start=1):
                record[f'new_diagnosis_{number}_raw_code'] = source_row[name]
            for number, name in enumerate(diagnosis_dates, start=1):
                record[f'diagnosis_{number}_date_raw'] = source_row[name]
            output_rows.append(record)
        summaries[setting] = {
            'first_presentations_n': len(cohort),
            'fallback_audit': {
                'new_first_special_code_n': int(new_first_special.sum()),
                'legacy_gallstone_new_first_special_n': int((old_gallstone & new_first_special).sum()),
                'legacy_gallstone_new_first_other_n': int(legacy_conflict.sum()),
                'of_conflicts_new_second_or_third_gallstone_n': int((legacy_conflict & new_later_gallstone).sum()),
                'new_first_blank_legacy_gallstone_new_later_present_n': int((new_first.eq('') & old_gallstone & new_codes.iloc[:, 1:].ne('').any(axis=1)).sum()),
            },
            'provisional_positive_n': int(positive.sum()),
            'provisional_control_n': int(control.sum()),
            'records_flagged_for_review_n': sum(any(bool(mask.loc[i]) for mask in flags.values()) for i in cohort.index),
            'review_reason_counts': dict(counts),
            'review_reason_counts_controls': {k: int((v & control).sum()) for k, v in flags.items()},
            'review_reason_counts_positives': {k: int((v & positive).sum()) for k, v in flags.items()},
            'controls_route_2_3_or_symptoms_7_n': int((control & (cohort.FU_einde_hoe.isin(['2', '3']) | cohort.FU_einde_persisterende_klachten.eq('7'))).sum()),
        }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    review = pd.DataFrame(output_rows)
    assert not review.duplicated(['Participant Id', 'setting']).any()
    review.to_csv(OUTPUT / 'outcome_ascertainment_review.csv', index=False, encoding='utf-8-sig')
    saved = pd.read_csv(OUTPUT / 'outcome_ascertainment_review.csv', dtype=str, encoding='utf-8-sig', keep_default_na=False)
    for setting, details in summaries.items():
        rows = saved.loc[saved.setting.eq(setting)]
        assert len(rows) == details['records_flagged_for_review_n']
        reasons = Counter(reason for text in rows.review_reasons for reason in text.split('; '))
        assert dict(reasons) == details['review_reason_counts']
    validation['output_roundtrip_counts_match'] = True
    validation['duplicate_setting_participant_rows_n'] = int(saved.duplicated(['Participant Id', 'setting']).sum())
    summary = {
        'source': str(SOURCE.relative_to(ROOT)), 'source_sha256': source_hash,
        'outcome_rule_module': str(STEP1.relative_to(ROOT)),
        'outcome_rule_module_sha256': hashlib.sha256(STEP1.read_bytes()).hexdigest(),
        'outcome_rule_changed': False, 'outcome_imputed': False,
        'review_scope': 'Existing primary positive/control records for ascertainment and dates; diagnosis-code conflicts across all first presentations. Review flags are not exclusions.',
        'interpretation': INTERPRETATION, 'code_mapping_caveat': CODE_CAVEAT,
        'crf_reference': 'OPERATE data/empty_crf.pdf, PDF pages 730-732 (printed pages 729-731).',
        'date_logic': 'Explicit dd-mm-yyyy or yyyy-mm-dd parsing; special/sentinel/unparseable dates unavailable. Supporting evidence duplicates the existing rule and is not truncated at FUend; contradictions are flagged.',
        'settings': summaries, 'validation': validation,
    }
    (OUTPUT / 'outcome_ascertainment_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == source_hash
    print(json.dumps({'review_rows': len(review), 'rows_by_setting': review.groupby('setting').size().to_dict(), 'outcomes_match_existing_step1': True, 'output_counts_verified': True, 'source_unchanged': True}, ensure_ascii=False))


if __name__ == '__main__':
    main()
