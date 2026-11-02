#!/usr/bin/env python3
"""Build protected clinical review queues from corrected v3 saved outputs.

No raw data edits, imputation, model fitting, correction or new exclusions.
ASA, urine SG and unresolved-unit eGFR are held field-wide: valid-looking ASA
observations are held pending coding review, not declared erroneous. Invalid
urine pH cells are held individually. Diastolic BP198 is review-only.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/first_visit_2026-10-27_v3'
QUARANTINE = {'classification_coding_review_required', 'measurement_scale_or_unit_review_required'}


def read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def numeric(value):
    text = str(value).strip().replace(',', '.')
    if not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)', text):
        return None
    return float(text)


def main():
    queue, inputs, settings = [], [], {}
    for setting in ('SEH', 'POLI'):
        clean_path = OUT / f'cohort_{setting}_cleaned.csv.gz'
        state_path = OUT / f'cohort_{setting}_analysis_missing_status.csv.gz'
        dict_path = OUT / f'variable_dictionary_analysis_{setting}.csv'
        clean, states, dictionary = map(read, (clean_path, state_path, dict_path))
        inputs.extend((clean_path, state_path, dict_path))
        if not clean.participant_id.is_unique or not states.participant_id.is_unique:
            raise ValueError('Review queue requires unique first-episode identifiers')
        states = states.set_index('participant_id').loc[clean.participant_id].reset_index()
        field_counts = Counter()
        source_counts = Counter()
        all_scopes = Counter()
        for _, d in dictionary.iterrows():
            target, source = d.variable, d.source_variable
            if target not in clean or target not in states or source not in clean:
                raise ValueError(f'Missing dictionary-linked column: {setting}/{target}')
            for i in clean.index:
                status, value = states.at[i, target], clean.at[i, source]
                number = numeric(value)
                reason = action = scope = ''
                if status in QUARANTINE:
                    scope = 'field_wide_hold'
                    action = 'source_owner_review_before_any_target_or_auxiliary_use'
                    if source == 'b_asa':
                        reason = 'ASA coding unresolved; all observations held, including valid-looking values; not all values declared wrong'
                    elif re.fullmatch(r'index_urine(?:_repeat\d+)?_sg', source):
                        reason = 'Urine specific-gravity scales mixed; no automatic scale conversion'
                    elif 'egfr' in source:
                        reason = 'eGFR unit kl/1.73m2 unresolved; no assumed unit correction'
                    else:
                        reason = 'Saved analytical cell explicitly quarantined pending clinical review'
                    if clean.at[i, target] != '':
                        raise ValueError(f'Quarantined analytical field contains a value: {target}')
                elif re.fullmatch(r'index_urine(?:_repeat\d+)?_ph', source) and number is not None and (number < 0 or number > 14):
                    reason = 'Urine pH outside physical range0–14; source retained without guessed decimal correction'
                    scope = 'individual_cell_hold'
                    action = 'source_owner_review_before_analytical_use'
                    if clean.at[i, target] != '' or status != 'numeric_plausibility_review':
                        raise ValueError('Invalid urine pH not held in corrected analytical representation')
                elif source in ('seh_lo_drr', 'poli_lo_drr') and number == 198:
                    reason = 'Diastolic blood pressure198 mmHg requires clinical/source review; no automatic exclusion or correction'
                    scope = 'review_only_no_exclusion'
                    action = 'review_original_record_value_retained_under_current_rules'
                if not reason:
                    continue
                field_counts[source] += 1
                source_counts[source] += int(bool(value.strip()))
                all_scopes[scope] += 1
                if not value.strip():
                    continue
                queue.append(dict(participant_id=clean.at[i,'participant_id'], setting=setting,
                                  analytical_variable=target, source_variable=source,
                                  original_source_names=d.source_names,
                                  source_value=value, source_present=bool(value.strip()),
                                  source_status=states.at[i,source] if source in states else 'not_available',
                                  analytical_value=clean.at[i,target], analytical_status=status,
                                  unit=d.unit, review_scope=scope, reason=reason, required_action=action,
                                  correction_applied=False,
                                  eligibility_changed_by_queue=False))
        settings[setting] = dict(first_visit_records=len(clean),
                                 saved_hold_or_review_cells_by_source=dict(field_counts),
                                 actionable_nonempty_source_cells_by_source=dict(source_counts),
                                 blank_source_held_cells_by_source={k:field_counts[k]-source_counts[k] for k in field_counts},
                                 saved_hold_or_review_scope_counts=dict(all_scopes))
    columns = ['participant_id','setting','analytical_variable','source_variable','original_source_names',
               'source_value','source_present','source_status','analytical_value','analytical_status','unit',
               'review_scope','reason','required_action','correction_applied','eligibility_changed_by_queue']
    frame = pd.DataFrame(queue, columns=columns)
    frame.to_csv(OUT / 'clinical_review_queue.csv', index=False)
    asa = frame.loc[frame.source_variable.eq('b_asa'), 'source_value'].map(numeric)
    asa_valid_looking = asa.isin([1,2,3,4,5,6])
    summary = dict(scope='Protected clinical review queue from corrected v3 saved outputs; no changes to analyses/data',
                   settings=settings, queued_cells=len(frame),
                   queued_nonempty_source_cells=int(frame.source_present.sum()) if len(frame) else 0,
                   review_scope_counts=frame.review_scope.value_counts().to_dict(),
                   ASA_nonempty_held_observations=len(asa),
                   ASA_valid_looking_class_values_held_not_declared_wrong=int(asa_valid_looking.sum()),
                   ASA_other_observed_values_requiring_coding_review=int((~asa_valid_looking).sum()),
                   interpretation=['Queue contains actionable nonempty-source cells only; field-wide held/blank totals are separate in this summary.',
                                   'Blank-source quarantined cells remain fully traceable in saved status matrices.',
                                   'Valid-looking ASA observations are held because coding is unresolved, not declared incorrect.',
                                   'Source values and original source/missingness states remain traceable.',
                                   'DiastolicBP198 is review-only and causes no new exclusion.',
                                   'No guessed decimal, scale or unit conversions are applied.'],
                   input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT / 'clinical_review_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('input_sha256','script_sha256','settings','interpretation')}))


if __name__ == '__main__':
    main()
