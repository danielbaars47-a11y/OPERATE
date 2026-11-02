#!/usr/bin/env python3
"""Versioned full first-visit cleaning and cell-level MI specification.

Preserves the September source snapshots and all initial fields. No unverified
category meanings or unperformed laboratory measurements are manufactured.
Outcome candidates are reproducible but explicitly NOT a frozen endpoint until
the current diagnostic/procedure mappings have been checked.
"""
from pathlib import Path
import importlib.util
import hashlib
import json
import re
import shutil
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PREV = ROOT / 'results/first_visit_2026-09-30'
OUT = ROOT / 'results/first_visit_2026-10-27'
SOURCE = ROOT / 'OPERATE data/Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
SOURCE_HASH = '833b9842d96274ac56fb73998189a07d4e2ed03375de1b9a5ff4082287585911'
CHECKED = 'historical_mapping_current_question_and_options_concordant'
GENUINE = {'asked_unknown', 'not_described', 'not_recorded', 'measurement_failed'}
OBSERVED = {'present_code_or_text', 'numeric_observed', 'derived_index_date_measurement'}

def read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)

def write(frame, name):
    p = OUT / name
    frame.to_csv(p, index=False, na_rep='', compression={'method':'gzip','mtime':0} if str(p).endswith('.gz') else None)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def module(filename, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'analysis' / filename)
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj

base = module('03_prepare_operate_first_visit.py', 'operate_clean_v1')
old = base.old

def endpoint(raw, setting, clean):
    """Revised clinical rule, retaining explicit unresolved mapping status."""
    raw = raw.set_index('Participant Id', drop=False).loc[clean.participant_id].reset_index(drop=True)
    raw = raw.apply(lambda s:s.str.strip())
    idx = base.date_series(raw[f'{setting}_datum'])
    end = base.date_series(raw.FU_einde_datum)
    def dated(frame, dates, column, condition):
        evidence_date = base.date_series(frame[column])
        return condition.fillna(False) & dates.notna() & end.notna() & evidence_date.ge(dates) & evidence_date.le(end)
    old.evidence_on_or_after_index = dated
    # Uncertain/negated narrative is not a definitive diagnosis label. Keep it
    # in the source and audit, but do not let a substring create an event.
    narrative_uncertain = pd.Series(False,index=raw.index)
    rule_raw = raw.copy()
    for col in old.DX_TEXT_COLUMNS:
        uncertain = raw[col].str.contains(r'(?:geen|niet|uitgesloten|verdenking|mogelijk|\?).{0,35}(?:galste|cholecyst|cholelith)',case=False,regex=True)
        narrative_uncertain |= uncertain
        rule_raw.loc[uncertain,col] = ''
    out = old.derive_outcome(rule_raw, setting)
    positive = out.outcome_definition_met_cholecystolithiasis.eq(1) | out.outcome_definition_met_cholecystitis.eq(1)
    labelled = out.outcome_candidate_cholecystolithiasis.eq(1) | out.outcome_candidate_cholecystitis.eq(1)
    other = out.outcome_candidate_other_gallstone_related.eq(1)
    assessed = raw.FU_einde_diagnose.isin(['0','1'])
    conflict = idx.isna() | (end.notna() & end.lt(idx))
    dateflags = []
    for col in ['definitieve_diagnose_datum','definitieve_diagnose_datum_2','Definiteve_diagnose_datum_3']:
        dt = base.date_series(raw[col])
        flag = dt.notna() & ((idx.notna() & dt.lt(idx)) | (end.notna() & dt.gt(end)))
        conflict |= flag
        dateflags.append(flag)
    conflict |= raw.FU_einde_diagnose.eq('0') & labelled
    new = raw[old.DX_NEW_COLUMNS]
    legacy_conflict = raw.FU_einde_diagnose_definitief_ordinaal.isin(['1','3']) & new.iloc[:,0].ne('') & ~new.isin(['5','6']).any(axis=1)
    conflict |= legacy_conflict
    conflict |= narrative_uncertain & labelled
    candidate = pd.Series('',index=raw.index,dtype='string')
    candidate.loc[assessed & ~labelled & ~other & ~conflict] = '0'
    candidate.loc[positive & ~conflict] = '1'
    prior = raw.VG_cholecystectomie.isin(['2','3','4'])
    prior_unknown = ~raw.VG_cholecystectomie.isin(['1','2','3','4'])
    duplicate = clean.duplicate_review_required.eq('1')
    reason = pd.Series('',index=raw.index,dtype='string')
    reason.loc[~assessed & ~positive] = 'outcome_assessment_missing'
    reason.loc[labelled & ~positive] = 'definitive_label_without_verifiable_in_episode_evidence'
    reason.loc[other & ~positive & ~labelled] = 'other_possible_gallstone_diagnosis_review'
    reason.loc[conflict] = 'outcome_chronology_or_source_conflict'
    reason.loc[prior_unknown] = 'prior_cholecystectomy_history_unresolved'
    reason.loc[duplicate] = 'possible_duplicate_person_identity_unresolved'
    reason.loc[prior] = 'prior_cholecystectomy_excluded_both_case_and_control'
    clinically_eligible = candidate.ne('') & reason.eq('')
    result = pd.DataFrame({'participant_id':clean.participant_id,'analysis_setting':setting,
        'index_date':idx.dt.strftime('%Y-%m-%d').fillna(''),
        'followup_end_date':end.dt.strftime('%Y-%m-%d').fillna(''),
        'followup_rule':'recorded_CRF134_end_after_3_months_without_contact_not_index_plus_90_days',
        'candidate_endpoint_before_exclusions':candidate,
        'candidate_endpoint_after_exclusions':candidate.where(clinically_eligible,''),
        'gallstone_related_disease':'',
        'outcome_mapping_status':'pending_current_diagnostic_and_procedure_code_verification',
        'outcome_assessment_recorded':assessed.astype(int),
        'uncertain_or_negated_diagnosis_narrative_flag':narrative_uncertain.astype(int),
        'unclear_or_no_definitive_diagnosis_included_as_candidate_zero':(candidate.eq('0') & raw.FU_einde_diagnose.eq('0')).astype(int),
        'prior_cholecystectomy_excluded':prior.astype(int),
        'prior_cholecystectomy_history_unknown':prior_unknown.astype(int),
        'identity_review_required':duplicate.astype(int),
        'clinical_eligibility_pending_code_verification':clinically_eligible.astype(int),
        'analysis_eligible':0,'review_reason':reason})
    for c in out:
        if c.startswith('outcome_'): result['audit_'+c] = out[c]
    return result

def block_for(name):
    if name.startswith('index_lab'): return 'laboratory'
    if '_lo_' in name: return 'physical_examination'
    if '_a_' in name: return 'symptoms'
    return 'baseline'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert digest(SOURCE)==SOURCE_HASH, 'Source changed; inspect before reusing prior cleaning.'
    receipt = json.loads((PREV/'first_visit_metadata.json').read_text())
    for name, info in receipt['outputs'].items():
        assert digest(PREV/name)==info['sha256'], f'Previous output changed: {name}'
    raw = pd.read_csv(SOURCE,sep=';',dtype=str,keep_default_na=False,na_filter=False)
    checked = {x['source_name']:x for x in json.loads((OUT/'codebook/reviewed_mappings.json').read_text())}
    crf_pages = json.loads((OUT/'codebook/crf_pages.json').read_text())
    # This export contains no SEH discharge timestamp or visit identifier for tests.
    # A next-day result cannot automatically be attributed to the index visit.
    overnight=[]; summaries={}; availability=[]; allspec=[]
    source_hashes={str(SOURCE.relative_to(ROOT)):digest(SOURCE)}
    for setting in ['SEH','POLI']:
        clean=read(PREV/f'cohort_{setting}_cleaned.csv.gz')
        states=read(PREV/f'cohort_{setting}_analysis_missing_status.csv.gz')
        dictionary=read(PREV/f'variable_dictionary_{setting}.csv')
        dictionary=dictionary[dictionary.output_location.eq('cleaned_cohort')].copy()
        assert clean.participant_id.equals(states.participant_id) and clean.columns.equals(states.columns)
        original_columns=list(clean)
        outcomes=endpoint(raw,setting,clean)
        spec=[]; inputs={'participant_id':clean.participant_id,'gallstone_related_disease':outcomes.gallstone_related_disease}
        wheres={'participant_id':clean.participant_id}; eligible={'participant_id':clean.participant_id}
        eligible_rows=~outcomes.prior_cholecystectomy_excluded.eq(1) & ~outcomes.prior_cholecystectomy_history_unknown.eq(1) & ~outcomes.identity_review_required.eq(1)
        metadata=dictionary.set_index('variable')
        original_map=dict(zip(dictionary.original_or_derived_name,dictionary.variable))
        # Append labelled and analytical companions; source codes remain intact.
        added_metadata=[]
        for original, info in checked.items():
            if original not in original_map or original in {'episodes','SEH_1stecontact','FU_einde_diagnose','VG_cholecystectomie'}: continue
            src=original_map[original]; codes=json.loads(info['codes_json'])
            substantive={k:v for k,v in codes.items() if not k.startswith('-')}
            unknown={k for k,v in substantive.items() if v.strip().lower() in {'onbekend','onduidelijk','niet omschreven'}}
            name=src+'_analysis'; label=src+'_label'
            clean[label]=clean[src].map(substantive).fillna('')
            states[label]=states[src]
            clean[name]=clean[src].where(clean[src].isin(set(substantive)-unknown),'')
            states[name]=states[src]
            states.loc[clean[src].isin(unknown),name]='asked_unknown'
            valid_codes=set(substantive)-unknown
            applicable=eligible_rows & states[name].isin(GENUINE|{'blank_or_dependency_not_triggered','present_code_or_text'})
            method='logreg' if len(valid_codes)==2 else 'polyreg'
            for n,kind in [(name,'categorical_analysis_companion'),(label,'value_label_companion')]:
                item=metadata.loc[src].to_dict(); item.update(variable=n,original_or_derived_name=original,
                    source_columns=json.dumps([original]),storage=kind,companion_of=src,value_labels=json.dumps(substantive),
                    note='Field-specific question and options checked against current CRF; numeric mapping provenance 2024 SPSS.')
                added_metadata.append(item)
            inputs[name]=clean[name]
            wheres[name]=applicable & clean[name].eq('')
            eligible[name]=applicable
            spec.append(dict(variable=name,method=method,block=block_for(src),levels_json=json.dumps(sorted(valid_codes)),
                reason='Verified unconditional question; explicit unknown treated as missing, source preserved',verification_state=CHECKED))
            # Sex is a complete auxiliary if available; keep a canonical alias.
            if original=='Geslacht': inputs['geslacht']=clean[name]
        dictionary=pd.concat([dictionary,pd.DataFrame(added_metadata)],ignore_index=True)
        numeric=metadata[metadata.numeric_companion.eq('numeric')]
        for name, info in numeric.iterrows():
            source=info.companion_of; unit=info.unit
            reason=''; method='pmm'; applicable=pd.Series(False,index=clean.index)
            if re.search(r'_repeat\d+_',name): reason='Repeated test retained; earliest availability/order not established'
            elif any(t in unit for t in ['not_specified','unresolved']): reason='Measurement unit unverified; no imputation or conversion'
            elif name.endswith(('_bmi_numeric','_twl_numeric')) or name=='start_bmi_numeric': reason='Derived quantity; do not impute independently of components'
            elif name.startswith('index_lab_'):
                # Analyte-level explicit failed/unknown result establishes an attempted
                # measurement; a blank panel field does NOT establish test ordering.
                link=clean.get('index_lab_link_status',pd.Series('',index=clean.index))
                linked=link.isin(['matched','multiple_same_date_records'])
                applicable=eligible_rows & linked & (states[name].eq('numeric_observed') | states[name].isin(GENUINE))
                reason='Only analyte-specific failed/unknown entries within concordant index panel; blank assay ordering unresolved'
            elif name in {'startgewicht_numeric','lengte_numeric'}:
                applicable=eligible_rows & states[name].isin(GENUINE|{'numeric_observed','blank_or_dependency_not_triggered'})
                reason='Unconditional baseline measurement; preserve not-done/not-applicable'
            elif '_lo_' in name:
                # Explicit unknown/failed measurement only; blank vital fields may
                # represent examinations not undertaken (especially telephone POLI).
                applicable=eligible_rows & (states[name].eq('numeric_observed') | states[name].isin(GENUINE))
                reason='Explicit missing measurement only; blank physical-exam performance unresolved'
            else: reason='Applicability not established'
            if not applicable.any(): method=''
            inputs[name]=clean[name]
            wheres[name]=applicable & clean[name].eq('')
            eligible[name]=applicable
            spec.append(dict(variable=name,method=method,block=block_for(name),levels_json='',reason=reason,
                verification_state='measured_numeric_unit_checked' if method else 'not_eligible'))
        for name in ['age_years_at_index','years_since_bariatric_surgery']:
            inputs[name]=clean[name]
        # Safe deterministic companions; never overwrite reported BMI.
        height=pd.to_numeric(clean.lengte_numeric,errors='coerce')/100
        weight=pd.to_numeric(clean.startgewicht_numeric,errors='coerce')
        clean['baseline_bmi_calculated']=weight/height.pow(2)
        states['baseline_bmi_calculated']=np.where(clean.baseline_bmi_calculated.notna(),'derived_from_observed_components','components_missing_not_independently_imputed')
        dictionary=pd.concat([dictionary,pd.DataFrame([dict(variable='baseline_bmi_calculated',setting=setting,
            original_or_derived_name='baseline_bmi_calculated',source_columns='["Startgewicht", "Lengte"]',
            role='derived_baseline_measurement',output_location='cleaned_cohort',storage='numeric',unit='kg/m2',
            note='Startgewicht/(Lengte/100)^2; reported BMI preserved')])],ignore_index=True).fillna('')
        for name in ['prior_cholecystectomy_excluded','identity_review_required','analysis_eligible']:
            clean[name]=outcomes[name].astype(str); states[name]='administrative_or_selection'
            dictionary=pd.concat([dictionary,pd.DataFrame([dict(variable=name,setting=setting,role='administrative_or_selection',
                output_location='cleaned_cohort',note='Not a predictor; never imputed')])],ignore_index=True).fillna('')
        rawsub=raw.set_index('Participant Id').loc[clean.participant_id].reset_index(drop=True)
        index=base.date_series(rawsub[f'{setting}_datum'])
        groups=base.build_groups(list(raw.columns))
        for modality, instances in groups.items():
            for group in instances:
                test_date=base.date_series(rawsub[group['date']])
                mask=test_date.eq(index+pd.Timedelta(days=1)) & index.notna()
                for i in rawsub.index[mask]: overnight.append(dict(participant_id=clean.participant_id[i],setting=setting,
                    modality=modality,source_instance=group['selector'],index_date=index[i].strftime('%Y-%m-%d'),
                    test_date=test_date[i].strftime('%Y-%m-%d'),status='next_day_candidate_no_verified_visit_end_or_link',
                    action='Retained in raw/master; not automatically included or imputed'))
        assessed_names={s['variable'] for s in spec}
        for _,item in dictionary.iterrows():
            n=item.variable
            if n not in assessed_names:
                spec.append(dict(variable=n,method='',block=block_for(n),levels_json='',
                    reason=('Administrative/derived/text/source representation or unverified coding/applicability; retained, not stochastically imputed'),
                    verification_state='not_eligible'))
        specframe=pd.DataFrame(spec)
        for col in clean:
            if col in {'participant_id','analysis_setting'}: continue
            vals=clean[col].fillna('').astype(str); st=states[col]
            counts=st[vals.eq('')].value_counts().to_dict()
            availability.append(dict(setting=setting,variable=col,n_patients=len(clean),observed_n=int(vals.ne('').sum()),
                observed_pct=round(float(vals.ne('').mean()*100),3),missing_reason_counts_json=json.dumps(counts),
                eligible_missing_n=int(pd.Series(wheres.get(col,False),index=clean.index).sum())))
        write(clean,f'cohort_{setting}_cleaned.csv.gz');write(clean,f'cohort_{setting}_cleaned.csv')
        write(states,f'cohort_{setting}_analysis_missing_status.csv.gz')
        write(dictionary,f'variable_dictionary_{setting}.csv');write(outcomes,f'outcomes_{setting}_review.csv')
        write(specframe,f'imputation_spec_{setting}.csv');write(pd.DataFrame(inputs),f'imputation_input_{setting}.csv.gz')
        write(pd.DataFrame(wheres),f'where_{setting}.csv.gz');write(pd.DataFrame(eligible),f'eligibility_{setting}.csv.gz')
        for suffix in ['before_timing_filter','source_missing_status']:
            shutil.copy2(PREV/f'cohort_{setting}_{suffix}.csv.gz',OUT/f'cohort_{setting}_{suffix}.csv.gz')
        assert list(clean.columns[:len(original_columns)])==original_columns
        assert len(clean)==len(rawsub) and clean.participant_id.is_unique
        orig=read(PREV/f'cohort_{setting}_cleaned.csv.gz')
        assert clean[original_columns].equals(orig), 'Existing source/analytical cells changed unexpectedly'
        assert states.columns.equals(clean.columns)
        assert outcomes.loc[outcomes.prior_cholecystectomy_excluded.eq(1),'candidate_endpoint_after_exclusions'].eq('').all()
        assert outcomes.loc[outcomes.prior_cholecystectomy_excluded.eq(1),'analysis_eligible'].eq(0).all()
        summaries[setting]=dict(first_presentations=len(clean),columns=len(clean.columns),
            prior_cholecystectomy_excluded_both_groups=int(outcomes.prior_cholecystectomy_excluded.sum()),
            identity_review_required=int(outcomes.identity_review_required.sum()),
            candidate_events_not_frozen=int(outcomes.candidate_endpoint_after_exclusions.eq('1').sum()),
            candidate_nonevents_not_frozen=int(outcomes.candidate_endpoint_after_exclusions.eq('0').sum()),
            no_definitive_diagnosis_retained_as_candidate_zero=int(outcomes.unclear_or_no_definitive_diagnosis_included_as_candidate_zero.sum()),
            target_fields_with_eligible_missing=sum(bool(pd.Series(v).any()) for k,v in wheres.items() if k!='participant_id'),
            eligible_missing_cells=sum(int(pd.Series(v).sum()) for k,v in wheres.items() if k!='participant_id'))
    for name in ['first_episodes_unassigned_setting.csv.gz','cohort_flow.csv','quality_flags.csv','timing_summary.csv']:
        shutil.copy2(PREV/name,OUT/name)
    write(pd.DataFrame(overnight),'overnight_visit_linkage_review.csv')
    write(pd.DataFrame(availability),'availability_all_variables.csv')
    configuration=dict(version='2026-10-27',outcome_mapping_status='pending_current_code_verification',
        clinical_imputation_run=False,m=50,maxit=30,maxit_if_unstable=100,
        primary_priority='clean_and_impute_both_settings_no_prediction_modelling',
        prior_cholecystectomy_policy='exclude_both_case_and_control_retain_clean_master',
        required_before_clinical_MI=['Current diagnosis category and procedure code mappings',
            'Resolve identity review for records to be included in primary analysis'],
        note='A populated candidate endpoint is not a validated outcome. No clinical MICE runs until the outcome mapping is validated.')
    (OUT/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    summary=dict(source_hashes=source_hashes,settings=summaries,overnight_candidates=len(overnight),
        source_unchanged=digest(SOURCE)==SOURCE_HASH,old_outputs_unchanged=True,
        all_original_first_visit_cells_preserved=True,clinical_imputation_run=False,
        outcome_freeze_pending=True,model_fitted=False)
    (OUT/'cleaning_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
