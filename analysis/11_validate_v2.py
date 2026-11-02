#!/usr/bin/env python3
"""Independent read-only validation of the 2026-10-27 cleaning release.
Writes only statistical_validation.json. Synthetic endpoint tests do not rerun
cleaning or MICE and do not establish validity of unverified source mappings.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/first_visit_2026-10-27'
PREV=ROOT/'results/first_visit_2026-09-30'

def read(p): return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    checks=[]
    def check(name,condition):
        checks.append({'check':name,'passed':bool(condition)})
    spec=importlib.util.spec_from_file_location('operate_v2_review',ROOT/'analysis/09_prepare_first_visit_v2.py')
    prep=importlib.util.module_from_spec(spec); spec.loader.exec_module(prep)
    check('source_sha256_unchanged',digest(prep.SOURCE)==prep.SOURCE_HASH)
    receipt=json.loads((PREV/'first_visit_metadata.json').read_text())
    check('previous_release_hashes_unchanged',all(digest(PREV/n)==i['sha256'] for n,i in receipt['outputs'].items()))
    mappings=json.loads((OUT/'codebook/reviewed_mappings.json').read_text())
    mapping_by_source={x['source_name']:x for x in mappings}
    for setting in ['SEH','POLI']:
        original=read(PREV/f'cohort_{setting}_cleaned.csv.gz')
        clean=read(OUT/f'cohort_{setting}_cleaned.csv.gz')
        states=read(OUT/f'cohort_{setting}_analysis_missing_status.csv.gz')
        inp=read(OUT/f'imputation_input_{setting}.csv.gz')
        where=read(OUT/f'where_{setting}.csv.gz')
        eligible=read(OUT/f'eligibility_{setting}.csv.gz')
        outcome=read(OUT/f'outcomes_{setting}_review.csv')
        dictionary=read(OUT/f'variable_dictionary_{setting}.csv')
        check(setting+'_all_original_columns_and_values_preserved',set(original)<=set(clean) and original.equals(clean[list(original)]))
        check(setting+'_patient_identity_order_preserved',all(clean.participant_id.equals(x.participant_id) for x in [states,inp,where,eligible,outcome]))
        check(setting+'_unique_patient_record_ids',clean.participant_id.is_unique)
        check(setting+'_status_alignment',list(states)==list(clean))
        check(setting+'_final_outcomes_unfrozen_unimputed',outcome.gallstone_related_disease.eq('').all() and inp.gallstone_related_disease.eq('').all() and outcome.analysis_eligible.eq('0').all())
        prior=outcome.prior_cholecystectomy_excluded.eq('1')
        duplicate=outcome.identity_review_required.eq('1')
        check(setting+'_prior_chole_both_groups_excluded',outcome.loc[prior,'candidate_endpoint_after_exclusions'].eq('').all() and outcome.loc[prior,'analysis_eligible'].eq('0').all())
        check(setting+'_duplicate_pending_excluded',outcome.loc[duplicate,'candidate_endpoint_after_exclusions'].eq('').all() and outcome.loc[duplicate,'analysis_eligible'].eq('0').all())
        for col in where:
            if col=='participant_id': continue
            w=where[col].str.lower().isin(['true','1'])
            e=eligible[col].str.lower().isin(['true','1'])
            check(setting+'_'+col+'_where_empty_eligible',(~w | (e & inp[col].eq(''))).all())
            check(setting+'_'+col+'_no_excluded_person_imputation',(~w | ~(prior|duplicate|outcome.prior_cholecystectomy_history_unknown.eq('1'))).all())
            allowed={'asked_unknown','not_described','not_recorded','measurement_failed'}
            if not col.startswith('index_lab') and '_lo_' not in col:
                allowed.add('blank_or_dependency_not_triggered')
            check(setting+'_'+col+'_no_structural_imputation',states.loc[w,col].isin(allowed).all())
            if col.startswith('index_lab'):
                check(setting+'_'+col+'_blank_analytes_not_imputed',~(w & states[col].eq('blank_or_dependency_not_triggered')).any())
        for _,row in dictionary[dictionary.storage.eq('categorical_analysis_companion')].iterrows():
            col=row.variable; src=row.companion_of
            mapping=mapping_by_source[row.original_or_derived_name]
            codes=json.loads(mapping['codes_json'])
            unknown={k for k,v in codes.items() if v.strip().lower() in {'onbekend','onduidelijk','niet omschreven'}}
            valid={k for k in codes if not k.startswith('-')}-unknown
            expected=clean[src].where(clean[src].isin(valid),'')
            check(setting+'_'+col+'_source_code_not_reordered',clean[col].equals(expected))
        for suffix in ['before_timing_filter','source_missing_status']:
            check(setting+'_'+suffix+'_snapshot_unchanged',digest(OUT/f'cohort_{setting}_{suffix}.csv.gz')==digest(PREV/f'cohort_{setting}_{suffix}.csv.gz'))
    overnight=read(OUT/'overnight_visit_linkage_review.csv')
    check('overnight_candidates_not_assumed_index_visit',overnight.status.eq('next_day_candidate_no_verified_visit_end_or_link').all())
    config=json.loads((OUT/'configuration.json').read_text())
    check('clinical_imputation_blocked_pending_current_mapping',config['clinical_imputation_run'] is False and config['outcome_mapping_status']=='pending_current_code_verification')
    check('no_imputed_clinical_outputs_published',not list(OUT.glob('imputed_*.csv.gz')))
    # Synthetic cases exercise the actual endpoint implementation. Header only;
    # no participant source values are copied into these constructed records.
    names=list(pd.read_csv(prep.SOURCE,sep=';',nrows=0).columns)
    cases=[('at_index','2020-01-01','1','1','5'),('at_end','2020-06-01','1','1','5'),
           ('after_end','2020-06-02','1','1','5'),('before_index','2019-12-31','1','1','5'),
           ('prior_positive','2020-01-01','1','3','5'),('prior_negative','','0','3',''),
           ('unclear_negative','','0','1',''),('uncertain_text_negative','','0','1',''),('no_assessment','','','1',''),
           ('label_unverified','','1','1','5'),('procedure_without_label','2020-01-01','1','1','')]
    rows=[]
    surgery=prep.old.SURGERY_EVIDENCE[0]
    for name,date,assessed,prior,dx in cases:
        row=dict.fromkeys(names,''); row.update({'Participant Id':name,'SEH_datum':'2020-01-01','FU_einde_datum':'2020-06-01','FU_einde_diagnose':assessed,'VG_cholecystectomie':prior,prep.old.DX_NEW_COLUMNS[0]:dx})
        if name=='uncertain_text_negative': row[prep.old.DX_TEXT_COLUMNS[0]]='mogelijk cholecystolithiasis'
        if date: row.update({surgery[0]:'1',surgery[1]:date,surgery[2]:'2'})
        rows.append(row)
    raw=pd.DataFrame(rows)
    clean=pd.DataFrame({'participant_id':raw['Participant Id'],'duplicate_review_required':'0'})
    result=prep.endpoint(raw,'SEH',clean).set_index('participant_id')
    for name in ['at_index','at_end']:
        check('synthetic_'+name+'_included',result.loc[name,'candidate_endpoint_after_exclusions']=='1')
    for name in ['after_end','before_index','label_unverified']:
        check('synthetic_'+name+'_unverified_not_zero',result.loc[name,'candidate_endpoint_before_exclusions']=='')
    for name in ['prior_positive','prior_negative']:
        check('synthetic_'+name+'_excluded',result.loc[name,'candidate_endpoint_after_exclusions']=='' and result.loc[name,'prior_cholecystectomy_excluded']==1)
    check('synthetic_unclear_final_diagnosis_zero',result.loc['unclear_negative','candidate_endpoint_after_exclusions']=='0')
    check('synthetic_uncertain_text_assessed_no_definitive_diagnosis_zero',result.loc['uncertain_text_negative','candidate_endpoint_after_exclusions']=='0')
    check('synthetic_absent_assessment_missing',result.loc['no_assessment','candidate_endpoint_before_exclusions']=='')
    check('synthetic_procedure_alone_not_positive',result.loc['procedure_without_label','candidate_endpoint_before_exclusions']!='1')
    payload={'review_role':'Separate Statistical Reviewer','scope':'Independent release invariants and synthetic endpoint tests; no clinical MICE/model execution','implementation_sha256':digest(ROOT/'analysis/09_prepare_first_visit_v2.py'),'passed':all(x['passed'] for x in checks),'checks_n':len(checks),'failed_checks':[x for x in checks if not x['passed']],'checks':checks,'limitations':['Current-release numerical category/procedure mappings remain unverified.','Synthetic checks validate the provisional rule implementation, not source-code meanings.','No outcome freeze or statistical validity certification.']}
    (OUT/'statistical_validation.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({k:payload[k] for k in ['passed','checks_n','failed_checks']},indent=2))
    return payload['passed']
if __name__=='__main__': raise SystemExit(0 if run() else 1)
