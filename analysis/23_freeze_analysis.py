#!/usr/bin/env python3
"""Freeze saved v3 no-adjudication person eligibility and diagnostic timing flags.

No refitting, imputation, correction, endpoint reclassification or new exclusions.
Same-day performed tests support an index-workup description, never proof that
results were available before assessment. Strict presentation models must omit
unproven diagnostic workup predictors; assess omission as a sensitivity for an
explicitly defined workup model. Source v3 files remain unchanged.
"""
from pathlib import Path
import hashlib
import json
from collections import Counter
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V3 = ROOT / 'results/first_visit_2026-10-27_v3'
OUT = ROOT / 'results/analysis_2026-10-27_no_adjudication'
PAIR = ROOT / 'results/cleaned_2026-09-30/candidate_duplicate_first_episodes.csv'
GOOD_LINK = {'matched','multiple_same_date_records'}
MODALITIES = ['lab','urine','feces','xboz','ct','echo','xslik','gastro','colon','mrcp','ercp','mri']


def read(path):
    return pd.read_csv(path,dtype=str,keep_default_na=False)


def main():
    OUT.mkdir(exist_ok=True)
    overnight = read(V3/'overnight_visit_linkage_review.csv')
    clinical = read(V3/'clinical_review_queue.csv')
    pairs = read(PAIR)
    duplicate_ids = set(pairs.participant_id_a)|set(pairs.participant_id_b)
    clinical_ids = set(clinical.participant_id)
    overnight_ids = set(overnight.participant_id)
    records,timing,summary,inputs = [],[],{},[PAIR,V3/'overnight_visit_linkage_review.csv',V3/'clinical_review_queue.csv']
    modeldata = {}
    for setting in ('SEH','POLI'):
        paths = [V3/f'outcomes_{setting}.csv',V3/f'cohort_{setting}_cleaned.csv.gz',
                 V3/f'cohort_{setting}_analysis_missing_status.csv.gz',V3/f'variable_dictionary_analysis_{setting}.csv']
        outcome,clean,states,dictionary = map(read,paths)
        inputs.extend(paths)
        assert outcome.participant_id.is_unique and clean.participant_id.is_unique
        clean=clean.set_index('participant_id').loc[outcome.participant_id]
        states=states.set_index('participant_id').loc[outcome.participant_id]
        modeldata[setting]=(clean,states)
        reliable_lab_fields = dictionary.loc[
            dictionary.role.eq('same_day_laboratory') & dictionary.type.isin(['numeric','calculation']) &
            dictionary.unit.ne('') & ~dictionary.unit.str.contains('unresolved|not_specified|kl/1.73',case=False,regex=True),
            'variable'].tolist()
        for _,r in outcome.iterrows():
            pid=r.participant_id;eligible=r.analysis_eligible=='1'
            if eligible != (r.gallstone_related_disease in ('0','1')):
                raise ValueError('Saved eligibility and final endpoint disagree')
            if eligible and (pid in duplicate_ids or r.prior_cholecystectomy_excluded=='1'):
                raise ValueError('Saved primary cohort contains unresolved identity or prior-chole exclusion')
            flags=[]
            if r.review_reason:flags.append(r.review_reason)
            if pid in duplicate_ids:flags.append('possible_duplicate_person_identity_unresolved')
            if pid in overnight_ids:flags.append('next_day_test_visit_linkage_unverified')
            if pid in clinical_ids:flags.append('clinical_value_or_scale_review_source_retained_fields_withheld')
            records.append(dict(participant_id=pid,setting=setting,finalY=r.gallstone_related_disease,
                                primaryeligibility=int(eligible),exclusionreason='' if eligible else r.review_reason,
                                existingopenflags=json.dumps(sorted(set(flags)))))
            c=clean.loc[pid];st=states.loc[pid]
            for modality in MODALITIES:
                prefix='index_'+modality
                link=c.get(prefix+'_link_status','not_recorded')
                n=pd.to_numeric(c.get(prefix+'_same_day_record_count',''),errors='coerce')
                date_performed=link in GOOD_LINK
                observed=0
                if modality=='lab':
                    observed=sum(c.get(field,'')!='' and st.get(field,'')=='numeric_observed'
                                 for field in reliable_lab_fields)
                timing.append(dict(participant_id=pid,setting=setting,primaryeligibility=int(eligible),
                                   modality=modality,existing_link_status=link,
                                   same_day_record_count=int(n) if pd.notna(n) else 0,
                                   same_day_plus_performance_supported=int(date_performed),
                                   test_time_recorded=int(bool(c.get(prefix+'_time',''))),
                                   verified_unit_observed_lab_cells=observed if modality=='lab' else '',
                                   reliable_index_workup_lab_available=int(date_performed and observed>0) if modality=='lab' else '',
                                   availability_before_assessment_proven=0,
                                   strict_presentation_predictor_use=0,
                                   workup_use_requires_explicit_landmark=1))
        primary=outcome.analysis_eligible.eq('1')
        summary[setting]=dict(first_visit_records=len(outcome),unique_ids=outcome.participant_id.nunique(),
                              primary_eligible=int(primary.sum()),events=int(outcome.gallstone_related_disease.eq('1').sum()),
                              nonevents=int(outcome.gallstone_related_disease.eq('0').sum()),excluded=int((~primary).sum()),
                              exclusion_reasons=outcome.loc[~primary,'review_reason'].value_counts().to_dict(),
                              duplicate_flag_records=int(outcome.participant_id.isin(duplicate_ids).sum()))
    unassigned_path=V3/'first_episodes_unassigned_setting.csv.gz'
    unassigned=read(unassigned_path);inputs.append(unassigned_path)
    for pid in unassigned.participant_id:
        flags=['first_presentation_setting_unresolved']
        if pid in duplicate_ids:flags.append('possible_duplicate_person_identity_unresolved')
        records.append(dict(participant_id=pid,setting='unassigned',finalY='',primaryeligibility=0,
                            exclusionreason='first_presentation_setting_unresolved',existingopenflags=json.dumps(flags)))
    cohort=pd.DataFrame(records);time=pd.DataFrame(timing)
    assert len(cohort)==1593 and cohort.participant_id.is_unique
    assert int(cohort.primaryeligibility.sum())==1130
    assert int(cohort.participant_id.isin(duplicate_ids).sum())==8
    assert not cohort.loc[cohort.participant_id.isin(duplicate_ids),'primaryeligibility'].any()
    spec_path=OUT/'model_specification.json'
    if not spec_path.exists():
        raise ValueError('Reviewed model_specification.json is required before freezing model populations')
    spec=json.loads(spec_path.read_text());inputs.append(spec_path)
    # Comparators intentionally share the full model's population. Evaluate the
    # richest predictor set for each population flag, rather than broadening its
    # comparator population when a predictor is omitted from the comparator fit.
    flag_specs={}
    for model in spec['models']:
        flag=model['population_flag']
        if flag not in flag_specs or len(model['predictors'])>len(flag_specs[flag]['predictors']):
            flag_specs[flag]=model
    for flag,model in flag_specs.items():
        clean,states=modeldata[model['setting']]
        eligible=cohort.set_index('participant_id').loc[clean.index,'primaryeligibility'].eq(1)
        for predictor in model['predictors']:
            source=predictor['source']
            if isinstance(source,list):
                available=pd.Series(False,index=clean.index)
                for field in source:
                    observed=clean[field].isin(predictor.get('levels',['0','1']))
                    gate=predictor.get('source_link_gates',{}).get(field)
                    if gate:
                        observed &= clean[gate['column']].isin(gate['allowed'])
                    available |= observed
            else:
                values=clean[source]
                if predictor.get('recode'):
                    available=values.isin(predictor['recode'])
                elif predictor['type']=='binary':
                    available=values.isin(predictor['levels'])
                else:
                    numbers=pd.to_numeric(values,errors='coerce')
                    available=numbers.notna() & ~numbers.isin([float('inf'),float('-inf')])
                if predictor.get('method'):
                    available |= values.eq('') & states[source].isin(predictor.get('eligible_missing_status',[]))
            eligible &= available
        for field,allowed in model.get('required_link_status',{}).items():
            eligible &= clean[field].isin(allowed)
        admitted=set(clean.index[eligible])
        cohort[flag]=cohort.participant_id.isin(admitted).astype(int)
    model_counts=[]
    for model in spec['models']:
        mask=cohort[model['population_flag']].eq(1)
        n=int(mask.sum());events=int((mask&cohort.finalY.eq('1')).sum())
        clean,_=modeldata[model['setting']]
        time_source=next(p['source'] for p in model['predictors'] if p['name']=='years_since_surgery')
        missing=int(clean.loc[cohort.loc[mask,'participant_id'],time_source].eq('').sum())
        if (n,events,missing)!=(model['expected_n'],model['expected_events'],model['expected_imputable_time_cells']):
            raise ValueError(f"Model population mismatch {model['id']}: observed{n}/{events}/{missing}")
        row=dict(model_id=model['id'],population_flag=model['population_flag'],setting=model['setting'],
                 n=n,events=events,nonevents=n-events,imputable_time_cells=missing)
        model_counts.append(row)
        cohort.loc[mask,['participant_id','setting','finalY',model['population_flag']]].to_csv(
            OUT/f"model_population_{model['id']}.csv",index=False)
    cohort[['participant_id','setting','primaryeligibility',*flag_specs]].to_csv(OUT/'model_population_flags.csv',index=False)
    cohort.to_csv(OUT/'frozen_analysis_cohort.csv',index=False)
    time.to_csv(OUT/'index_diagnostic_timing_flags.csv',index=False)
    summary['unassigned']=dict(first_visit_records=len(unassigned),unique_ids=unassigned.participant_id.nunique(),
                               primary_eligible=0,duplicate_flag_records=int(unassigned.participant_id.isin(duplicate_ids).sum()))
    primary_time=time[time.primaryeligibility.eq(1)]
    timing_counts=[]
    for (setting,modality),group in primary_time.groupby(['setting','modality']):
        timing_counts.append(dict(setting=setting,modality=modality,n_primary=len(group),
                                  same_day_performed_supported=int(group.same_day_plus_performance_supported.sum()),
                                  link_status_counts=group.existing_link_status.value_counts().to_dict(),
                                  reliable_index_workup_lab_available=int(pd.to_numeric(group.reliable_index_workup_lab_available,errors='coerce').fillna(0).sum()),
                                  test_time_recorded=int(group.test_time_recorded.sum()),
                                  before_assessment_availability_proven=0))
    result=dict(person_freeze='Saved v3 endpoint/eligibility only; no new adjudication or person exclusions',
                settings=summary,total_first_episodes=len(cohort),total_primary=1130,
                duplicate_pairs=4,duplicate_ids_excluded=8,
                counting_limit='Unique study IDs; unconfirmed identity pairs are excluded, not presumed different people.',
                primary_diagnostic_timing=timing_counts,
                model_population_counts=model_counts,
                timing_policy=['Same date plus performed flag supports index-workup association only.',
                               'Neither test time nor same-day date proves availability before clinical assessment.',
                               'Strict presentation prediction excludes diagnostic workup predictors with unproven availability.',
                               'An index-workup model needs an explicit landmark and sensitivity omitting labs/imaging.',
                               'Next-day unverified encounters remain withheld; no colleague adjudication assumed.'],
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'frozen_analysis_cohort_summary.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['input_sha256','script_sha256','primary_diagnostic_timing','timing_policy']}))


if __name__=='__main__':
    main()
