#!/usr/bin/env python3
"""Release v3: explicit Castor metadata, dependency-aware typing, MI inputs.

Patient source remains CSV2026; XLSX2024 supplies versioned metadata only.
Original first-visit cells retained, new analytical companions are authoritative
typed values. Never impute nonperformance, unsupported coding or uncertain timing.
"""
from pathlib import Path
from functools import lru_cache
import importlib.util
import json
import re
import shutil
import hashlib
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/first_visit_2026-10-27_v3'
PREV=ROOT/'results/first_visit_2026-09-30'
V2=ROOT/'results/first_visit_2026-10-27'
spec=importlib.util.spec_from_file_location('v2',ROOT/'analysis/09_prepare_first_visit_v2.py')
v2=importlib.util.module_from_spec(spec);spec.loader.exec_module(v2)
base=v2.base
GENUINE={'asked_unknown','not_described','not_recorded','measurement_failed'}
BLOCKED_PARTS=('timing','after_index','reference_standard','postoperative_30day_window_incomplete','date_missing','surgery_after_index')
UNKNOWN_LABELS={'onbekend','onduidelijk','niet omschreven','niet beschreven','unknown','asked but unknown'}
ELIGIBLE_ROLES={'baseline_or_first_visit_source_field','historical_information_with_timing_gate',
    'same_day_laboratory','same_day_diagnostic_information','same_day_questionnaire'}

def read(p):return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
def write(df,name):df.to_csv(OUT/name,index=False,na_rep='',compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def normunit(x):
    s=x.strip().lower().replace('µ','u').replace('μ','u').replace(' ','')
    return {'kilo':'kg','kilogram':'kg','graden':'degree_C','gradencelcius':'degree_C','gradencelsius':'degree_C',
        'mmol/l':'mmol/L','umol/l':'umol/L','mg/l':'mg/L','u/l':'U/L','beats/min':'beats/min',
        'mmhg':'mmHg','%':'percent','kg/m2':'kg/m2','kg/m²':'kg/m2'}.get(s,s)

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    assert sha(v2.SOURCE)==v2.SOURCE_HASH
    mappings=json.loads((OUT/'codebook/mappings.json').read_text())
    raw=pd.read_csv(v2.SOURCE,sep=';',dtype=str,keep_default_na=False,na_filter=False).apply(lambda s:s.str.strip())
    # Every coded item consulted by the inherited outcome rule has explicit metadata.
    critical=list(v2.old.DX_NEW_COLUMNS)+['FU_einde_diagnose_definitief_ordinaal','FU_einde_diagnose','VG_cholecystectomie']
    critical += [c for row in v2.old.SURGERY_EVIDENCE for c in [row[0],row[2],row[4]]]
    critical += [row[1] for row in v2.old.ECHO_CHOLECYSTITIS_EVIDENCE+v2.old.CT_CHOLECYSTITIS_EVIDENCE]
    critical += [c for row in v2.old.ERCP_EVIDENCE for c in [row[1],row[3]]]
    assert all(mappings[c]['verified_crosswalk'] for c in critical),[c for c in critical if not mappings[c]['verified_crosswalk']]
    for c in v2.old.DX_NEW_COLUMNS:
        assert mappings[c]['codes']['5']=='Cholecystolithiasis' and mappings[c]['codes']['6']=='Cholecystitis'
    for row in v2.old.SURGERY_EVIDENCE:
        assert 'cholecystectomie' in mappings[row[2]]['codes']['2'].lower()
        assert 'cholecystolithiotomie' in mappings[row[2]]['codes']['3'].lower()
    for row in v2.old.ECHO_CHOLECYSTITIS_EVIDENCE+v2.old.CT_CHOLECYSTITIS_EVIDENCE:
        assert mappings[row[1]]['codes']['1'].lower()=='yes'
    groups=base.build_groups(list(raw.columns))
    selectors={g['selector']:dict(g['mapping']) for gs in groups.values() for g in gs}
    totals={}; available=[]; applicability_audit=[]; dictionaries={}
    for setting in ['SEH','POLI']:
        source_clean=read(PREV/f'cohort_{setting}_cleaned.csv.gz')
        source_states=read(PREV/f'cohort_{setting}_analysis_missing_status.csv.gz')
        dictionary=read(PREV/f'variable_dictionary_{setting}.csv')
        dictionary=dictionary[dictionary.output_location.eq('cleaned_cohort')].copy()
        clean=source_clean.copy();states=source_states.copy()
        r=raw.set_index('Participant Id',drop=False).loc[clean.participant_id].reset_index(drop=True)
        outcome=v2.endpoint(raw,setting,clean)
        outcome['gallstone_related_disease']=outcome.candidate_endpoint_after_exclusions
        outcome['outcome_mapping_status']='validated_versioned_castor_metadata_crosswalk'
        outcome['analysis_eligible']=outcome.gallstone_related_disease.ne('').astype(int)
        # Missing-outcome people are kept in a separate unconditioned MI stratum.
        person_ok=~outcome.prior_cholecystectomy_excluded.eq(1)&~outcome.prior_cholecystectomy_history_unknown.eq(1)&~outcome.identity_review_required.eq(1)
        @lru_cache(maxsize=None)
        def applicability(src, trail=()):
            """Three-valued applicability (1 yes,0 no,-1 unresolved)."""
            meta=mappings.get(src,{})
            if not meta.get('verified_crosswalk') or src in trail:return pd.Series(-1,index=r.index,dtype=int)
            if not meta.get('dependency_parent'):return pd.Series(1,index=r.index,dtype=int)
            parent=meta.get('dependency_parent_source','')
            if not parent or parent not in r:return pd.Series(-1,index=r.index,dtype=int)
            parent_app=applicability(parent,trail+(src,))
            values=r[parent];pm=mappings.get(parent,{})
            valid=values.ne('')&~values.isin(v2.base.master.MISSING)
            operator=meta.get('dependency_operator');expected=meta.get('dependency_value','')
            if pm.get('type')=='checkbox':
                # Expanded checkbox dependency is valid only for its own option.
                if expected!=pm.get('selected_option_group_code'):return pd.Series(-1,index=r.index,dtype=int)
                expected='1'
            if operator=='==': hit=values.eq(expected)
            elif operator=='!=': hit=values.ne(expected)
            elif operator in ('>','>=','<','<='):
                a=pd.to_numeric(values,errors='coerce');b=pd.to_numeric(expected,errors='coerce')
                valid &= a.notna() & pd.notna(b)
                hit={'>':a.gt,'>=':a.ge,'<':a.lt,'<=':a.le}[operator](b)
            else:return pd.Series(-1,index=r.index,dtype=int)
            out=pd.Series(-1,index=r.index,dtype=int)
            out.loc[parent_app.eq(0)]=0
            known=parent_app.eq(1)&valid
            out.loc[known]=hit.loc[known].astype(int)
            return out
        columns_by_original=dict(zip(dictionary.original_or_derived_name,dictionary.variable))
        additions={};status_additions={};rows=[];specs=[];input_cols={};where={};elig={}
        numerical_companions={x.companion_of:x for _,x in dictionary[dictionary.numeric_companion.eq('numeric')].iterrows()}
        allblocks=set()
        for order,(_,d) in enumerate(dictionary.iterrows()):
            col=d.variable; original=d.original_or_derived_name
            if col not in clean:continue
            if d.numeric_companion or d.role not in ELIGIBLE_ROLES:continue
            source_options=json.loads(d.source_columns)
            if not source_options:continue
            verified=[s for s in source_options if mappings.get(s,{}).get('verified_crosswalk')]
            if not verified:continue
            representative=mappings[verified[0]]
            typ=representative['type']
            if typ not in {'radio','checkbox','numeric','slider','textarea','string','calculation'}:continue
            # Determine source instance per person, never borrow a later result.
            chosen=pd.Series('',index=r.index,dtype='string')
            if original in r:chosen[:]=original
            else:
                match=re.fullmatch(r'(index_[a-z]+(?:_repeat\d+)?)_(.+)',original)
                if match:
                    prefix,field=match.groups();selector_col=prefix.lower()+'_source_instance'
                    if selector_col in clean:
                        for selector in clean[selector_col].unique():
                            src=selectors.get(selector,{}).get(field,'')
                            if src:chosen.loc[clean[selector_col].eq(selector)]=src
                if len(source_options)==1:chosen.loc[chosen.eq('')]=source_options[0]
            verified_row=chosen.isin(verified)
            app=pd.Series(-1,index=r.index,dtype=int)
            for src in chosen[verified_row].unique():
                use=chosen.eq(src);app.loc[use]=applicability(src).loc[use]
            # Require common type and code meaning before combining repeated sources.
            signature=lambda m:(m['type'],json.dumps(m.get('codes',{}),sort_keys=True))
            if any(signature(mappings[s])!=signature(representative) for s in verified):continue
            name=col+'_analysis';value=clean[col].copy();state=states[col].copy()
            unit=normunit(representative.get('units',''))
            reason='';method='';levels=[];formula=''
            block='form_'+representative['form_number']
            allblocks.add(block)
            temporal=state.str.contains('|'.join(BLOCKED_PARTS),regex=True)
            safe_app=app.eq(1)&verified_row&~temporal
            # Always retain source cells. Analytical companions withhold hidden branches.
            state.loc[app.eq(0)&~temporal]='not_applicable_dependency'
            state.loc[app.eq(-1)&~temporal]='dependency_unresolved'
            value.loc[~safe_app]=''
            unknown_codes=[]
            if typ=='calculation':
                reason='Derived quantity: recompute from components, never impute independently'
                continue
            elif typ in {'textarea','string'}:
                # Free text is retained; exact NO only recoded where CRF says so.
                if 'niet omschreven' in representative.get('info','').lower():
                    missing=value.str.strip().str.upper().eq('NO')
                    value.loc[missing]='';state.loc[missing]='not_described'
                else:continue
                reason='Narrative text; never stochastically imputed'
            elif typ=='radio':
                codes=representative['codes']
                unknown_codes=[k for k,v in codes.items() if v.strip().lower() in UNKNOWN_LABELS]
                valid={k:v for k,v in codes.items() if k not in unknown_codes}
                missing=clean[col].isin(unknown_codes)&safe_app
                state.loc[missing]='asked_unknown';value.loc[missing]=''
                bad=value.ne('')&~value.isin(valid)
                state.loc[bad]='category_code_unresolved';value.loc[bad]=''
                levels=sorted(valid)
                method='logreg' if len(levels)==2 else 'polyreg' if len(levels)>2 else ''
                reason='Explicit Castor option codes; missing-answer categories separately recorded'
                additions[col+'_label']=clean[col].map(codes).fillna('')
                status_additions[col+'_label']=states[col]
                rows.append(dict(variable=col+'_label',source_variable=col,source_names=json.dumps(source_options),
                    type='label',role=d.role,method='',note='Source labels; not imputed'))
            elif typ=='checkbox':
                # 0 is not a negative clinical observation if the whole group is unfilled.
                srcbase=verified[0].split('#')[0]
                siblings=[s for s in r if s.split('#')[0]==srcbase and '#' in s]
                answered=r[siblings].eq('1').any(axis=1) if siblings else pd.Series(False,index=r.index)
                unanswered=safe_app&~answered
                value.loc[unanswered]='';state.loc[unanswered]='checkbox_group_unanswered'
                bad=value.ne('')&~value.isin(['0','1']);value.loc[bad]='';state.loc[bad]='category_code_unresolved'
                levels=['0','1'];reason='Multiselect structure retained; no independent option imputations'
            else:
                # Preserve exact observed numbers, inequalities and parse/range issues separately.
                if col in numerical_companions:
                    prior=numerical_companions[col]
                    value=clean[prior.variable].copy();state=states[prior.variable].copy()
                    if not unit and not any(t in prior.unit for t in ['not_specified','unresolved']):unit=normunit(prior.unit)
                    value.loc[~safe_app]=''
                    state.loc[app.eq(0)&~temporal]='not_applicable_dependency'
                    state.loc[app.eq(-1)&~temporal]='dependency_unresolved'
                else:
                    parsed=[];newstate=[]
                    for v,st in zip(value,state):
                        if not v:parsed.append('');newstate.append(st);continue
                        m=base.NUM.fullmatch(v.strip())
                        if not m:parsed.append('');newstate.append('numeric_parse_unresolved');continue
                        if m.group(1):parsed.append('');newstate.append('censored_measurement_bound_retained');continue
                        number=float(m.group(2).replace(',','.'))
                        bounds=[pd.to_numeric(representative.get(k,''),errors='coerce') for k in ['minimum','maximum']]
                        invalid=(pd.notna(bounds[0]) and number<bounds[0]) or (pd.notna(bounds[1]) and number>bounds[1])
                        parsed.append('' if invalid else str(number));newstate.append('numeric_plausibility_review' if invalid else 'numeric_observed')
                    value=pd.Series(parsed,index=r.index);state=pd.Series(newstate,index=r.index)
                method='pmm';reason='Exact numeric measurement; censored/unresolved values not filled'
            # Blanks in applicable mandatory clinical questions can be missing.
            # For tests, performance of a panel does not prove each analyte was ordered.
            blank_allowed=representative.get('required')=='1' and d.role in {'baseline_or_first_visit_source_field','historical_information_with_timing_gate'}
            if d.role=='same_day_laboratory':
                blank_allowed=False
                if not unit or 'unresolved' in unit:method='';reason='Laboratory unit unresolved; retain measured source, no imputation'
            if d.role in {'same_day_diagnostic_information','same_day_questionnaire'}:blank_allowed=False
            if typ in {'checkbox','textarea','string'}:method=''
            if re.search(r'_repeat\d+_',original):method='';reason='Repeat result retained; primary initial measurement only in MI'
            reasons=set(GENUINE)
            if blank_allowed:reasons.add('blank_or_dependency_not_triggered')
            if original=='B_ASA':
                # Metadata declares numeric, but the classification contains values
                # such as 163/173. No invented correction or continuous PMM model.
                # Preserve source cells; withhold this analytical representation
                # both as an imputation target and as an auxiliary until adjudicated.
                value[:]='';state[:]='classification_coding_review_required'
                safe_app[:]=False;method=''
                reason='ASA classification quarantined: implausible source codes; source retained, no target or auxiliary use'
            if re.fullmatch(r'index_urine(?:_repeat\d+)?_sg',col) or ('egfr' in col and unit=='kl/1.73m2'):
                value[:]='';state[:]='measurement_scale_or_unit_review_required'
                safe_app[:]=False;method=''
                reason='Unresolved measurement scale/unit; source retained without guessed conversion, no target or auxiliary use'
            if re.fullmatch(r'index_urine(?:_repeat\d+)?_ph',col):
                ph=pd.to_numeric(value,errors='coerce')
                invalid=ph.lt(0)|ph.gt(14)
                value.loc[invalid]='';state.loc[invalid]='numeric_plausibility_review'
                safe_app.loc[invalid]=False
            missing=value.eq('')&state.isin(reasons)&safe_app&person_ok
            observed=value.ne('')&safe_app&person_ok
            # Selection-defining history never imputed to manufacture eligibility.
            if original=='VG_cholecystectomie':method='';reason='Defines exclusion; history not imputed'
            additions[name]=value;status_additions[name]=state
            specrow=dict(variable=name,method=method,block=block,levels_json=json.dumps(levels) if levels else '',
                reason=reason,verification_state='current_machine_dictionary_verified',auxiliary_rank=order,
                auxiliary_for_blocks=json.dumps([block]),source_metadata_version='2024 Castor metadata crosschecked against 2026 CSV/CRF')
            specs.append(specrow)
            rows.append(dict(variable=name,source_variable=col,source_names=json.dumps(source_options),
                type=typ,role=d.role,unit=unit,method=method,note=reason,question=representative['question'],
                form=representative['form'],option_codes_json=representative['codes_json'],
                dependency_parent=representative['dependency_parent'],dependency_operator=representative['dependency_operator'],
                dependency_value=representative['dependency_value'],metadata_status=representative['verification_state']))
            if typ not in {'textarea','string'}:
                input_cols[name]=value;where[name]=missing if method else pd.Series(False,index=r.index)
                elig[name]=(observed|missing) if method else observed
            applicability_audit.append(dict(setting=setting,variable=name,applicable_n=int(app.eq(1).sum()),
                not_applicable_n=int(app.eq(0).sum()),unresolved_n=int(app.eq(-1).sum()),
                candidate_missing_n=int(missing.sum()),method=method,reason=reason))
        clean=pd.concat([clean,pd.DataFrame(additions)],axis=1)
        states=pd.concat([states,pd.DataFrame(status_additions)],axis=1)
        assert not clean.columns.duplicated().any()
        for col in ['prior_cholecystectomy_excluded','identity_review_required','analysis_eligible']:
            clean[col]=outcome[col].astype(str);states[col]='administrative_or_selection'
        # One consistent derived BMI. Reported BMI and numerical companion retained.
        if {'startgewicht_analysis','lengte_analysis'} <= set(clean):
            clean['baseline_bmi_calculated']=pd.to_numeric(clean.startgewicht_analysis,errors='coerce')/(pd.to_numeric(clean.lengte_analysis,errors='coerce')/100)**2
            states['baseline_bmi_calculated']=np.where(clean.baseline_bmi_calculated.notna(),'derived_from_observed_components','components_missing_not_independently_imputed')
        input_frame=pd.DataFrame({'participant_id':clean.participant_id,'gallstone_related_disease':outcome.gallstone_related_disease,
            'age_years_at_index':clean.age_years_at_index,'years_since_bariatric_surgery':clean.years_since_bariatric_surgery,
            'geslacht':clean.get('geslacht_analysis',clean.geslacht),**input_cols})
        allspec=pd.DataFrame(specs)
        # Explicitly assess every retained field, including source/administrative representations.
        represented=set(allspec.variable)
        remainder=[dict(variable=c,method='',block='',levels_json='',reason='Source, derived, administrative, narrative or unverified field; retained without stochastic imputation',verification_state='not_eligible',auxiliary_rank=99999,auxiliary_for_blocks='[]') for c in clean if c not in represented]
        allspec=pd.concat([allspec,pd.DataFrame(remainder)],ignore_index=True)
        dictionary['metadata_type']='';dictionary['metadata_status']=''
        for i,d in dictionary.iterrows():
            ss=json.loads(d.source_columns);mm=[mappings[s] for s in ss if s in mappings]
            if mm:
                dictionary.at[i,'metadata_type']=mm[0].get('type','')
                dictionary.at[i,'metadata_status']=mm[0].get('verification_state','')
        enriched=pd.DataFrame(rows);enriched['setting']=setting
        for col in clean:
            if col in ['participant_id','analysis_setting']:continue
            values=clean[col].fillna('').astype(str); st=states[col]
            available.append(dict(setting=setting,variable=col,n=len(clean),observed_n=int(values.ne('').sum()),
                observed_pct=round(100*values.ne('').mean(),3),missing_status_counts=json.dumps(st[values.eq('')].value_counts().to_dict()),
                eligible_missing_n=int(where[col].sum()) if col in where else 0))
        write(clean,f'cohort_{setting}_cleaned.csv.gz');write(clean,f'cohort_{setting}_cleaned.csv')
        write(states,f'cohort_{setting}_analysis_missing_status.csv.gz')
        write(outcome,f'outcomes_{setting}.csv')
        write(dictionary,f'variable_dictionary_source_{setting}.csv');write(enriched,f'variable_dictionary_analysis_{setting}.csv')
        write(input_frame,f'imputation_input_{setting}.csv.gz');write(allspec,f'imputation_spec_{setting}.csv')
        write(pd.DataFrame({'participant_id':clean.participant_id,**where}),f'where_{setting}.csv.gz')
        write(pd.DataFrame({'participant_id':clean.participant_id,**elig}),f'eligibility_{setting}.csv.gz')
        assert source_clean.equals(clean[source_clean.columns])
        assert clean.participant_id.is_unique and list(clean)==list(states)
        for c in where:assert (~where[c]|(input_frame[c].eq('')&elig[c])).all()
        totals[setting]=dict(first_presentations=len(clean),full_clean_columns=len(clean.columns),
            confirmed_events=int(outcome.gallstone_related_disease.eq('1').sum()),nonevents=int(outcome.gallstone_related_disease.eq('0').sum()),
            prior_chole_excluded=int(outcome.prior_cholecystectomy_excluded.sum()),duplicate_identity_pending=int(outcome.identity_review_required.sum()),
            outcome_unresolved_or_excluded=int(outcome.gallstone_related_disease.eq('').sum()),
            fields_assessed=len(allspec),typed_analytical_companions=len(rows),
            fields_with_eligible_missing=sum(int(x.any()) for x in where.values()),eligible_missing_cells=sum(int(x.sum()) for x in where.values()))
    write(pd.DataFrame(available),'availability_all_variables.csv')
    write(pd.DataFrame(applicability_audit),'applicability_audit.csv')
    for name in ['first_episodes_unassigned_setting.csv.gz','cohort_flow.csv','quality_flags.csv','timing_summary.csv']:
        shutil.copy2(PREV/name,OUT/name)
    shutil.copy2(V2/'overnight_visit_linkage_review.csv',OUT/'overnight_visit_linkage_review.csv')
    cfg=dict(version='2026-10-27_v3',outcome_mapping_status='validated',clinical_imputation_run=False,m=50,maxit=30,maxit_if_unstable=100,
        source_patient_data=str(v2.SOURCE.relative_to(ROOT)),metadata_version='2024-06-10 checked against CSV2026 and CRF2026',
        prior_cholecystectomy_policy='exclude_both_cases_controls_retain_fullclean',
        unresolved_identity_policy='retain_fullclean_exclude_MI_and_primary_analysis_until_identity_review',
        lab_missingness_policy='only explicitly missing/failed analytes; blank assays not assumed ordered',
        predictor_scope='all first-visit fields assessed; no outcome association screening; no prediction model fitted')
    (OUT/'configuration.json').write_text(json.dumps(cfg,indent=2)+'\n')
    result=dict(settings=totals,source_sha256=sha(v2.SOURCE),source_unchanged=sha(v2.SOURCE)==v2.SOURCE_HASH,
        original_first_visit_cells_preserved=True,outcome_metadata_fields_verified=len(set(critical)),
        all_patients_retained=True,imputation_run=False)
    (OUT/'cleaning_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
