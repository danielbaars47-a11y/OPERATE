#!/usr/bin/env python3
"""Separate statistical release audit; no cleaning, imputation or model reruns.
Checks source protection, analytical applicability and completed MI draws.
Run with analysis/.venv/bin/python; writes statistical_validation.json only.
"""
from pathlib import Path
import hashlib,importlib.util,json,csv,gzip
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/first_visit_2026-10-27_v3'
PREV=ROOT/'results/first_visit_2026-09-30'
def read(p):return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def equal(a,b):
    exact=a.eq(b); an=pd.to_numeric(a,errors='coerce');bn=pd.to_numeric(b,errors='coerce')
    return exact|(an.notna()&bn.notna()&np.isclose(an,bn,rtol=1e-12,atol=1e-12))
def main():
    checks=[];pending=[];summaries={};distribution=[]
    def ck(n,x):checks.append({'check':n,'passed':bool(x)})
    sp=importlib.util.spec_from_file_location('prep_v3_audit',ROOT/'analysis/14_prepare_first_visit_v3.py');p=importlib.util.module_from_spec(sp);sp.loader.exec_module(p)
    ck('patient_source_unchanged',sha(p.v2.SOURCE)==p.v2.SOURCE_HASH)
    receipt=json.loads((PREV/'first_visit_metadata.json').read_text())
    ck('previous_release_unchanged',all(sha(PREV/n)==d['sha256'] for n,d in receipt['outputs'].items()))
    mappings=json.loads((OUT/'codebook/mappings.json').read_text())
    raw=read(p.v2.SOURCE) if False else pd.read_csv(p.v2.SOURCE,sep=';',dtype=str,keep_default_na=False).set_index('Participant Id')
    for setting in ['SEH','POLI']:
        old=read(PREV/f'cohort_{setting}_cleaned.csv.gz');clean=read(OUT/f'cohort_{setting}_cleaned.csv.gz');states=read(OUT/f'cohort_{setting}_analysis_missing_status.csv.gz')
        inp=read(OUT/f'imputation_input_{setting}.csv.gz');w=read(OUT/f'where_{setting}.csv.gz');e=read(OUT/f'eligibility_{setting}.csv.gz');y=read(OUT/f'outcomes_{setting}.csv')
        d=read(OUT/f'variable_dictionary_analysis_{setting}.csv');spec=read(OUT/f'imputation_spec_{setting}.csv').set_index('variable')
        ck(setting+'_original_cells_preserved',old.equals(clean[list(old)]))
        ck(setting+'_all_records_aligned',all(clean.participant_id.equals(x.participant_id) for x in [states,inp,w,e,y]))
        ck(setting+'_unique_record_ids',clean.participant_id.is_unique)
        ck(setting+'_all_clean_columns_have_status',list(clean)==list(states))
        excluded=y.prior_cholecystectomy_excluded.eq('1')|y.prior_cholecystectomy_history_unknown.eq('1')|y.identity_review_required.eq('1')
        ck(setting+'_excluded_cases_controls_unlabelled',y.loc[excluded,'gallstone_related_disease'].eq('').all()&y.loc[excluded,'analysis_eligible'].eq('0').all())
        ck(setting+'_outcome_input_agrees',inp.gallstone_related_disease.equals(y.gallstone_related_disease))
        ck(setting+'_outcome_not_target','gallstone_related_disease' not in w)
        if 'b_asa_analysis' in inp:
            asa=pd.to_numeric(inp.b_asa_analysis,errors='coerce')
            ck(setting+'_ASA_analytical_quarantine',inp.b_asa_analysis.eq('').all() and spec.loc['b_asa_analysis','method']=='')
        quarantined=[c for c in inp if c=='b_asa_analysis' or c=='index_urine_sg_analysis' or c.startswith('index_lab_egfr')]
        for c in quarantined:
            ck(setting+'_'+c+'_quarantined_no_analytic_use',inp[c].eq('').all() and spec.loc[c,'method']=='')
        if 'index_urine_ph_analysis' in inp:
            ph=pd.to_numeric(inp.index_urine_ph_analysis,errors='coerce').dropna()
            ck(setting+'_urine_pH_analytical_domain',ph.between(0,14).all())
        predictor_path=OUT/f'imputation_predictor_audit_{setting}.csv'
        if predictor_path.exists():
            pa=read(predictor_path)
            ck(setting+'_quarantined_fields_not_predictors',not (pa.predictor.isin(quarantined)&pa.status.eq('selected')).any())
        allowed={'asked_unknown','not_described','not_recorded','measurement_failed','blank_or_dependency_not_triggered'}
        for c in w:
            if c=='participant_id':continue
            wm=w[c].str.lower().isin(['true','1']);em=e[c].str.lower().isin(['true','1'])
            ck(setting+'_'+c+'_eligible_empty_only',(~wm|(em&inp[c].eq('')&~excluded)).all())
            ck(setting+'_'+c+'_no_structural_fill',states.loc[wm,c].isin(allowed).all())
            if c.startswith('index_lab'):
                ck(setting+'_'+c+'_no_blank_assays',~(wm&states[c].eq('blank_or_dependency_not_triggered')).any())
            if c.startswith('rfu_'):
                ck(setting+'_'+c+'_no_undated_historical_fill',not wm.any())
        # Recheck actual reference-standard imaging performance parent gates.
        rr=raw.loc[y.participant_id].reset_index(drop=True)
        for date,finding in p.v2.old.ECHO_CHOLECYSTITIS_EVIDENCE+p.v2.old.CT_CHOLECYSTITIS_EVIDENCE:
            parent=mappings[finding]['dependency_parent_source']
            ck(setting+'_'+finding+'_no_positive_finding_without_performance',not (rr[finding].eq('1')&~rr[parent].eq('1')).any())
        ck(setting+'_positive_has_supported_matching_endpoint',y.loc[y.gallstone_related_disease.eq('1'),'candidate_endpoint_after_exclusions'].eq('1').all())
        # Analytic branch masking must remain separate from preserved source codes.
        for _,a in d[d.type.isin(['radio','checkbox','numeric','slider'])].iterrows():
            c=a.variable
            ck(setting+'_'+c+'_hidden_or_unresolved_branches_empty',clean.loc[states[c].isin(['not_applicable_dependency','dependency_unresolved']),c].eq('').all())
        files=[OUT/f'imputed_{setting}_{i:02}.csv.gz' for i in range(1,51)]
        if not all(f.exists() for f in files) or not (OUT/f'imputed_{setting}_long.csv.gz').exists():
            pending.append(setting+'_50_complete_MI_datasets_not_yet_available');continue
        totalfilled=[];fieldcounts={};filledmasks=[]
        fullpath=OUT/f'cohort_{setting}_all_imputations.csv.gz'
        fullchunks=pd.read_csv(fullpath,dtype=str,keep_default_na=False,na_filter=False,chunksize=len(clean)) if fullpath.exists() else None
        if fullchunks is None:pending.append(setting+'_full_cohort_MI_merge_not_yet_available')
        for i,f in enumerate(files,1):
            cc=read(f);ck(f'{setting}_MI{i}_record_alignment',cc.participant_id.equals(inp.participant_id));ck(f'{setting}_MI{i}_imp_id',cc['.imp'].eq(str(i)).all());ck(f'{setting}_MI{i}_columns',set(cc)==set(inp)|{'.imp'})
            changepath=OUT/f'imputed_changed_{setting}_{i:02}.csv.gz';ch=read(changepath) if changepath.exists() else None
            if changepath.exists():
                with gzip.open(changepath,'rt') as stream:header=next(csv.reader(stream))
                ck(f'{setting}_MI{i}_change_mask_unique_headers',len(header)==len(set(header)) and header.count('participant_id')==1)
            count=0;filledmask={}
            for c in inp:
                same=equal(inp[c],cc[c]) if c!='participant_id' else inp[c].eq(cc[c])
                mask=w[c].str.lower().isin(['true','1']) if c in w and c!='participant_id' else pd.Series(False,index=inp.index)
                ck(f'{setting}_MI{i}_{c}_protected_cells',same.loc[~mask].all())
                changed=~same;fill=inp[c].eq('')&cc[c].ne('');count+=int(fill.sum());filledmask[c]=fill.to_numpy()
                if ch is not None and c in ch and c!='participant_id':ck(f'{setting}_MI{i}_{c}_change_mask',changed.equals(ch[c].str.lower().isin(['true','1'])))
                if c in spec.index and spec.loc[c,'method']:
                    lev=spec.loc[c,'levels_json']
                    if lev:ck(f'{setting}_MI{i}_{c}_valid_categories',cc.loc[fill,c].isin(json.loads(lev)).all())
                    elif spec.loc[c,'method']=='pmm':
                        donors=pd.to_numeric(inp[c],errors='coerce').dropna().to_numpy();v=pd.to_numeric(cc.loc[fill,c],errors='coerce').to_numpy()
                        ck(f'{setting}_MI{i}_{c}_observed_donor_support',all(np.isclose(z,donors,rtol=1e-12,atol=1e-12).any() for z in v))
                if c in spec.index and spec.loc[c,'method'] and mask.any():
                    method=spec.loc[c,'method']
                    eligible_rows=e[c].str.lower().isin(['true','1'])
                    samples=[('imputed',i,cc.loc[fill,c])]
                    if i==1:samples.insert(0,('observed',0,inp.loc[eligible_rows&inp[c].ne(''),c]))
                    for kind,draw,vals in samples:
                        if method=='pmm':
                            numbers=pd.to_numeric(vals,errors='coerce').dropna()
                            record={'setting':setting,'variable':c,'method':method,'sample':kind,'imp':draw,'n':len(numbers),'mean':numbers.mean(),'sd':numbers.std(),'category':'','count':'','proportion':''}
                            for q,label in [(0,'min'),(.05,'p05'),(.25,'p25'),(.5,'median'),(.75,'p75'),(.95,'p95'),(1,'max')]:record[label]=numbers.quantile(q) if len(numbers) else np.nan
                            distribution.append(record)
                        else:
                            levels=json.loads(spec.loc[c,'levels_json'])
                            for level in levels:
                                n=int(vals.eq(level).sum())
                                distribution.append({'setting':setting,'variable':c,'method':method,'sample':kind,'imp':draw,'n':len(vals),'category':level,'count':n,'proportion':n/len(vals) if len(vals) else np.nan})
                if i==1:fieldcounts[c]=int(fill.sum())
            if i==1:filledmasks=filledmask
            else:ck(f'{setting}_MI{i}_same_filled_cells_across_draws',all(np.array_equal(v,filledmasks[c]) for c,v in filledmask.items()))
            totalfilled.append(count)
            if fullchunks is not None:
                full=next(fullchunks).reset_index(drop=True)
                ck(f'{setting}_MI{i}_full_patient_alignment',full.participant_id.equals(clean.participant_id))
                ck(f'{setting}_MI{i}_full_draw_number',full['.imp'].eq(str(i)).all())
                retained=[c for c in clean if c not in w and c!='baseline_bmi_calculated']
                ck(f'{setting}_MI{i}_full_source_cells_preserved',full[retained].equals(clean[retained]))
                for c in w:
                    if c=='participant_id':continue
                    ck(f'{setting}_MI{i}_{c}_full_canonical_agreement',equal(full[c],cc[c]).all())
                ck(f'{setting}_MI{i}_full_outcome_agreement',full.gallstone_related_disease.equals(y.gallstone_related_disease))
                bmi=pd.to_numeric(full.startgewicht_analysis,errors='coerce')/(pd.to_numeric(full.lengte_analysis,errors='coerce')/100)**2
                got=pd.to_numeric(full.baseline_bmi_calculated,errors='coerce')
                ck(f'{setting}_MI{i}_passive_bmi',((bmi.isna()&got.isna())|np.isclose(bmi,got,rtol=1e-12,atol=1e-12)).all())
                ck(f'{setting}_MI{i}_passive_bmi_status',full.loc[got.isna(),'baseline_bmi_passive_status'].eq('components_missing_not_derived').all())
        long=read(OUT/f'imputed_{setting}_long.csv.gz');ck(setting+'_long_draw_counts',len(long)==50*len(inp) and long.groupby('.imp').size().eq(len(inp)).all());ck(setting+'_long_unique_ids_each_draw',not long.duplicated(['participant_id','.imp']).any())
        summaries[setting]={'filled_cells_each_draw':totalfilled,'fields_with_released_draws':sum(n>0 for n in fieldcounts.values()),'eligible_cells':int(w.drop(columns='participant_id').apply(lambda s:s.str.lower().isin(['true','1'])).sum().sum())}
    if distribution:pd.DataFrame(distribution).to_csv(OUT/'imputation_distribution_diagnostics.csv',index=False,na_rep='')
    payload={'review_role':'Separate Statistical Reviewer','implementation_sha256':sha(ROOT/'analysis/14_prepare_first_visit_v3.py'),'passed':all(x['passed'] for x in checks),'complete':not pending,'pending':pending,'checks_n':len(checks),'failed_checks':[x for x in checks if not x['passed']],'MI_summary':summaries,'checks':checks,'limitations':['Validation establishes protected-cell/applicability and release invariants, not missing-at-random assumptions or clinical endpoint sensitivity.','Unreliable targets can remain missing; MI is not a fully filled source dataset.','No prediction model fitted or reviewed.','ASA, mixed-scale urine specific gravity and unresolved-unit eGFR remain analytically quarantined pending source review; original source cells are preserved.','US/CT evidence helper omits explicit performed-parent gate; current data contain no contradictory positive findings.','Distribution summaries compare observed eligible values with imputed values; conditional differences alone do not establish an imputation defect.']}
    (OUT/'statistical_validation.json').write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps({k:payload[k] for k in ['passed','complete','pending','checks_n','failed_checks','MI_summary']},indent=2));return payload['passed']
if __name__=='__main__':raise SystemExit(0 if main() else 1)
