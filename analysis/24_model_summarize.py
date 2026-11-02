#!/usr/bin/env python3
"""Summarize saved model outputs and paired same-population bootstrap contrasts.
No new models, imputations, or associations are selected here.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]/'results/analysis_2026-10-27_no_adjudication'
models=json.loads((ROOT/'model_specification.json').read_text())['models']
performance=[];coefficients=[]
for m in models:
    ident=m['id'];summary=json.loads((ROOT/f'{ident}_summary.json').read_text())
    p=pd.read_csv(ROOT/f'{ident}_performance.csv');p.insert(0,'model',ident)
    for key in ['n','events','predictor_df','bootstrap_success','bootstrap_failed']:p[key]=summary[key]
    performance.append(p)
    c=pd.read_csv(ROOT/f'{ident}_coefficients.csv');c.insert(0,'model',ident)
    c['coefficient_type']='MI-averaged penalized coefficient; bootstrap percentile interval'
    meta={p['name']:p for p in m['predictors']}
    c['predictor_transform']=c.term.map(lambda term:meta.get(term,{}).get('transform','intercept'))
    c['source_scale']=c.term.map(lambda term:meta.get(term,{}).get('scale_description','binary code' if term!='(Intercept)' else 'baseline log odds'))
    c['exp_coefficient_interpretation']=c.term.map(lambda term:'baseline odds, not an odds ratio' if term=='(Intercept)' else 'penalized odds ratio per one transformed predictor unit')
    c['term_interpretation']=c.term.replace({'RYGB':'Documented RYGB code group versus other documented surgery codes'})
    coefficients.append(c)
pd.concat(performance).to_csv(ROOT/'model_performance_summary.csv',index=False)
pd.concat(coefficients).to_csv(ROOT/'model_coefficient_summary.csv',index=False)
comparisons=[]
for base,extended in [('SEH_clinical_on_lab','SEH_lab'),('SEH_imaging_baseline','SEH_imaging'),('SEH_imaging_broad_link_sensitivity_baseline','SEH_imaging_broad_link_sensitivity')]:
    binput=pd.read_csv(ROOT/f'model_input_{base}.csv',dtype={'participant_id':str})
    einput=pd.read_csv(ROOT/f'model_input_{extended}.csv',dtype={'participant_id':str})
    assert binput.participant_id.equals(einput.participant_id) and binput.outcome.equals(einput.outcome)
    b=pd.read_csv(ROOT/f'{base}_bootstrap_metrics.csv');e=pd.read_csv(ROOT/f'{extended}_bootstrap_metrics.csv')
    paired=b.merge(e,on='replicate',suffixes=('_base','_extended'),validate='one_to_one')
    bp=pd.read_csv(ROOT/f'{base}_performance.csv').set_index('metric');ep=pd.read_csv(ROOT/f'{extended}_performance.csv').set_index('metric')
    for metric in bp.index:
        optimism=(paired[f'train_{metric}_extended']-paired[f'test_{metric}_extended'])-(paired[f'train_{metric}_base']-paired[f'test_{metric}_base'])
        testdelta=paired[f'test_{metric}_extended']-paired[f'test_{metric}_base']
        apparent=ep.loc[metric,'apparent']-bp.loc[metric,'apparent']
        comparisons.append(dict(base=base,extended=extended,metric=metric,n=len(binput),events=int(binput.outcome.sum()),paired_valid_bootstraps=len(paired),apparent_difference=apparent,mean_paired_optimism_difference=optimism.mean(),paired_optimism_corrected_difference=apparent-optimism.mean(),raw_original_evaluation_difference_p025=testdelta.quantile(.025),raw_original_evaluation_difference_p975=testdelta.quantile(.975),paired_optimism_difference_p025=optimism.quantile(.025),paired_optimism_difference_p975=optimism.quantile(.975),interval_note='Raw bootstrap distributions; not external-validation confidence limits'))
pd.DataFrame(comparisons).to_csv(ROOT/'same_population_model_comparisons.csv',index=False)
print('Summarized eight models and three prespecified same-population comparisons; no refitting.')
