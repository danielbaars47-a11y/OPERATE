"""Observed-data cohort description; no adjudication or outcome imputation."""
from pathlib import Path
import pandas as pd
import numpy as np
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/first_visit_2026-10-27_v3'
O=ROOT/'results/analysis_2026-10-27_no_adjudication';O.mkdir(parents=True,exist_ok=True)
def read(p):return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
availability=[];descriptive=[];flow=[];hashes={}
for setting in ['SEH','POLI']:
 clean=read(P/f'cohort_{setting}_cleaned.csv.gz')
 y=read(P/f'outcomes_{setting}.csv')
 assert clean.participant_id.equals(y.participant_id) and clean.participant_id.is_unique
 eligible=y.gallstone_related_disease.isin(['0','1'])&y.analysis_eligible.eq('1')
 assert not (eligible & y.prior_cholecystectomy_excluded.eq('1')).any()
 assert not (eligible & y.identity_review_required.eq('1')).any()
 d=read(P/f'variable_dictionary_analysis_{setting}.csv')
 st=read(P/f'cohort_{setting}_analysis_missing_status.csv.gz')
 assert clean.participant_id.equals(st.participant_id)
 flow.extend([{'setting':setting,'stage':'first_presentations','n':len(clean)},
              {'setting':setting,'stage':'eligible_known_outcome','n':int(eligible.sum())},
              {'setting':setting,'stage':'gallstone_events','n':int((eligible&y.gallstone_related_disease.eq('1')).sum())},
              {'setting':setting,'stage':'nonevents','n':int((eligible&y.gallstone_related_disease.eq('0')).sum())}])
 for flag in ['prior_cholecystectomy_excluded','prior_cholecystectomy_history_unknown','identity_review_required']:
  flow.append({'setting':setting,'stage':flag,'n':int(y[flag].eq('1').sum())})
 for c in clean:
  if c in ['participant_id','analysis_setting']:continue
  observed=clean[c].ne('')
  for pop,mask in [('all_first_presentations',pd.Series(True,index=clean.index)),('eligible_known_outcome',eligible)]:
   missing=(mask&~observed)
   availability.append(dict(setting=setting,population=pop,variable=c,n=int(mask.sum()),observed=int((mask&observed).sum()),observed_percent=100*(mask&observed).sum()/mask.sum(),missing_reason_counts=json.dumps(st.loc[missing,c].value_counts().to_dict())))
 fields=[]
 for r in d.itertuples():
  if r.type in ['label','textarea','string']:continue
  fields.append((r.variable,r.question,r.type,r.unit,r.option_codes_json,r.source_names))
 for c,label,unit in [('age_years_at_index','Leeftijd op indexdatum','jaren'),('years_since_bariatric_surgery','Tijd sinds bariatrische operatie','jaren'),('baseline_bmi_calculated','BMI vóór bariatrische chirurgie, berekend uit startgewicht en lengte','kg/m2')]:
  fields.append((c,label,'numeric',unit,'{}','[]'))
 for v,label,typ,unit,codes,source_names in fields:
  if v not in clean:continue
  z=clean[v];obs=z.ne('');code_map=json.loads(codes) if codes else {}
  for group,mask in [('Alle',eligible),('Galsteenuitkomst',eligible&y.gallstone_related_disease.eq('1')),('Geen galsteenuitkomst',eligible&y.gallstone_related_disease.eq('0'))]:
   row=dict(setting=setting,group=group,variable=v,label=label,unit=unit,type=typ,total_n=int(mask.sum()),observed_n=int((mask&obs).sum()),missing_n=int((mask&~obs).sum()))
   if typ in ['numeric','slider']:
    values=pd.to_numeric(z[mask&obs],errors='raise')
    assert np.isfinite(values).all()
    q=values.quantile([.25,.5,.75]) if len(values) else pd.Series([np.nan]*3,index=[.25,.5,.75])
    descriptive.append(dict(row,category='',category_label='',n='',pct_among_observed='',median=q.loc[.5],q1=q.loc[.25],q3=q.loc[.75],minimum=values.min(),maximum=values.max()))
   else:
    counts=z[mask&obs].value_counts()
    if not len(counts):descriptive.append(dict(row,category='',category_label='Geen beschikbare waarnemingen',n=0,pct_among_observed='',median='',q1='',q3='',minimum='',maximum=''))
    for code,n in counts.items():
     if typ=='checkbox':lab={'0':'Niet geselecteerd binnen beantwoorde groep','1':'Geselecteerd'}.get(code,code)
     else:lab=code_map.get(code,code)
     descriptive.append(dict(row,category=code,category_label=lab,n=int(n),pct_among_observed=100*n/(mask&obs).sum(),median='',q1='',q3='',minimum='',maximum=''))
 for file in [f'cohort_{setting}_cleaned.csv.gz',f'outcomes_{setting}.csv',f'cohort_{setting}_analysis_missing_status.csv.gz']:
  hashes[file]=hashlib.sha256((P/file).read_bytes()).hexdigest()
pd.DataFrame(flow).to_csv(O/'cohort_flow_observed.csv',index=False)
pd.DataFrame(availability).to_csv(O/'availability_observed_all_fields.csv',index=False)
pd.DataFrame(descriptive).to_csv(O/'descriptive_all_verified_fields.csv',index=False)
core=['age_years_at_index','geslacht_analysis','years_since_bariatric_surgery','b_operatie_type_analysis','baseline_bmi_calculated','index_lab_crp_analysis','index_lab_bilirubin_total_analysis','index_lab_af_analysis','index_lab_alat_analysis','seh_a_waaro_in_de_rechter_bovenbuik_analysis','poli_ac_waar_in_de_rechter_bovenbuik_analysis','seh_a_misselijk_analysis','poli_a_misselijk_analysis','index_echo_cholecystitis_analysis']
table=pd.DataFrame(descriptive);table[table.variable.isin(core)].to_csv(O/'table1_observed_core.csv',index=False)
status={'source_version':'2026-10-27_v3','analysis_authorization':'User requested analyses without further adjudication on 2026-10-27','source_hashes':hashes,'outcomes_imputed':False,'observed_only_descriptive':True,'unassigned_first_presentations':len(read(P/'first_episodes_unassigned_setting.csv.gz')),'notes':['All availability and Table 1 summaries use original observed analytical values, not imputed data.','Counts for exclusion flags can overlap; they must not be summed.','No univariate association testing or outcome-driven predictor screening was performed.','Unresolved values and source fields remain retained in v3 but are not silently converted to analysis values.']}
(O/'descriptive_provenance.json').write_text(json.dumps(status,indent=2));print(pd.DataFrame(flow).to_string(index=False))
