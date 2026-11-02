#!/usr/bin/env python3
"""Build frozen model-specific observed inputs; never reuse block-imputed v3.
Selection derives solely from prespecified population flags, never associations.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/analysis_2026-10-27_no_adjudication'
SOURCE=ROOT/'results/first_visit_2026-10-27_v3'
def read(p):return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
def main():
    spec=json.loads((OUT/'model_specification.json').read_text())
    flags=read(OUT/'model_population_flags.csv')
    frozen=read(OUT/'frozen_analysis_cohort.csv')
    flags=flags.merge(frozen[['participant_id','setting','finalY']],on=['participant_id','setting'],how='left',validate='one_to_one')
    assert not flags.duplicated(['setting','participant_id']).any()
    audit=[];receipt=[]
    for model in spec['models']:
        setting=model['setting'];flag=model['population_flag']
        use=flags.setting.eq(setting)&flags[flag].isin(['1','True','true'])
        f=flags.loc[use].copy()
        clean=read(SOURCE/f'cohort_{setting}_cleaned.csv.gz').set_index('participant_id')
        states=read(SOURCE/f'cohort_{setting}_analysis_missing_status.csv.gz').set_index('participant_id')
        assert clean.index.is_unique
        c=clean.loc[f.participant_id].reset_index();st=states.loc[f.participant_id].reset_index()
        ycol='finalY' if 'finalY' in f else 'outcome'
        frame=pd.DataFrame({'participant_id':f.participant_id.to_numpy(),'outcome':pd.to_numeric(f[ycol]).to_numpy()})
        for p in model['predictors']:
            source=p['source']
            if isinstance(source,list):
                parts=[]
                for col in source:
                    val=c[col].copy()
                    gate=p.get('source_link_gates',{}).get(col)
                    if gate:val=val.where(c[gate['column']].isin(gate['allowed']),'')
                    assert val.isin(['','0','1']).all()
                    parts.append(val)
                z=pd.concat(parts,axis=1)
                val=pd.Series(np.where(z.eq('1').any(axis=1),'1',np.where(z.eq('0').any(axis=1),'0','')))
            else:
                val=c[source].copy()
                if p.get('recode'):
                    assert val.isin(['',*p['recode']]).all(),(model['id'],source)
                    val=val.map(p['recode']).fillna('')
                if val.eq('').any():
                    assert p['method']=='pmm',(model['id'],p['name'],'missing nonimputable predictor')
                    assert st.loc[val.eq(''),source].isin(p['eligible_missing_status']).all(),(model['id'],source,'protected missing state')
            if p['type']=='binary':assert val.isin(['','0','1']).all()
            numeric=pd.to_numeric(val.replace('',np.nan),errors='raise')
            assert np.isfinite(numeric.dropna()).all()
            if p['transform']=='log1p':
                assert numeric.ge(0).all();numeric=np.log1p(numeric)
            frame[p['name']]=numeric.to_numpy()
            audit.append(dict(model=model['id'],predictor=p['name'],source=json.dumps(source),transform=p['transform'],observed_n=int(numeric.notna().sum()),missing_n=int(numeric.isna().sum()),minimum=numeric.min(),maximum=numeric.max()))
        assert len(frame)==model['expected_n'],(model['id'],len(frame),model['expected_n'])
        assert int(frame.outcome.sum())==model['expected_events']
        assert int(frame.drop(columns=['participant_id','outcome']).isna().sum().sum())==model['expected_imputable_time_cells']
        frame.to_csv(OUT/f"model_input_{model['id']}.csv",index=False,na_rep='')
        receipt.append(dict(model=model['id'],n=len(frame),events=int(frame.outcome.sum()),missing_cells=int(frame.isna().sum().sum()),columns=list(frame)))
    pd.DataFrame(audit).to_csv(OUT/'model_input_availability.csv',index=False)
    (OUT/'model_input_receipt.json').write_text(json.dumps({'models':receipt,'source':'observed cleaned v3; no existing MI reused','specification_sha256':hashlib.sha256((OUT/'model_specification.json').read_bytes()).hexdigest()},indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
