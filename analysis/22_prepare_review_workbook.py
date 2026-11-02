"""Prepare review-only workbook data. Never updates clinical data or imputation."""
from pathlib import Path
import json,re,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/first_visit_2026-10-27_v3';R=ROOT/'results/review_overview_2026-10-27'
OUT=ROOT/'outputs/01a0f3d7-80be-7ae0-97c0-484d081c78dc';OUT.mkdir(parents=True,exist_ok=True)
def read(p):return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
raw=pd.read_csv(ROOT/'OPERATE data/Buikpijn_na_bariatrische_chirurg_export_20260622.csv',sep=';',dtype=str,keep_default_na=False,na_filter=False).set_index('Participant Id',drop=False)
mappings=json.loads((P/'codebook/mappings.json').read_text())
assert raw.index.is_unique
clean={s:read(P/f'cohort_{s}_cleaned.csv.gz').set_index('participant_id',drop=False) for s in ['SEH','POLI']}
outcomes={s:read(P/f'outcomes_{s}.csv').set_index('participant_id') for s in clean}
dicts={s:read(P/f'variable_dictionary_source_{s}.csv').set_index('variable') for s in clean}
adicts={s:read(P/f'variable_dictionary_analysis_{s}.csv').set_index('variable') for s in clean}
def exclusion(s,i):
 o=outcomes[s].loc[i]
 return 'Eerdere cholecystectomie' if str(o.prior_cholecystectomy_excluded)=='1' else ''
def orig(s,i,col,candidates=None):
 d=clean[s].loc[i]
 if candidates is None:
  candidates=json.loads(dicts[s].loc[col,'source_columns']) if col in dicts[s].index else []
 if not candidates:
  candidates=[c for c in raw if c.lower()==col.lower()]
 if len(candidates)==1 and candidates[0] in raw:
  c=candidates[0];return c,raw.at[i,c]
 match=re.fullmatch(r'(index_[a-z]+(?:_repeat\d+)?)_(.+)',col)
 if match:
  selector=d.get(match[1]+'_source_instance','')
  hits=[c for c in candidates if selector and c.startswith(selector+'_')]
  if len(hits)==1:return hits[0],raw.at[i,hits[0]]
 return col,d.get(col,'')
def row(s,i,issue,var,value,question,context='',kind='Dossiercontrole',source=''):
 return dict(study_id=i,setting=s,issue=issue,variables=var,original_value=value if value!='' else '[leeg]',question=question,context=context,priority=kind,existing_exclusion=exclusion(s,i),source_reference=source or 'Castor CSV 20260622: '+var)
rows=[]
q=read(P/'clinical_review_queue.csv')
for r in q.itertuples():
 s=r.setting;i=r.participant_id;c=r.source_variable
 src,val=orig(s,i,c,json.loads(r.original_source_names))
 if c=='b_asa':
  bad=pd.to_numeric(val,errors='coerce') not in range(1,7)
  issue='ASA: afwijkende code' if bad else 'ASA: veldafspraak G01'
  ask='Controleer de juiste ASA-klasse in de bron; leg geen vermoede correctie vast.' if bad else 'Volg besluit G01. Alleen een individuele correctie invullen als het dossier daartoe aanleiding geeft.'
  ctx='Veldbreed achtergehouden. Een waarde 1–6 is op zichzelf niet als fout aangemerkt.'
  kind='Dossiercontrole' if bad else 'Veldafspraak'
 elif 'egfr' in c:
  issue='eGFR: veldafspraak G02';ask='Volg besluit G02 over eenheid/methode. Noteer eventuele individuele uitzondering.';ctx='Geregistreerde eenheid kl/1.73m2 is onopgelost.';kind='Veldafspraak'
 elif c.endswith('_sg'):
  issue='Urine-SG: veldafspraak G03';ask='Volg besluit G03 over de gebruikte schaal. Noteer eventuele individuele uitzondering.';ctx='Zowel 1.xxx als 10xx komt voor; nog geen omzetting toegepast.';kind='Veldafspraak'
 elif c.endswith('_ph'):
  issue='Urine-pH buiten bereik';ask='Controleer de oorspronkelijke uitslag. Geef juiste waarde of onbekend; geen decimaal raden.';ctx='Analytische waarde achtergehouden.';kind='Dossiercontrole'
 else:
  issue='Diastolische druk 198: bevestiging';ask='Bevestig 198 mmHg of geef broncorrectie. De systolische druk is hoger.';ctx='Controlepunt; oorspronkelijke waarde is behouden.';kind='Dossiercontrole'
 rows.append(row(s,i,issue,src,val,ask,ctx,kind))
# Other numeric source flags; dates are consolidated below.
f=read(P/'quality_flags.csv')
for r in f.itertuples(index=False,name=None):
 rec=dict(zip(f.columns,r));s=rec['setting'];i=rec['Participant Id'];v=rec['variable']
 if not (v.startswith('SEH_LO_') or v.startswith('index_lab_')):continue
 col=v.lower();src,val=orig(s,i,col)
 rows.append(row(s,i,'Afwijkende meetwaarde',src,val,'Controleer originele uitslag, eenheid en decimale notatie. Vul de juiste waarde in of bevestig dat deze onbekend blijft.',f"Signaal uit plausibiliteitscontrole: {rec['value']}. Geen correctie geraden."))
ids=read(R/'additional_review_patient_ids.csv')
for r in ids[ids.issue.eq('dates_needing_source_review')].itertuples():
 s=r.setting;i=r.participant_id
 cols=['Geboortedatum','B_operatie_datum','SEH_datum' if s=='SEH' else 'POLI_datum']
 # Actual source names are taken from the stored source dictionary.
 pairs=[orig(s,i,c) for c in ['geboortedatum','b_operatie_datum','seh_datum' if s=='SEH' else 'poli_datum']]
 text='\n'.join(c+' = '+(v if v else '[leeg]') for c,v in pairs)
 rows.append(row(s,i,'Basisdatums controleren','; '.join(c for c,v in pairs),text,'Bevestig geboorte-, bariatrische operatie- en indexdatum. Geef bij correctie telkens veldnaam en YYYY-MM-DD.', 'Ontbrekende of onlogische leeftijd/operatie-interval; meerdere velden kunnen betrokken zijn.'))
cohort=json.loads((R/'cohort_review_rows.json').read_text())
labels={'definitive_label_without_verifiable_in_episode_evidence':'Diagnoselabel zonder verifieerbaar bewijs','first_presentation_setting_unresolved':'Eerste contactsetting onbekend','other_possible_gallstone_diagnosis_review':'Andere mogelijke galsteendiagnose','outcome_assessment_missing':'Uitkomstregistratie ontbreekt','outcome_chronology_or_source_conflict':'Tegenstrijdige diagnose of datums','possible_duplicate_person_identity_unresolved':'Mogelijke dubbele patiënt','prior_cholecystectomy_history_unresolved':'Eerdere cholecystectomie onbekend','overnight_next_day_linkage':'Onderzoek op volgende kalenderdag'}
for cr in cohort:
 cr['issue']=labels[cr['issue']]
 if cr['priority']=='Bezoekkoppeling':
  lines=cr['original_value'].split('\n')
  # This question concerns the visit linkage, so retain its exact source dates
  # and performance flag, without repeating the full test panel.
  cr['original_value']='; '.join(lines[:3])
  cr['variables']='; '.join(cr['variables'].split('; ')[:3])
 else:
  # Present the constructed field/value list compactly; no value is corrected.
  cr['original_value']=cr['original_value'].replace('\n','; ')
rows=cohort+rows
optional=[]
for s in clean:
 w=read(P/f'where_{s}.csv.gz').set_index('participant_id');im=read(P/f'imputed_{s}_01.csv.gz').set_index('participant_id')
 for c in w:
  if c not in im:continue
  mask=w[c].str.lower().eq('true')&im[c].eq('')
  for i in w.index[mask]:
   a=adicts[s].loc[c];src,val=orig(s,i,a.source_variable,json.loads(a.source_names))
   m=mappings.get(src,{})
   label=m.get('codes',{}).get(val,'') if m.get('verified_crosswalk') else ''
   context=str(a.get('question',''))+' | Geen betrouwbare imputatie vrijgegeven.'
   if label:context+=' Oorspronkelijk antwoordlabel: '+label+'.'
   optional.append(row(s,i,'Optioneel: ontbrekend brongegeven',src,val,'Alleen aanvullen als de werkelijk gemeten/gevraagde waarde in de bron beschikbaar is. Anders: onbekend blijft.',context,'Optionele bronaanvulling'))
assert len(optional)==1217+645
globalrows=[
 ['G01','ASA-classificatie','B_ASA',815,'Codes 1, 2, 3, 4 en afwijkend 163, 168, 165, 173','Bevestig definitie en codering. Vier afwijkende codes vragen daarnaast individuele broncontrole.','','','',''],
 ['G02','eGFR-eenheid en methode','Lab*_EGFR; eventuele MDRD-velden',757,'Bronunit: kl/1.73m2','Bevestig werkelijke eenheid, methode en omgang met rapportagegrenzen. Geef omzetregel indien aantoonbaar nodig.','','','',''],
 ['G03','Schaal urine-soortelijk gewicht','Urine*_SG',472,'Zowel 1.005–1.035 als 1005–1035 komt voor','Bevestig betekenis van beide schalen, geldige omzetregel en uitzonderingen.','','','','']]
meta=read(R/'unverified_recorded_variables.csv')
for n,(v,g) in enumerate(meta.groupby('variable',sort=True),1):
 sources=sorted(set(x for z in g.source_columns for x in json.loads(z)))
 numeric_codes=set()
 for src in sources:
  if src in raw:
   numeric_codes.update(x for x in raw[src].unique() if re.fullmatch(r'-?\d+(?:[.,]\d+)?',x))
 codes=', '.join(sorted(numeric_codes)) if len(numeric_codes)<=20 else 'Meer dan 20 numerieke bronwaarden; betekenis/type verifiëren.'
 nids=sum(int(z) for z in g.unique_study_ids)
 globalrows.append([f'M{n:03d}','Aanvullende velddefinitie', '; '.join(sources) or v,nids,codes or 'Vraagtype/tekstveld nog verifiëren','Voor dit veld: vraagtekst, type, antwoordcodes en afhankelijkheden bevestigen. Geen afzonderlijke dossiercontrole voor alle betrokken IDs nodig.','','','',''])
# Stable row identifiers preserve corrections through sorting and later linkage.
headers=['Regel-ID','Studie-ID','Cohort','Beoordelingspunt','Variabele(n)','Oorspronkelijke waarde (tekst)','Juiste waarde / beslissing','Toelichting / geraadpleegde bron','Beoordelaar','Beoordelingsdatum','Status','Gevraagde controle','Context','Bestaande uitsluiting','Bronverwijzing','Soort controle']
def matrix(items,prefix):
 result=[];seen=set()
 for r in sorted(items,key=lambda r:(r['setting'],r['study_id'],r['issue'],r['variables'],r['original_value'])):
  ident=prefix+hashlib.sha256(json.dumps(r,sort_keys=True,ensure_ascii=False).encode()).hexdigest()[:12]
  assert ident not in seen;seen.add(ident)
  result.append([ident,r['study_id'],r['setting'],r['issue'],r['variables'],r['original_value'],'','','','','',r['question'],r.get('context',''),r.get('existing_exclusion',''),r['source_reference'],r['priority']])
 return result
main=matrix(rows,'P');opt=matrix(optional,'A')
payload=dict(headers=headers,patienten=main,bronaanvulling=opt,veldheaders=['Afspraak-ID','Onderwerp','Variabelen','Betrokken studie-IDs','Huidige registratie','Te beantwoorden vraag','Juiste definitie / beslissing','Toelichting / bron','Beoordelaar','Beoordelingsdatum'],veldafspraken=globalrows)
assert all(x[1] in raw.index for x in main+opt)
assert all(all(v=='' for v in x[6:11]) for x in main+opt)
assert all(len(str(v))<32767 for x in main+opt+globalrows for v in x)
(OUT/'review_workbook_data.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2))
summary={'patient_rows':len(main),'patient_unique_ids':len({x[1] for x in main}),'optional_rows':len(opt),'optional_unique_ids':len({x[1] for x in opt}),'field_agreements':len(globalrows),'source_clinical_rows':len(q),'source_cohort_rows':len(cohort),'max_source_text_length':max(len(x[5]) for x in main+opt),'all_response_cells_blank':True}
(OUT/'preparation_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
