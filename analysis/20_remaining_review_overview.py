"""Read-only counts for remaining adjudication; output is separate from frozen v3."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/first_visit_2026-10-27_v3'
O=ROOT/'results/review_overview_2026-10-27';O.mkdir(exist_ok=True)
def read(p):return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
q=read(P/'clinical_review_queue.csv');counts=[];details=[];metadata=[];remaining=[]
def add(issue,setting,ids,variables):
 ids=set(ids);counts.append(dict(issue=issue,setting=setting,unique_study_ids=len(ids),variables='; '.join(sorted(set(variables)))))
 details.extend(dict(issue=issue,setting=setting,participant_id=i) for i in sorted(ids))
for s in ['SEH','POLI']:
 z=q[q.setting.eq(s)]
 for issue,pattern in [('ASA_field','^b_asa$'),('urine_specific_gravity','^index_urine.*_sg$'),('egfr_unit','^index_lab.*egfr'),('urine_pH','^index_urine.*_ph$'),('blood_pressure_198_review_only','^seh_lo_drr$')]:
  zz=z[z.source_variable.str.contains(pattern,regex=True)];add(issue,s,zz.participant_id,zz.source_variable)
 asa=z[z.source_variable.eq('b_asa')];bad=asa[~pd.to_numeric(asa.source_value,errors='coerce').isin(range(1,7))]
 add('ASA_outside_classes',s,bad.participant_id,bad.source_variable)
 d=read(P/f'cohort_{s}_cleaned.csv.gz');st=read(P/f'cohort_{s}_analysis_missing_status.csv.gz');di=read(P/f'variable_dictionary_source_{s}.csv')
 f=read(P/'quality_flags.csv');f=f[f.setting.eq(s)]
 for issue,mask in [('other_vitals_source_check',f.variable.str.startswith('SEH_LO_')),('other_laboratory_source_check',f.variable.str.startswith('index_lab_')),('date_plausibility',f.variable.isin(['age_years_at_index','years_since_bariatric_surgery','B_operatie_datum']))]:
  zz=f[mask];add(issue,s,zz['Participant Id'],zz.variable)
 for c in ['index_date_iso','age_years_at_index','years_since_bariatric_surgery']:
  add('missing_or_withheld_'+c,s,d.loc[d[c].eq(''),'participant_id'],[c])
 date_ids=set(f.loc[f.variable.isin(['age_years_at_index','years_since_bariatric_surgery','B_operatie_datum']),'Participant Id'])
 for c in ['index_date_iso','age_years_at_index','years_since_bariatric_surgery']:
  date_ids.update(d.loc[d[c].eq(''),'participant_id'])
 add('dates_needing_source_review',s,date_ids,['index_date_iso','geboortedatum','b_operatie_datum'])
 valid_roles=['baseline_or_first_visit_source_field','historical_information_with_timing_gate','same_day_laboratory','same_day_diagnostic_information','same_day_questionnaire']
 cand=di[di.role.isin(valid_roles)&di.numeric_companion.eq('')&di.metadata_status.isin(['unresolved','provisional'])]
 affected=set()
 for r in cand.itertuples():
  if r.variable not in d:continue
  if r.variable+'_analysis' in d:continue
  has=d[r.variable].ne('') & st[r.variable].isin(['present_code_or_text','numeric_observed'])
  if not has.any():continue
  ids=set(d.loc[has,'participant_id']);affected|=ids
  metadata.append(dict(setting=s,variable=r.variable,source_columns=r.source_columns,role=r.role,metadata_type=r.metadata_type,metadata_status=r.metadata_status,unique_study_ids=len(ids)))
 add('unverified_metadata_with_recorded_values',s,affected,[r['variable'] for r in metadata if r['setting']==s])
 w=read(P/f'where_{s}.csv.gz');imp=read(P/f'imputed_{s}_01.csv.gz');assert d.participant_id.equals(imp.participant_id)
 ids=set();cells=0
 for c in w:
  if c=='participant_id':continue
  mask=w[c].str.lower().eq('true')&imp[c].eq('')
  if mask.any():
   ids|=set(d.loc[mask,'participant_id']);cells+=int(mask.sum())
   remaining.append(dict(setting=s,variable=c,patients=int(mask.sum())))
 add('eligible_missing_unreleased',s,ids,[x['variable'] for x in remaining if x['setting']==s])
 counts[-1]['cells']=cells
pd.DataFrame(counts).to_csv(O/'additional_review_counts.csv',index=False)
pd.DataFrame(details).to_csv(O/'additional_review_patient_ids.csv',index=False)
pd.DataFrame(metadata).to_csv(O/'unverified_recorded_variables.csv',index=False)
pd.DataFrame(remaining).to_csv(O/'remaining_missing_by_variable.csv',index=False)
print(pd.DataFrame(counts).drop(columns='variables').to_string(index=False))
print('Metadata fields by setting',pd.DataFrame(metadata).groupby('setting').size().to_dict())
# Human-facing overview; no data corrections or changes to the released MI.
a=json.loads((O/'cohort_adjudication_summary.json').read_text())
rows=[]
def point(title,seh,poli,unknown,variables,needed,note=''):
 rows.append(dict(beoordelingspunt=title,SEH_studie_IDs=seh,POLI_studie_IDs=poli,setting_onbekend_studie_IDs=unknown,variabelen=variables,benodigde_informatie=needed,toelichting=note))
point('Mogelijke dubbele personen',1,6,1,'Participant Id; Episode_hoeveel; Episode_nummer_*','Bevestig per paar dezelfde of verschillende persoon; bij dezelfde persoon welke episode de eerste is en welke records samen horen.','4 paren / 8 studie-IDs; nog geen 8 bewezen unieke personen.')
point('Onbekende eerste setting',0,0,28,'SEH_1stecontact; POLI_1stecontact; SEH_datum; POLI_datum','Bevestig eerste contactsetting en indexdatum vanuit het dossier.','27 beide vlaggen leeg; 1 beide nee; 1 ID overlapt de mogelijke dubbelen.')
point('Onbekende eerdere cholecystectomie',3,3,0,'VG_cholecystectomie; operatiedatum','Was de galblaas al verwijderd vóór indexpresentatie? Zo ja, datum of bevestiging van de volgorde.')
issue_labels={
'outcome_assessment_missing':('Ontbrekende uitkomstregistratie','Aanvullen of follow-upbeoordeling is verricht, definitieve diagnose(s), bewijs en afsluitdatum; niet automatisch als negatief coderen.'),
'outcome_chronology_or_source_conflict':('Tegenstrijdige uitkomst/datum','Bevestig juiste indexdatum, diagnose- en follow-updatum en welk diagnoseveld leidend is op basis van het dossier.'),
'other_possible_gallstone_diagnosis_review':('Andere mogelijke galsteendiagnose','Adjudiceer per geval of dit voldoet aan de afgesproken cholecystolithiasis/cholecystitis-uitkomst. Een bredere uitkomstdefinitie vergt een expliciete protocolkeuze.'),
'definitive_label_without_verifiable_in_episode_evidence':('Diagnoselabel zonder verifieerbaar bewijs','Lever datum en uitslag van kwalificerende echo/CT/operatie/ERCP binnen de indexepisode, of corrigeer het definitieve diagnoselabel.')}
for code,(title,ask) in issue_labels.items():
 vals={z['setting']:z['unique_participant_ids'] for z in a['terminal_review_reason_counts'] if z['issue']==code}
 point(title,vals.get('SEH',0),vals.get('POLI',0),0,'; '.join(a['variables_by_issue'][code]),ask)
point('Onderzoek op volgende kalenderdag',248,30,0,'index_date; test_date; source_instance; modality; start/einde bezoek','Bevestig via bezoektijden of bezoek-ID of het onderzoek bij hetzelfde indexbezoek hoort en wanneer de uitslag beschikbaar was.','428 onderzoeken bij 278 IDs. Voor primaire analyse na bekende eerdere cholecystectomie: 196 SEH / 20 POLI. Geen automatische uitbreiding met één dag.')
lookup={(x['issue'],x['setting']):x for x in counts}
def clinical(code,title,variables,ask,note=''):
 point(title,lookup[(code,'SEH')]['unique_study_ids'],lookup[(code,'POLI')]['unique_study_ids'],0,variables,ask,note)
clinical('ASA_field','ASA-classificatie','B_ASA / b_asa','Bevestig schaal/codering en herstel de vier waarden 163,168,165,173 vanuit de bron. Bevestig dat de overige 1–4-waarden echte ASA-klassen zijn.','815 aanwezige waarden veldbreed achtergehouden; 4 afwijkende codes (2 SEH/2 POLI); 811 waarden niet als fout aangemerkt.')
clinical('egfr_unit','eGFR-eenheid','index_lab_egfr; index_lab_repeat*_egfr; index_lab_egfrmdrd','Bevestig de werkelijke eenheid, berekeningsmethode en eventuele verslaggrenzen; verklaar de Castor-unit kl/1.73m2.','757 IDs met een aanwezige bronwaarde; herhaalde metingen tellen niet als extra patiënten.')
clinical('urine_specific_gravity','Urine-soortelijk gewicht','index_urine_sg; index_urine_repeat*_sg','Bevestig of 1005–1035 staat voor 1.005–1.035, of andere schaal/codering. Geef een geldige omzetregel en uitzonderingen.')
clinical('urine_pH','Urine-pH','index_urine_ph','Controleer de twee waarden 60 in de bron en geef de juiste waarde of bevestig dat deze onbekend blijft.')
clinical('other_vitals_source_check','Overige vitale waarden','SEH_LO_AF; SEH_LO_DRR; SEH_LO_Sat; SEH_LO_Temp; SEH_LO_pols','Controleer geregistreerde waarden en decimalen/eenheden in het dossier: AF 1; DRR 774; saturatie 00; temperatuur 3.2/377/62; pols 870.','7 SEH-IDs; niet automatisch decimalen herstellen.')
clinical('blood_pressure_198_review_only','Hoge diastolische bloeddruk: controle zonder uitsluiting','SEH_LO_DRR','Bevestig of 198 mmHg juist is; de systolische waarde ligt hoger. Indien correct blijft de meting behouden.','Geen vastgestelde fout; waarde is niet uitgesloten.')
clinical('other_laboratory_source_check','Overige laboratoriumwaarden','index_lab_MCV; index_lab_kalium; index_lab_calcium; index_lab_glucose','Controleer originele uitslagen en eenheden voor MCV 12.7, kalium 47, calcium 46 en glucose 83. Geef juiste waarden of bevestig bronwaarde.')
clinical('dates_needing_source_review','Ontbrekende/onlogische basisdatums','index_date_iso; geboortedatum; B_operatie_datum; age_years_at_index; years_since_bariatric_surgery','Bevestig geboorte-, bariatrische operatie- en indexdatum, inclusief welke operatie als bariatrische voorgeschiedenis geldt.','Union ontbrekende/achtergehouden afgeleide leeftijden/intervallen en datumplausibiliteitsflags; overlap met uitkomstconflicten mogelijk.')
clinical('unverified_metadata_with_recorded_values','Aanvullende definities van nog niet geverifieerde bronvelden','Zie unverified_recorded_variables.csv; onder andere M_* medicatievelden en B_compl_30dgn_*','Voor de benoemde velden: aanvulling met huidige vraagtekst, antwoordcodes en afhankelijkheden. Waar dit uit bestaande documenten verifieerbaar is, kan de koppeling zonder patiëntdossier worden opgelost.','342 SEH / 318 POLI bronvelden met opgeslagen waarden en zonder geverifieerde analytische variant. Opgeslagen 0-codes tellen mee; dit zijn geen dossiercontroles voor alle patiënten.')
clinical('eligible_missing_unreleased','Legitieme missings die niet betrouwbaar konden worden geïmputeerd','Zie remaining_missing_by_variable.csv','Alleen aanvullende werkelijk waargenomen brongegevens als die beschikbaar zijn. Anders ontbrekend laten; model-specifieke aanpak later met statisticus bepalen.','1217 cellen bij 504 SEH-IDs; 645 cellen bij 407 POLI-IDs. Dit zijn geen bewezen registratie- of imputatiefouten.')
pd.DataFrame(rows).to_csv(O/'beoordelingspunten_met_vragen.csv',index=False)
