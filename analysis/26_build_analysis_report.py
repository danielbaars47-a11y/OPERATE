"""Readable local report from completed, independently reviewed analysis outputs."""
from pathlib import Path
import json,html,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/analysis_2026-10-27_no_adjudication'
results=json.loads((P/'model_results.json').read_text())
spec=json.loads((P/'model_specification.json').read_text())
review=json.loads((P/'post_execution_statistical_review.json').read_text())
assert review['passed'] and review['complete']
assert len(results)==len(spec['models'])==8
labels={'SEH_clinical':'SEH: klinische kenmerken','SEH_lab':'SEH: kliniek + laboratorium (primair)','SEH_imaging':'SEH: uitbreiding met inflammatoire beeldvorming','POLI_clinical':'Poli: klinische kenmerken','SEH_imaging_broad_link_sensitivity':'SEH: ruimere koppeling beeldvorming (sensitiviteit)','SEH_imaging_baseline':'SEH: labmodel binnen strikte beeldvormingsgroep','SEH_imaging_broad_link_sensitivity_baseline':'SEH: labmodel binnen ruime beeldvormingsgroep','SEH_clinical_on_lab':'SEH: klinisch model binnen laboratoriumgroep'}
rows=[]
for r in results:
 perf=pd.read_csv(P/f"{r['id']}_performance.csv").set_index('metric')
 assert r['bootstrap_success']>=180 and r['outcome_in_imputation'] is False
 rows.append(dict(model=r['id'],omschrijving=labels[r['id']],n=r['n'],events=r['events'],predictor_parameters=r['predictor_df'],AUC_apparent=perf.loc['AUC','apparent'],AUC_corrected=perf.loc['AUC','optimism_corrected'],Brier_corrected=perf.loc['Brier','optimism_corrected'],calibration_intercept_corrected=perf.loc['calibration_intercept','optimism_corrected'],calibration_slope_corrected=perf.loc['calibration_slope','optimism_corrected'],bootstrap_success=r['bootstrap_success'],bootstrap_requested=r['bootstrap_requested']))
summary=pd.DataFrame(rows);summary.to_csv(P/'performance_overview.csv',index=False)
paired=[]
for ext,base in [('SEH_lab','SEH_clinical_on_lab'),('SEH_imaging','SEH_imaging_baseline'),('SEH_imaging_broad_link_sensitivity','SEH_imaging_broad_link_sensitivity_baseline')]:
 a=pd.read_csv(P/f'{ext}_predictions.csv',dtype={'participant_id':str});b=pd.read_csv(P/f'{base}_predictions.csv',dtype={'participant_id':str})
 assert a.participant_id.equals(b.participant_id) and a.outcome.equals(b.outcome)
 ee=summary.set_index('model').loc[ext];bb=summary.set_index('model').loc[base]
 be=pd.read_csv(P/f'{ext}_bootstrap_metrics.csv').set_index('replicate');bbm=pd.read_csv(P/f'{base}_bootstrap_metrics.csv').set_index('replicate')
 common=be.index.intersection(bbm.index);delta={}
 for met in ['AUC','Brier']:
  ep=pd.read_csv(P/f'{ext}_performance.csv').set_index('metric').loc[met,'apparent'];bp=pd.read_csv(P/f'{base}_performance.csv').set_index('metric').loc[met,'apparent']
  do=(be.loc[common,'train_'+met]-be.loc[common,'test_'+met])-(bbm.loc[common,'train_'+met]-bbm.loc[common,'test_'+met])
  delta[met]=ep-bp-do.mean()
 paired.append(dict(extension=ext,baseline=base,n=len(a),events=int(a.outcome.sum()),valid_pairs=len(common),delta_corrected_AUC=delta['AUC'],delta_corrected_Brier=delta['Brier']))
def table(headers,rows):return '<table><thead><tr>'+''.join('<th>'+html.escape(str(v))+'</th>' for v in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in r)+'</tr>' for r in rows)+'</tbody></table>'
metric_rows=[[labels[r.model],r.n,r.events,r.predictor_parameters,f'{r.AUC_apparent:.3f}',f'{r.AUC_corrected:.3f}',f'{r.Brier_corrected:.3f}',f'{r.calibration_slope_corrected:.2f}',f'{r.bootstrap_success}/{r.bootstrap_requested}'] for r in summary.itertuples()]
comparisons=[[labels[r['extension']],r['n'],r['valid_pairs'],f"{r['delta_corrected_AUC']:+.3f}",f"{r['delta_corrected_Brier']:+.3f}"] for r in paired]
desc=pd.read_csv(P/'table1_observed_core.csv',keep_default_na=False)
core=[]
for s in ['SEH','POLI']:
 for v in ['age_years_at_index','geslacht_analysis','years_since_bariatric_surgery','baseline_bmi_calculated','index_lab_bilirubin_total_analysis','index_lab_af_analysis','index_lab_alat_analysis','index_lab_crp_analysis']:
  z=desc[(desc.setting==s)&(desc.group=='Alle')&(desc.variable==v)]
  if z.empty:continue
  if v=='geslacht_analysis':z=z[z.category.astype(str).eq('0')]
  r=z.iloc[0]
  if r.type=='radio':value=f"{int(r['n'])}/{int(r.observed_n)} vrouw ({float(r.pct_among_observed):.1f}%)"
  elif r['median']=='':value='Geen waarnemingen'
  else:value=f"{float(r['median']):.1f} [{float(r.q1):.1f}–{float(r.q3):.1f}]"
  core.append([s,r.label,r.unit,value,f'{int(r.observed_n)}/{int(r.total_n)}'])
figs=''
for mid in ['SEH_lab','SEH_clinical','POLI_clinical','SEH_imaging']:
 figs+=f'<section><h3>{html.escape(labels[mid])}</h3><div class="figures"><figure><img src="{mid}_discrimination.png" alt="OOB ROC"/><figcaption>ROC op gemiddeld voorspelde kansen uit bootstrapmodellen waarin de patiënt niet in de training zat.</figcaption></figure><figure><img src="{mid}_calibration.png" alt="OOB calibratie"/><figcaption>OOB-calibratie in acht groepen; Wilson-intervallen rond geobserveerde fracties.</figcaption></figure></div></section>'
coefficient_sections=''
terms={'(Intercept)':'Intercept','age':'Leeftijd, per jaar','male':'Man versus vrouw','years_since_surgery':'Jaren sinds bariatrische operatie, per jaar','RYGB':'RYGB versus overige operatietypen','RUQ_pain':'Geregistreerde pijn rechtsboven','bilirubin':'Bilirubine, per eenheid log(1 + µmol/L)','ALP':'Alkalische fosfatase, per eenheid log(1 + U/L)','ALT':'ALAT, per eenheid log(1 + U/L)','CRP':'CRP, per eenheid log(1 + mg/L)','positive_index_CT_or_US_cholecystitis':'Positieve CT/echo voor cholecystitis'}
for mid in ['SEH_lab','POLI_clinical']:
 c=pd.read_csv(P/f'{mid}_coefficients.csv');cr=[]
 for r in c.itertuples():
  if r.term=='(Intercept)':continue
  cr.append([terms.get(r.term,r.term),f'{r.exp_coefficient:.2f}',f'{r.exp_lower:.2f}–{r.exp_upper:.2f}'])
 coefficient_sections+=f'<h3>{labels[mid]}</h3>'+table(['Predictor en schaal','Geschatte OR met ridge-penalisatie','Bootstrap-percentielinterval (95%)'],cr)
body=f'''<!doctype html><html lang="nl"><meta charset="utf-8"><title>OPERATE — analyses zonder aanvullende broncontrole</title><style>
body{{font:16px/1.55 Arial,sans-serif;color:#172b38;max-width:1200px;margin:40px auto;padding:0 24px;background:white}}h1{{font-size:30px}}h2{{font-size:23px;margin-top:36px}}h3{{font-size:18px}}table{{border-collapse:collapse;width:100%;font-size:14px;margin:16px 0}}th{{background:#eaf0f5;text-align:left;padding:10px;border-bottom:2px solid #34546c}}td{{padding:9px 10px;border-bottom:1px solid #dce3e8;vertical-align:top}}.notice{{border-left:4px solid #9c6d20;padding:12px 18px;background:#fff8e8}}.figures{{display:flex;gap:18px;flex-wrap:wrap}}figure{{margin:10px 0;max-width:440px}}img{{width:100%;height:auto}}figcaption{{font-size:13px;color:#475967}}a{{color:#155b87}}@media print{{body{{margin:0;font-size:11px}}section,table{{break-inside:avoid}}h2{{break-after:avoid}}}}</style>
<h1>OPERATE: analyses zonder aanvullende broncontrole</h1><p>27 oktober 2026. De huidige werkversie is geanalyseerd op verzoek van de onderzoeker. Geen aanvullende bronwaarden, diagnoses of identiteitsbeslissingen zijn ingevuld.</p>
<div class="notice"><b>Onderzoeksresultaten, geen gevalideerd klinisch beslisinstrument.</b> De modellen gelden voor de beschreven beschikbaarheidssubgroepen. Onopgeloste bronproblemen blijven bestaan. Interne validatie vervangt geen externe validatie.</div>
<h2>Cohort en uitkomst</h2><p>Van 1.593 eerste studie-episodes waren 789 SEH, 776 poli en 28 zonder vastgestelde setting. Na de vastgelegde uitsluitingen en het vereiste van een beoordeelbare uitkomst blijven 567 SEH-patiënten (88 events, 479 controles) en 563 polipatiënten (88 events, 475 controles). Personen met eerdere cholecystectomie worden als case én controle uitgesloten. Ook onopgeloste identiteit, onbekende relevante voorgeschiedenis en onopgeloste uitkomsten worden niet alsnog ingevuld.</p>
<p>Een event is de geregistreerde definitieve cholecystolithiasis/cholecystitis volgens de afgesproken bevestigingscriteria binnen de indexepisode. Onbekende uitkomsten worden niet als controle gecodeerd. Alle aantallen en beschikbaarheidspercentages zijn te herleiden in de bijgevoegde tabellen.</p>
<h2>Geobserveerde kenmerken</h2><p>Continue variabelen: mediaan [interkwartielafstand]. Alleen werkelijk beschikbare waarden, vóór imputatie. De noemer toont het aantal waarnemingen binnen de volledige uitkomstbekende groep; dit is niet altijd de modelspecifieke groep.</p>{table(['Setting','Kenmerk','Eenheid','Waarde','Beschikbaar'],core)}
<h2>Modellen en interne validatie</h2><p>De vooraf vastgelegde klinische kern bestaat uit leeftijd, geslacht, tijd sinds bariatrische chirurgie, RYGB versus overige operaties en geregistreerde pijn rechtsboven. Het primaire SEH-model voegt bilirubine, alkalische fosfatase, ALAT en CRP toe, elk als log(1 + waarde). Er is geen selectie op univariabele p-waarden of automatische backwardselectie uitgevoerd.</p>
<p>Alle modellen gebruiken ridge-logistische regressie met alpha = 0 en lambda = 1/n van de trainingssteekproef, met standaardisatie uitsluitend in de training. Per eindmodel zijn 50 imputaties gemaakt. Alleen 1–3 werkelijk ontbrekende waarden voor tijd sinds chirurgie zijn in deze geselecteerde groepen geïmputeerd. Overige geselecteerde predictoren moeten waargenomen zijn; structurele missings zijn nooit ingevuld.</p>
<p>Per model zijn 200 bootstrapsteekproeven gevraagd. Binnen elke steekproef worden imputatie (5 trekkingen) en regressie opnieuw uitgevoerd. Evaluatiepatiënten beïnvloeden het imputatiemodel niet en hun uitkomst wordt nooit aan de imputer aangeboden. De uitkomst is ook uit de trainingsimputatie weggelaten om de voorspelsituatie na te bootsen. Het eindmodel gebruikt dus modelspecifieke imputatie, niet de eerdere algemene blokimputaties.</p>
{table(['Model','N','Events','Parameters','AUC schijnbaar','AUC gecorrigeerd','Brier gecorrigeerd','Calibratiehelling gecorrigeerd','Geldige bootstraps'],metric_rows)}
<p>AUC is discriminatie; hoger is beter. Brier is de gemiddelde kwadratische voorspellingsfout; lager is beter. De ideale calibratiehelling is 1 en het ideale calibratie-intercept is 0. Gecorrigeerde waarden zijn aangepast voor gemiddeld bootstrapoptimisme. De afzonderlijke CSV-bestanden bevatten ook intercepten en optimismeverdelingen. De percentielgrenzen van het optimisme zijn geen betrouwbaarheidsinterval van de uiteindelijke modelprestatie.</p>
<h2>Vergelijkingen binnen dezelfde patiëntengroep</h2><p>Een verschil tussen AUC's uit verschillende populaties bewijst geen toegevoegde waarde. Hieronder staan uitsluitend vergelijkingen met opnieuw geschatte basismodellen op dezelfde patiënten.</p>{table(['Uitbreiding versus passend basismodel','N','Geldige paren','Verschil gecorrigeerde AUC','Verschil gecorrigeerde Brier'],comparisons)}
<h2>Regressiecoëfficiënten</h2><p>De onderstaande oddsratio's zijn gemiddelden van de gepenaliseerde coëfficiënten op de log-schaal, daarna geëxponentieerd. Intervallen zijn percentielintervallen uit de volledige bootstrapprocedure, geen Rubin-standaardfouten of causale effecten. Er worden geen p-waarden gebruikt om predictoren te selecteren.</p>{coefficient_sections}
<h2>Discriminatie en calibratiefiguren</h2><p>De figuren gebruiken per patiënt het gemiddelde van voorspellingen uit bootstrapmodellen waarin deze patiënt niet was opgenomen (out-of-bag, OOB). Dit zijn interne, gemiddeld voorspelde diagnostieken; het zijn niet de ROC- of calibratiecurven van één extern gevalideerd eindmodel.</p>{figs}
<h2>Beperkingen die bij deze analyse blijven gelden</h2><ul><li>Geen herstel van ASA, urine-soortelijk gewicht, onzekere eGFR-eenheden of overige onopgeloste bronwaarden. Deze zijn niet ongemerkt als modelpredictor gebruikt.</li><li>Het laboratoriummodel betreft informatie uit de eerste diagnostische beoordeling; een gelijke kalenderdatum bewijst niet dat alle uitslagen al vóór de eerste klinische beoordeling bekend waren.</li><li>De beeldvormingsuitbreiding betreft alleen betrouwbare cholecystitisbevindingen op CT/echo. De belangrijkste CT-steenvariabele heeft geen voldoende geverifieerde codering voor modelgebruik; dit is geen volledig galsteenbeeldvormingsmodel.</li><li>Beeldvorming kan deel uitmaken van de definitieve uitkomstbepaling. Daardoor dreigen incorporation- en verificatiebias. De strikte beeldvormingsgroep bevat slechts 38 events voor 10 parameters, met aanzienlijk overfittingsrisico ondanks penalisatie.</li><li>De ruimere beeldvormingsanalyse neemt expliciet gemarkeerde koppelconflicten mee als sensitiviteitsanalyse. Zij vervangt de strikte analyse niet.</li><li>De poliklinische analyse betreft uitsluitend patiënten met de vereiste gegevens uit de actuele klachtenepisode. Ontbrekende checkboxgroepen zijn niet tot 'nee' gemaakt of structureel geïmputeerd.</li><li>Beschikbaarheid en uitkomstbeoordeling kunnen selectief zijn; generaliseer deze prestaties niet zonder meer naar alle patiënten met buikpijn na bariatrische chirurgie.</li></ul>
<h2>Bestanden en methodebronnen</h2><p>De afzonderlijke nacontrole is afgerond: {review['checks_n']:,} controles geslaagd, inclusief herberekening vanuit de opgeslagen modellen en controle van imputatie, bootstraplidmaatschap en gepaarde vergelijkingen.</p><ul><li><a href="performance_overview.csv">Modelprestaties</a> en <a href="same_population_model_comparisons.csv">vergelijkingen op dezelfde populatie</a></li><li><a href="table1_observed_core.csv">Geobserveerde cohortkenmerken</a> en <a href="descriptive_all_verified_fields.csv">alle geverifieerde velden beschrijvend</a></li><li><a href="availability_observed_all_fields.csv">Beschikbaarheid van alle behouden velden</a></li><li><a href="frozen_analysis_cohort.csv">Bevroren inclusies en uitsluitingen</a></li><li><a href="model_specification.json">Vooraf vastgelegde modelspecificatie</a> en <a href="statistical_review.json">afzonderlijke statistische beoordeling</a></li><li>Per model: coëfficiënten, voorspellingen, 50 voltooide datasets, R-objecten, bootstrapstatus en figuren in PNG/PDF.</li><li><a href="https://amices.org/mice/reference/mice.html">MICE: imputatie en uitsluiten van evaluatierijen uit het imputatiemodel</a>; <a href="https://amices.org/mice/reference/pool.html">Pooling van schattingen in plaats van ingevulde datasets</a>.</li></ul></html>'''
(P/'OPERATE_analyserapport.html').write_text(body,encoding='utf-8')
manifest={f.name:{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in P.iterdir() if f.is_file() and f.name!='analysis_manifest.json'}
(P/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2))
print(summary[['model','n','events','AUC_corrected','calibration_slope_corrected','bootstrap_success']].to_string(index=False))
