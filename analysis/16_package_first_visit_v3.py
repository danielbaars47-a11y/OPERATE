#!/usr/bin/env python3
"""Assemble full v3 cohorts for every MI draw; preserve source and status files.
Only canonical where-permitted empty cells are filled. BMI is derived afresh
from the completed baseline weight and height, never stochastically imputed.
"""
from pathlib import Path
import gzip
import hashlib
import json
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/first_visit_2026-10-27_v3'
SETTINGS=('SEH','POLI')
def read(p):return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()
def align(reference,other):
    assert reference.participant_id.is_unique and other.participant_id.is_unique
    assert set(reference.participant_id)==set(other.participant_id)
    return other.set_index('participant_id').loc[reference.participant_id].reset_index()
def refresh_manifest():
    manifest={str(p.relative_to(OUT)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='package_manifest.json'}
    (OUT/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
def main():
    cfg=json.loads((OUT/'configuration.json').read_text())
    assert cfg['outcome_mapping_status']=='validated'
    expected=[OUT/f'imputed_{s}_{i:02d}.csv.gz' for s in SETTINGS for i in range(1,51)]
    expected += [OUT/f'mice_{s}.rds' for s in SETTINGS]
    missing=[p.name for p in expected if not p.exists()]
    if missing:raise RuntimeError('MI outputs incomplete; packaging not started: '+', '.join(missing[:5]))
    clean_summary=json.loads((OUT/'cleaning_validation.json').read_text())
    summary={}
    for setting in SETTINGS:
        clean=read(OUT/f'cohort_{setting}_cleaned.csv.gz')
        canonical=align(clean,read(OUT/f'imputation_input_{setting}.csv.gz'))
        where=align(clean,read(OUT/f'where_{setting}.csv.gz'))
        outcomes=align(clean,read(OUT/f'outcomes_{setting}.csv'))
        inventory=read(OUT/f'imputation_inventory_{setting}.csv')
        allowed={c:where[c].str.lower().isin(['true','1']) for c in where if c!='participant_id'}
        eligible_cells=sum(int(z.sum()) for z in allowed.values())
        filled_draws=[]; passive_draws=[]; variables_filled=set(); first_fill=None
        target=OUT/f'cohort_{setting}_all_imputations.csv.gz';temporary=target.with_suffix('.gz.tmp')
        with gzip.open(temporary,'wt',encoding='utf-8',newline='') as handle:
            for i in range(1,51):
                values=align(clean,read(OUT/f'imputed_{setting}_{i:02d}.csv.gz'))
                assert values['.imp'].eq(str(i)).all()
                assert values.gallstone_related_disease.equals(outcomes.gallstone_related_disease)
                frame=clean.copy(); filled=0
                for col,permit in allowed.items():
                    if col not in values:continue
                    original=canonical[col]
                    use=original.eq('') & values[col].ne('')
                    assert (~use|permit).all(),col
                    assert col in clean or not use.any(),col
                    if use.any():
                        assert clean.loc[use,col].eq('').all(),col
                        frame.loc[use,col]=values.loc[use,col]
                        variables_filled.add(col);filled+=int(use.sum())
                if first_fill is None:first_fill=filled
                assert filled==first_fill,'Different released cell counts between draws'
                weight=pd.to_numeric(frame.startgewicht_analysis,errors='coerce')
                height=pd.to_numeric(frame.lengte_analysis,errors='coerce')/100
                bmi=(weight/height.pow(2)).where((weight>0)&(height>0))
                bmi=bmi.where(np.isfinite(bmi))
                previous=pd.to_numeric(clean.baseline_bmi_calculated,errors='coerce')
                changed=(bmi.notna() & (previous.isna()|~np.isclose(bmi,previous,rtol=1e-12,atol=1e-12))) | (bmi.isna()&previous.notna())
                imputed_components=(clean.startgewicht_analysis.eq('')&frame.startgewicht_analysis.ne(''))|(clean.lengte_analysis.eq('')&frame.lengte_analysis.ne(''))
                frame['baseline_bmi_calculated']=bmi.map(lambda z:'' if pd.isna(z) else format(z,'.15g'))
                frame['baseline_bmi_passive_status']=np.select([bmi.isna(),imputed_components],['components_missing_not_derived','derived_from_completed_components'],default='derived_from_observed_components')
                frame['baseline_bmi_passive_changed']=changed.astype(int)
                frame['gallstone_related_disease']=outcomes.gallstone_related_disease
                frame['.imp']=i
                frame.to_csv(handle,index=False,header=i==1,na_rep='')
                filled_draws.append(filled);passive_draws.append(int(changed.sum()))
        temporary.replace(target)
        summary[setting]={
            'rows_per_imputation':len(clean),'imputations':50,'long_rows':len(clean)*50,
            'eligible_missing_cells':eligible_cells,'released_cells_per_imputation':first_fill,
            'unreleased_cells':eligible_cells-first_fill,'unique_variables_with_released_cells':len(variables_filled),
            'unique_variables_flagged_unreliable':inventory.loc[inventory.status.ne('imputed'),'variable'].nunique(),
            'imputed_target_strata':int(inventory.status.eq('imputed').sum()),
            'unreliable_target_strata':int(inventory.status.ne('imputed').sum()),
            'passive_bmi_changed_per_draw_min':min(passive_draws),'passive_bmi_changed_per_draw_max':max(passive_draws),
            'source_columns_retained':len(clean.columns),'assembled_columns':len(frame.columns)}
        print(setting,json.dumps(summary[setting]),flush=True)
    rows='\n'.join(f"| {s} | {v['rows_per_imputation']} | {v['eligible_missing_cells']} | {v['released_cells_per_imputation']} | {v['unreleased_cells']} | {v['unique_variables_with_released_cells']} | {v['unique_variables_flagged_unreliable']} |" for s,v in summary.items())
    cohort_text=' '.join(f"{setting}: {v['first_presentations']} eerste presentaties, {v['confirmed_events']} bevestigde events, {v['nonevents']} non-events en {v['outcome_unresolved_or_excluded']} uitgesloten of onopgeloste uitkomsten; {v['prior_chole_excluded']} eerdere cholecystectomieën en {v['duplicate_identity_pending']} records met mogelijke dubbele identiteit." for setting,v in clean_summary['settings'].items())
    unassigned=len(read(OUT/'first_episodes_unassigned_setting.csv.gz'))
    text=f'''# OPERATE — eerste presentatie, versie 27 oktober 2026 v3

De opgeschoonde bestanden en **50 afzonderlijke imputaties per setting** zijn beschikbaar. Er is geen predictiemodel geschat. Imputaties zijn geen vastgestelde metingen: behoud alle 50 datasets bij latere analyse en combineer onzekerheid passend.

| Setting | Personen per imputatie | Toegestane ontbrekende cellen | Vrijgegeven ingevulde cellen per imputatie | Niet vrijgegeven cellen | Variabelen met ingevulde cellen | Variabelen met onbetrouwbare strata |
|---|---:|---:|---:|---:|---:|---:|
{rows}

{cohort_text} Daarnaast blijven {unassigned} eerste episodes zonder vastgestelde setting apart beschikbaar.

Een variabele kan in een betrouwbaar stratum wel en in een ander stratum niet zijn ingevuld. De tellingen gelden daarom op cel- en stratumniveau; niet alle lege cellen zijn geschikt voor imputatie.

## Bestanden en koppeling

- [SEH: volledig opgeschoond](cohort_SEH_cleaned.csv.gz) · [alle 50 imputaties](cohort_SEH_all_imputations.csv.gz) · [MICE-objecten en diagnostiek](mice_SEH.rds)
- [POLI: volledig opgeschoond](cohort_POLI_cleaned.csv.gz) · [alle 50 imputaties](cohort_POLI_all_imputations.csv.gz) · [MICE-objecten en diagnostiek](mice_POLI.rds)
- [SEH: oorspronkelijke ontbreekredenen](cohort_SEH_analysis_missing_status.csv.gz) · [POLI: oorspronkelijke ontbreekredenen](cohort_POLI_analysis_missing_status.csv.gz)
- [Afzonderlijke statistische validatie](statistical_validation.json) · [Geobserveerde versus geïmputeerde verdelingen](imputation_distribution_diagnostics.csv)
- [Klinische broncontrolewachtrij](clinical_review_queue.csv) · [Samenvatting broncontrole](clinical_review_summary.json)
- [SEH: veldbeoordeling](imputation_spec_SEH.csv) · [POLI: veldbeoordeling](imputation_spec_POLI.csv)
- [SEH: uitkomsten en uitsluitingen](outcomes_SEH.csv) · [POLI: uitkomsten en uitsluitingen](outcomes_POLI.csv)
- [SEH: vrijgave per variabele/stratum](imputation_inventory_SEH.csv) · [POLI: vrijgave per variabele/stratum](imputation_inventory_POLI.csv)
- [SEH: predictorselectie en technische uitsluitingen](imputation_predictor_audit_SEH.csv) · [POLI: predictorselectie en technische uitsluitingen](imputation_predictor_audit_POLI.csv)

`participant_id` is de unieke koppelkolom binnen een setting; in het lange bestand is `participant_id` samen met `.imp` uniek. `.imp` loopt van 1 tot en met 50. De oorspronkelijke bronkolommen blijven behouden; imputaties vullen uitsluitend toegestane lege analytische `*_analysis`-cellen. De aparte bestanden `imputed_changed_{{setting}}_{{01..50}}.csv.gz` bevatten celgewijze wijzigingsindicatoren. De ontbreekredenbestanden beschrijven de oorspronkelijke toestand en worden niet overschreven door imputatie.

`baseline_bmi_calculated` is in iedere dataset consequent opnieuw berekend uit `startgewicht_analysis / (lengte_analysis/100)^2`. `baseline_bmi_passive_status` maakt onderscheid tussen geobserveerde, voltooide en ontbrekende componenten; `baseline_bmi_passive_changed` markeert verandering ten opzichte van de opgeschoonde afleiding. De gerapporteerde BMI blijft in de oorspronkelijke kolommen behouden. BMI is niet zelfstandig stochastisch geïmputeerd.

## Methode en grenzen

MICE is afzonderlijk uitgevoerd voor SEH en POLI: 50 imputaties, 30 iteraties, zo nodig verlengd tot maximaal 100. Klinische formulierblokken met gelijke toepasbaarheid worden gezamenlijk geïmputeerd. Bekende en ontbrekende uitkomsten vormen aparte strata; uitkomsten worden nooit geïmputeerd. Hulpvariabelen volgen vooraf vastgelegde klinische/CRF-volgorde, met maximaal 10 predictorvrijheidsgraden en conservatieve ondersteuning-/rangcontroles. Er is niet op uitkomstassociatie geselecteerd. Niet-ondersteunde of onvoldoende geconvergeerde doelvariabelen blijven in de vrijgegeven bestanden leeg; diagnostiek en eventuele trekkingen blijven in de RDS-bestanden inspecteerbaar. Deze technische controles bewijzen geen afwezigheid van imputatiebias.

Dit is een voorbereidende MI-oplevering. De afzonderlijke formulier-/toepasbaarheidsblokken vormen geen enkel gezamenlijk model dat alle verbanden tussen onvolledige variabelen uit verschillende blokken bewaart. De bestanden zijn daarom niet zonder meer geschikt voor elk willekeurig predictiemodel. Na vaststelling van het predictiemodel zijn modelspecifieke imputatie en interne validatie nodig.

De uitkomstcodes zijn gecontroleerd met de aangeleverde Castor-metadata (2024) en de export/CRF van 2026. Patiënten met eerdere cholecystectomie zijn uitgesloten als case én controle, maar blijven in de volledige bestanden. Mogelijke dubbele identiteit, onzekere voorgeschiedenis en andere uitsluitingen blijven zichtbaar; een uniek export-ID bewijst geen unieke persoon. Eerste episodes zonder vastgestelde setting staan apart. Nachtelijke onderzoeken en onzekere verrichtings-/aanvraagstatus worden niet automatisch aan een indexbezoek toegeschreven. Structurele missings, niet-verrichte onderzoeken, onzekere afhankelijkheden en niet-geverifieerde vrije-tekstcodes worden niet geforceerd ingevuld.

De klinische broncontrolewachtrij bevat 2.059 niet-lege broncellen voor beoordeling, grotendeels door veldbrede blokkades. Dit zijn geen 2.059 bewezen fouten. Ook op zichzelf plausibele ASA-waarden vallen onder de veldbrede quarantaine zolang de schaal/codering niet is opgehelderd.

CSV: UTF-8, komma als scheidingsteken. Lees ID's als tekst en behandel alleen lege cellen als ontbrekend; de letterlijke tekst `NA` kan een bronwaarde zijn. Deze oplevering vormt stap 1; latere modelontwikkeling en interne validatie zijn niet uitgevoerd. Gebruik de bron-/analysewoordenboeken en de actuele afzonderlijke statistische beoordeling voor interpretatie.

Reproduceerbare scripts: 13 (import Castor-metadata), 14 (v3 voorbereiding), 10 (MICE), 15 (afzonderlijke validatie), 16 (deze samenvoeging) en 18 (klinische broncontrolewachtrij). September- en v2-bestanden blijven behouden. De eerste v3-uitvoer is ingetrokken en uitsluitend voor herleidbaarheid bewaard in `first_visit_2026-10-27_v3_initial_withdrawn`; gebruik die versie niet voor analyse.

Bij de aanvullende controle zijn de analytische ASA-classificatie, urine-soortelijk-gewicht en eGFR met de onopgeloste bronunit `kl/1.73m2` in quarantaine geplaatst: analytische waarden blijven leeg en worden niet geïmputeerd of als hulpvariabele gebruikt. Urine-pH buiten 0–14 is eveneens analytisch leeggemaakt en uitgesloten van imputatie. De oorspronkelijke bronwaarden blijven behouden. Er zijn geen vermoedelijke decimalen, schaalcorrecties of ASA-categorieën geraden.

Ververs na de definitieve validatie en eventuele aanvullende rapporten alleen het integriteitsmanifest met `python 16_package_first_visit_v3.py --refresh-manifest` vanuit de analysemap, met de vastgelegde Python-omgeving. Deze optie voert geen imputatie of samenvoeging opnieuw uit.
'''
    (OUT/'LEESMIJ.md').write_text(text,encoding='utf-8')
    cfg['clinical_imputation_run']=True
    cfg['cleaning_validation_imputation_flag_note']='cleaning_validation.json records preparation-stage status before MICE; delivery_status.json records completed release'
    (OUT/'configuration.json').write_text(json.dumps(cfg,indent=2)+'\n')
    status={'cleaned_files_available':True,'clinical_imputation_completed':True,'outcome_mapping_status':cfg['outcome_mapping_status'],'prediction_model_fitted':False,'settings':summary}
    (OUT/'delivery_status.json').write_text(json.dumps(status,indent=2)+'\n')
    refresh_manifest()
if __name__=='__main__':
    if '--refresh-manifest' in sys.argv:refresh_manifest()
    else:main()
