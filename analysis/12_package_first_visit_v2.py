#!/usr/bin/env python3
"""Human handover and optional full-cohort assembly after validated clinical MI.

Does not create imputation placeholders. If MICE has not run, accurately report
the blocked stage. Later completed data are joined one-to-one and only eligible
originally empty analytical cells may change. BMI is passively recalculated.
"""
from pathlib import Path
import gzip
import hashlib
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/first_visit_2026-10-27'

def read(p): return pd.read_csv(p,dtype=str,keep_default_na=False,na_filter=False)

def main():
    config=json.loads((OUT/'configuration.json').read_text())
    summary=json.loads((OUT/'cleaning_validation.json').read_text())
    completed=False
    if config['outcome_mapping_status']=='validated':
        expected=[OUT/f'imputed_{s}_{i:02d}.csv.gz' for s in ['SEH','POLI'] for i in range(1,51)]
        if all(p.exists() for p in expected):
            for setting in ['SEH','POLI']:
                clean=read(OUT/f'cohort_{setting}_cleaned.csv.gz')
                where=read(OUT/f'where_{setting}.csv.gz')
                outcomes=read(OUT/f'outcomes_{setting}_review.csv')
                assert clean.participant_id.equals(outcomes.participant_id)
                with gzip.open(OUT/f'cohort_{setting}_all_imputations.csv.gz','wt',encoding='utf-8',newline='') as handle:
                    for i in range(1,51):
                        values=read(OUT/f'imputed_{setting}_{i:02d}.csv.gz')
                        assert clean.participant_id.equals(values.participant_id)
                        frame=clean.copy()
                        for col in values:
                            if col not in frame or col=='participant_id': continue
                            use=clean[col].eq('') & values[col].ne('')
                            if use.any():
                                assert col in where
                                assert where.loc[use,col].str.lower().isin(['true','1']).all()
                                frame.loc[use,col]=values.loc[use,col]
                        h=pd.to_numeric(frame.lengte_numeric,errors='coerce')/100
                        w=pd.to_numeric(frame.startgewicht_numeric,errors='coerce')
                        missing_bmi=clean.baseline_bmi_calculated.eq('')
                        frame.loc[missing_bmi,'baseline_bmi_calculated']=(w/h.pow(2)).where(missing_bmi).map(lambda x:'' if pd.isna(x) else str(x)).loc[missing_bmi]
                        frame['gallstone_related_disease']=outcomes.gallstone_related_disease
                        frame['.imp']=i
                        frame.to_csv(handle,index=False,header=i==1,na_rep='')
            completed=True
    summaryrows='\n'.join(f"| {s} | {x['first_presentations']} | {x['columns']} | {x['prior_cholecystectomy_excluded_both_groups']} | {x['identity_review_required']} |" for s,x in summary['settings'].items())
    status=('Klinische imputaties beschikbaar; controleer de veld- en celrapportages voor niet-imputeerbare variabelen.' if completed else
        '**De opgeschoonde bestanden zijn beschikbaar. Klinische multiple imputation is nog niet uitgevoerd.** De huidige CSV bevat geen code–antwoordwoordenboek. De definitieve-diagnose- en procedurecodes kunnen nog niet voldoende worden bevestigd; de imputatie met uitkomst als hulpvariabele wacht hierop. Er zijn geen fictieve imputatiebestanden gemaakt.')
    text=f'''# OPERATE — eerste presentatie, versie 27 oktober 2026

{status}

## Beschikbare gegevens

| Cohort | Eerste presentaties | Kolommen in opgeschoond bestand | Eerdere cholecystectomie: uitgesloten als case én controle | Mogelijke dubbele identiteit: broncontrole |
|---|---:|---:|---:|---:|
{summaryrows}

De volledige bestanden behouden deze patiënten. De 28 eerste episodes zonder vastgestelde setting blijven apart beschikbaar. Een uniek Castor-recordnummer bewijst niet automatisch een unieke persoon. Er zijn geen mogelijke dubbelen op basis van alleen demografische overeenkomsten verwijderd.

- [SEH: volledig opgeschoond bestand](cohort_SEH_cleaned.csv)
- [Poli: volledig opgeschoond bestand](cohort_POLI_cleaned.csv)
- [Beschikbaarheid van alle variabelen](availability_all_variables.csv)
- [SEH: variabelenwoordenboek](variable_dictionary_SEH.csv)
- [Poli: variabelenwoordenboek](variable_dictionary_POLI.csv)
- [SEH: voorlopige uitkomstbeoordeling](outcomes_SEH_review.csv)
- [Poli: voorlopige uitkomstbeoordeling](outcomes_POLI_review.csv)

## Wat is veranderd?

Alle oorspronkelijke eerste-bezoekvelden en waarden uit de gecontroleerde septemberversie zijn behouden. Toegevoegd zijn gelabelde en analytische variabelen voor gericht gecontroleerde categorische vragen, een berekende baseline-BMI, expliciete uitsluitingsmarkeringen, beschikbaarheidsrapporten en een celgewijze imputatiespecificatie. De ruwe bron en septemberresultaten zijn ongewijzigd.

Een cholecystectomie vóór presentatie sluit een patiënt uit als zowel case als controle. Een onduidelijke diagnose bij een beschikbare follow-upbeoordeling mag volgens het afgesproken eindpunt als 0 gelden: geen definitief galsteenlijden vastgesteld. Een werkelijk ontbrekende uitkomstbeoordeling wordt niet gelijkgesteld aan 0.

Follow-up gebruikt de geregistreerde afsluiting van CRF134 na drie maanden zonder contact, niet indexdatum + 90 dagen. Bevestigende onderzoeken en behandelingen tellen alleen binnen die episode. Chronologische of inhoudelijke tegenstrijdigheden blijven gemarkeerd. De `candidate_endpoint_*`-kolommen zijn uitsluitend voor controle; `gallstone_related_disease` blijft leeg tot de uitkomstcoderingen bevestigd zijn. Gebruik de voorlopige tellingen niet als definitieve onderzoeksresultaten.

## Wat moet nog worden bevestigd?

1. **De actuele antwoordcodes en velddefinities uit Castor.** De lege CRF toont vragen en antwoordteksten, niet hun opgeslagen numerieke codes. Historische SPSS-metadata zijn gedeeltelijk bruikbaar: 17 veldspecifieke vraag-/optiecontroles zijn uitgevoerd, met expliciete historische provenance. Bij zes andere velden bevat de huidige export vrije tekst waar vroeger categoriecodes stonden. Definitieve-diagnosevelden hebben bovendien andere namen. Daarom is automatische overdracht van alle oude labels onjuist.
2. **Patiëntidentiteit:** vier mogelijke dubbele paren en de onbekende setting van 28 eerste episodes. Zonder betrouwbare bronkoppeling wordt niets samengevoegd of verzonnen.
3. **Onderzoeksmoment en verrichting:** {summary['overnight_candidates']} onderzoekrecords op de volgende kalenderdag zijn als mogelijke nachtelijke voortzetting gemarkeerd. De export bevat geen bruikbare ontslagtijd of eenduidige bezoek-ID waarmee deze automatisch aan het indexbezoek kunnen worden toegewezen. Een lege labwaarde bewijst evenmin dat die specifieke bepaling was aangevraagd.

## Imputatie

De R-procedure is geïmplementeerd voor 50 imputaties per setting, 30 iteraties en verlenging tot maximaal 100 bij onvoldoende convergentie. Bescherming van oorspronkelijke waarden, structurele missings, uitkomst en ID is ingebouwd. PMM en binaire, nominale en ordinale methoden zijn op synthetische gegevens getest. Die testgegevens zijn **geen patiëntresultaten**.

De invoermatrices en een beoordeling van alle velden staan klaar. Niet-geverifieerde codes, vrije tekst, niet-verrichte onderzoeken en onzekere afhankelijkheden worden niet ingevuld. Een analyseklaar bestand hoeft daarom niet overal waarden te bevatten. Zodra de noodzakelijke uitkomstmetadata bevestigd zijn, wordt dezelfde procedure op de klinische gegevens uitgevoerd; per variabele wordt dan getoond wat wel en niet betrouwbaar kon worden geïmputeerd.

## Herleidbaarheid en gebruik

CSV-bestanden gebruiken UTF-8 en komma als scheidingsteken. Lees ze als tekst zonder automatische omzetting van de letterlijke tekst `NA` naar een ontbrekende waarde. Gebruik het woordenboek voor numerieke/categorische conversie; gebruik de bijbehorende missing-statusbestanden om de reden van een lege cel te bepalen.

Reproduceerbare verwerking: scripts 08 (codeboekcontrole), 09 (opschoning), 10 (MICE, met controle op uitkomstmetadata), 11 (onafhankelijke validatie) en 12 (oplevering). Het MICE-pakket is in de eigen OPERATE R-bibliotheek geïnstalleerd. Er is geen predictiemodel geschat.
'''
    (OUT/'LEESMIJ.md').write_text(text)
    (OUT/'delivery_status.json').write_text(json.dumps({'cleaned_files_available':True,'clinical_imputation_completed':completed,
        'outcome_mapping_status':config['outcome_mapping_status'],'prediction_model_fitted':False},indent=2)+'\n')
    manifest={str(p.relative_to(OUT)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
        for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='package_manifest.json'}
    (OUT/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'files':len(manifest),'clinical_imputation_completed':completed}))

if __name__=='__main__': main()
