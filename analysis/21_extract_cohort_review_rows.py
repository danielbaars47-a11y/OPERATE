#!/usr/bin/env python3
"""Extract protected colleague-review rows; never change cohorts or endpoints.

One row per ID/issue, except next-day linkage retains one row per test instance.
Raw source strings are preserved; verified metadata labels are appended. No names,
birth dates or unrelated clinical histories are exported. No workbook is authored.
"""
from pathlib import Path
import importlib.util
import json
import re
from collections import Counter
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V3 = ROOT / 'results/first_visit_2026-10-27_v3'
OUT = ROOT / 'results/review_overview_2026-10-27/cohort_review_rows.json'
SOURCE = ROOT / 'OPERATE data/Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
PAIR_SOURCE = ROOT / 'results/cleaned_2026-09-30/candidate_duplicate_first_episodes.csv'
IDENTITY = 'possible_duplicate_person_identity_unresolved'
PRIOR = 'prior_cholecystectomy_excluded_both_case_and_control'
DX = ['FU_einde_diagnose','FU_einde_diagnose_definitief_ordinaal',
      'FU_einde_definitieve_diagnose_1_DEF_1','FU_einde_diagnose_2',
      'FU_einde_diagnose_definitief_ordinaal_2','FU_einde_definitieve_diagnose_3',
      'FU_einde_definitieve_diagnose_3_ordinaal']
TEXT_DX = ['FU_einde_diagnose_definitief','FU_einde_def_diagnose_overig_1',
           'FU_einde_def_diagnose_overig_2','FU_einde_definitieve_diagnose_3_overig']
DATES = ['FU_einde_datum','definitieve_diagnose_datum','definitieve_diagnose_datum_2','Definiteve_diagnose_datum_3']
QUESTIONS = {
 'outcome_assessment_missing':'Wilt u de definitieve uitkomstbeoordeling aanvullen op basis van het dossier? Is er een diagnose gesteld, welke, wanneer en op welk bewijs?',
 'outcome_chronology_or_source_conflict':'Welke bezoek-, follow-up- en diagnosedata en diagnosecodes zijn correct? Wilt u de tegenstrijdigheid met het dossier oplossen?',
 'other_possible_gallstone_diagnosis_review':'Wat is de definitieve diagnose en het onderliggende bewijs? Valt deze diagnose volgens de afgesproken definitie onder het eindpunt?',
 'definitive_label_without_verifiable_in_episode_evidence':'Is er tussen het indexbezoek en einde follow-up gedateerd kwalificerend echo-, CT-, operatie- of ERCP-bewijs? Zo ja: welk bewijs en welke datum?',
 'prior_cholecystectomy_history_unresolved':'Heeft deze deelnemer vóór het indexbezoek een cholecystectomie gehad? Graag status, datum en bron bevestigen; niet afleiden of imputeren.',
 IDENTITY:'Zijn dit dezelfde persoon of verschillende personen? Welke registratie is de eerste episode en hoe moeten de records gekoppeld of behouden worden?',
 'first_presentation_setting_unresolved':'Wat was het eerste studiecontact: SEH of polikliniek, en op welke datum? Graag bevestigen met het oorspronkelijke bezoekverslag.',
 'overnight_next_day_linkage':'Hoort dit onderzoek van de volgende kalenderdag bij hetzelfde SEH-/polikliniekbezoek? Graag bezoekbegin/einde of expliciete episodekoppeling én beschikbaarheid op het voorspeltijdstip bevestigen.'}


def read(path):
    return pd.read_csv(path,dtype=str,keep_default_na=False)


def main():
    raw = pd.read_csv(SOURCE,sep=';',dtype=str,keep_default_na=False)
    if not raw['Participant Id'].is_unique:
        raise ValueError('Source participant IDs not unique')
    raw = raw.set_index('Participant Id',drop=False)
    metadata = json.loads((V3/'codebook/mappings.json').read_text())
    outcomes = {s:read(V3/f'outcomes_{s}.csv') for s in ('SEH','POLI')}
    setting_by_id = {r.participant_id:s for s,d in outcomes.items() for _,r in d.iterrows()}
    prior_ids = {r.participant_id for d in outcomes.values() for _,r in d.iterrows()
                 if r.prior_cholecystectomy_excluded=='1'}
    unassigned = read(V3/'first_episodes_unassigned_setting.csv.gz')
    first_ids = set(setting_by_id)|set(unassigned.participant_id)
    prior_ids.update(pid for pid in unassigned.participant_id
                     if raw.at[pid,'VG_cholecystectomie'].strip() in ('2','3','4'))
    rows = []

    def value(pid,col):
        if col not in raw:
            raise ValueError(f'Unknown raw field {col}')
        v = raw.at[pid,col]
        if v == '':
            return '[leeg]'
        m = metadata.get(col,{})
        label = m.get('codes',{}).get(v) if m.get('verified_crosswalk') else None
        if not label and m.get('verified_crosswalk'):
            label = m.get('missing_reason_codes',{}).get(v)
        return v + (f' [{label}]' if label else '')

    def add(pid,setting,issue,fields,context='',extra='',priority='Dossiercontrole',source=''):
        if pid not in first_ids:
            raise ValueError('Review record outside the saved first-episode cohort')
        fields = list(dict.fromkeys(fields))
        original = '\n'.join(f'{c}: {value(pid,c)}' for c in fields)
        if extra:
            original += '\n'+extra
        rows.append(dict(study_id=pid,setting=setting,issue=issue,variables='; '.join(fields),
                         original_value=original,context=context,question=QUESTIONS[issue],
                         source_reference=str(SOURCE.relative_to(ROOT))+'; '+source,
                         priority=priority,existing_exclusion='Ja: eerdere cholecystectomie, bestaande protocoluitsluiting' if pid in prior_ids else 'Nee: geen reeds vastgestelde uitsluiting wegens eerdere cholecystectomie'))

    for setting,d in outcomes.items():
        for _,r in d.iterrows():
            issue = r.review_reason
            if issue in ('',PRIOR,IDENTITY):
                continue
            fields = [f'{setting}_datum']+DATES
            if issue=='prior_cholecystectomy_history_unresolved':
                fields = [f'{setting}_datum','VG_cholecystectomie']
            elif issue=='other_possible_gallstone_diagnosis_review':
                fields += DX+TEXT_DX
            else:
                fields += DX+TEXT_DX
            context = 'Huidige v3-status: niet ingedeeld als uitkomst 0/1; deze rij vraagt broncontrole, geen automatische correctie.'
            if issue=='definitive_label_without_verifiable_in_episode_evidence':
                fields += [c for c in raw if re.fullmatch(r'(?:Echo[1-5]_(?:datum|cholecystitis)|Her_op(?:_datum|_soort)?(?:_[12])?|AO_overig_ERCP_(?:datum|afwijkend|diagnose)(?:_[123])?)',c)]
                context += ' Relevante datum-/procedure-/bevindingcodes zijn toegevoegd; volledige verslagen blijven in het dossier.'
            add(r.participant_id,setting,issue,fields,context,source=f'v3/outcomes_{setting}.csv: review_reason')

    pairs = read(PAIR_SOURCE)
    identity_seen = set()
    episode_fields = ['episodes','Episode_hoeveel']+[c for c in raw if c.startswith('Episode_nummer_')]
    for _,pair in pairs.iterrows():
        for own,other in [('a','b'),('b','a')]:
            pid,partner = pair[f'participant_id_{own}'],pair[f'participant_id_{other}']
            if pid in identity_seen:
                raise ValueError('ID belongs to more than one candidate pair; aggregate partners before exporting')
            identity_seen.add(pid)
            setting = setting_by_id.get(pid,'unassigned')
            fields = episode_fields+['SEH_1stecontact','SEH_datum','POLI_1stecontact','POLI_datum']
            add(pid,setting,IDENTITY,fields,
                context=f'Kandidaatpaar {pair.candidate_group_id}; gekoppeld kandidaat-studie-ID: {partner}. Identiteit niet bevestigd; niet automatisch samenvoegen of verwijderen.',
                extra=f'Kandidaat partner studie-ID: {partner}; paarregistratie: {pair.candidate_group_id}',
                source=str(PAIR_SOURCE.relative_to(ROOT)))

    for _,r in unassigned.iterrows():
        add(r.participant_id,'unassigned','first_presentation_setting_unresolved',
            ['episodes','SEH_1stecontact','SEH_datum','POLI_1stecontact','POLI_datum'],
            'Expliciete eerste episode, maar eerste zorgsetting niet bevestigd. Twee lege velden of twee No-antwoorden bepalen de setting niet.',
            source='v3/first_episodes_unassigned_setting.csv.gz')

    # Import helper definitions only; use static source-column group descriptors.
    spec = importlib.util.spec_from_file_location('review_group_definitions',ROOT/'analysis/03_prepare_operate_first_visit.py')
    module = importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    groups = {(mod,g['selector']):g for mod,gs in module.build_groups(list(raw.columns)).items() for g in gs}
    overnight = read(V3/'overnight_visit_linkage_review.csv')
    for _,r in overnight.iterrows():
        group = groups[(r.modality,r.source_instance)]
        result_fields = []
        for _,field in group['mapping']:
            m = metadata.get(field,{})
            if field in (group['selector'],group['date']):
                continue
            if m.get('type') in ('string','textarea') or re.search(r'verslag|informatie|_info|naam|arts|chirurg|bevindingen',field,re.I):
                continue
            raw_value = raw.at[r.participant_id,field]
            if m.get('type') not in ('radio','checkbox','dropdown','date','time') and not re.fullmatch(
                    r'\s*(?:[<>]=?\s*)?[+-]?(?:\d+(?:[.,]\d*)?|[.,]\d+)\s*',raw_value):
                continue
            if raw_value != '':
                result_fields.append(field)
        fields = [f'{r.setting}_datum',group['selector'],group['date']]+result_fields[:35]
        context = (f'Onderzoek: {r.modality}; broninstantie: {r.source_instance}; indexdatum: {r.index_date}; '
                   f'onderzoeksdatum: {r.test_date}. Volgende kalenderdag; koppeling en bezoekeinde nog niet bevestigd. '
                   'Geen automatische opname of imputatie; uitsluitend aanwezige compacte gestructureerde uitslagen toegevoegd.')
        if len(result_fields)>35:
            context += f' Overige gestructureerde uitslagvelden ({len(result_fields)-35}) staan in de bron.'
        add(r.participant_id,r.setting,'overnight_next_day_linkage',fields,context,
            priority='Bezoekkoppeling',source='v3/overnight_visit_linkage_review.csv; '+r.source_instance)

    counts = Counter(r['issue'] for r in rows)
    if len(identity_seen)!=8 or len(unassigned)!=28 or len(overnight)!=428 or len(rows)!=507:
        raise ValueError(f'Unexpected review row counts: {dict(counts)}')
    nonovernight = [(r['study_id'],r['issue']) for r in rows if r['priority']!='Bezoekkoppeling']
    if len(nonovernight)!=len(set(nonovernight)):
        raise ValueError('Duplicate ID/issue review rows')
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rows,ensure_ascii=False,indent=2))
    print(json.dumps(dict(rows=len(rows),unique_study_ids=len({r['study_id'] for r in rows}),
                         issue_counts=dict(counts),priority_counts=dict(Counter(r['priority'] for r in rows))),ensure_ascii=False))


if __name__=='__main__':
    main()
