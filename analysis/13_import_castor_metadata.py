#!/usr/bin/env python3
"""Import explicit Castor 2024 metadata and validate its 2026 source crosswalk.

Reads only metadata workbook sheets (never Study results participant rows), the
2022 variable list, current CRF text and aggregate current column code coverage.
Exact source names take precedence. Aliases require independent historical SAV
question/section corroboration or an explicitly listed reviewed rename.
No numerical codes are inferred from option display order or prevalence.
"""
from pathlib import Path
from collections import defaultdict
import hashlib
import importlib.util
import json
import re
import unicodedata
import zipfile
import xml.etree.ElementTree as ET
import pandas as pd
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'OPERATE data'
XLSX = DATA / 'OPERATE_Excelexport_10062024.xlsx'
SOURCE = DATA / 'Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
TXT = DATA / 'Variabelen lijst (zonder dubbelen).txt'
DOCX = DATA / 'WErkdiagnoses_definitieve diagnoses.docx'
OUT = ROOT / 'results/first_visit_2026-10-27_v3/codebook'
V2 = ROOT / 'results/first_visit_2026-10-27/codebook'
MISSING = {'-99': 'Not done', '-98': 'Asked but unknown', '-97': 'Not asked',
           '-96': 'Not applicable', '-95': 'Measurement failed'}

# Explicit first-visit renames reviewed against the unchanged CRF forms16/17,
# 22/23/24. Keep checkbox suffixes verbatim and never join by display order.
REVIEWED_SUFFIXES = {
    'waar': 'buikpijn_waarO', 'waar_O': 'buikpijn_waar_anders',
    'waar_text': 'buikpijn_waar', 'soort': 'buikpijn_soortO',
    'soort_text': 'buikpijn_soort', 'hoelang': 'buikpijn_hoelangO',
    'hoelang_O': 'buikpijn_hoelang_anders', 'hoelang_text': 'buikpijn_hoelang',
    'uitstraling': 'buikpijn_uitstralingO', 'uitstraling_O': 'buikpijn_uitstraling_anders',
    'uitstraling_text': 'buikpijn_uitstraling', 'continu': 'buikpijn_COA',
    'continu_text': 'buikpijn_COA_info', 'houding': 'buikpijn_houdingO',
    'houding_O': 'buikpijn_houding_anders', 'houding_text': 'buikpijn_houding',
    'intake': 'buikpijn_SMIO', 'intake_text': 'buikpijn_SMI',
    'vervoerspijn': 'buikpijn_vervoerspijn', 'bewegingsdrang': 'buikpijn_Bewegingsdrang',
    'ACNES': 'buikpijn_ACNES_RTS', 'passageklachten': 'passageklachten2',
    'intakeglad': 'intake',
}
BASELINE_RENAMES = {
    'B_operatie_waar_a':'Waar_operatie_anders',
    'B_operatie_vensters_info':'Operatie_vensters_info',
    'B_operatie_bijzonderheden':'Operatie_bijzonderheden',
    'B_operatie_bijz_info':'Operatie_bijzonderheden_uitleg',
    'B_operatie_complicaties30dgn':'complicaties_operatie',
    'B_compl_30dgn_lekkage':'Compl_30dagen_lekkage',
    'B_compl_30dgn_bloeding':'Compl_30dagen_bloeding',
    'B_compl_30dgn_stenose':'Compl_30dagen_stenose',
    'B_compl_30dgn_overig':'Compl_30dagen_overig',
    'B_compl_30dgn_anders':'Compl_30dagen_overig_anders',
    'B_compl_30dgn_hoeveel':'Complicaties_hoeveel',
    'B_compl_30dgn_info':'Complicatie_omschrijving',
    'B_compl_30dgn_IC':'Compl_30dagen_IC',
    'B_compl_30dgn_chi':'Compli_chirugisch',
    'B_AO_30dgn':'AO_1ste_30dagen',
    'B_SEH_30dgn':'SEH_1ste_30dagen',
    'B_SEH_30dgn_info':'SEH_1ste_30dagen_omschrijving',
    'B_heropname_30dgn':'Heropname_1ste_30dagen',
    'B_heropname_30dgn_info':'Heropname_1ste_30dagen_omschrijving',
    'B_heroperatie_30dgn':'Heroperatie_1ste_30dagen',
    'B_heroperatie_30dgn_info':'Heroperatie_1ste_30dagen_omschrijving',
    'VG_overige_abdom_OKs':'VG_basis_abdom_operaties_overig',
    'VG_abdom_operatie_anders':'VG_basis_abdom_operaties_overig_anders',
}


def reviewed_rename(name, by_name):
    base, sep, option = name.partition('#')
    target = ''
    echo_match = re.fullmatch(r'Echo([1-5])(?:_(datum|cholecystitis))?', base)
    if echo_match:
        number = int(echo_match.group(1))
        finding = echo_match.group(2)
        target = 'Echo_basis' + ('_' + finding if finding else '')
        target += '' if number == 1 else '_' + str(number-1)
    if base in BASELINE_RENAMES:
        target = BASELINE_RENAMES[base]
    if base == 'POLI_1stecontact':
        target = 'POLI_1ste_studiecontact'
    if base.startswith('SEH_A_'):
        suffix = base[len('SEH_A_'):]
        target = 'SEH_A_buikpijn_' + suffix
        if suffix == 'bewegingsdrang':
            target = 'SEH_A_buikpijn_Bewegingsdrang'
    for prefix in ('POLI_AC_', 'POLI_A_'):
        if base.startswith(prefix):
            suffix = base[len(prefix):]
            target = prefix + REVIEWED_SUFFIXES.get(suffix, '[NO_REVIEWED_ALIAS]')
            if prefix == 'POLI_AC_':
                target = target.replace('buikpijn_COA', 'buikpijn_COAC')
                target = target.replace('buikpijn_ACNES_RTS', 'buikpijn_ACCNES_RTS')
                if suffix == 'soort_O':
                    target = prefix + 'buikpijn_soort_anders'
    for prefix, old in (('SEH_LO_', 'LO_SEH_'), ('POLI_LO_', 'LO_POLI_')):
        if base.startswith(prefix):
            suffix = base[len(prefix):]
            special = {'temp':'Temp','sat':'Sat','pols':'pols','ACNES':'tekenen_ACNES',
                       'ACNES2':'ACNES_gevoel','soepel':'SB','overige_info':'overige_informatie',
                       'overig':'overige_informatie'}
            target = old + special.get(suffix, suffix)
    if base == 'POLI_LO_percussie':
        target = '24.1.11|Percussie'
    if target:
        target += sep + option
    if target in by_name and len(by_name[target]) == 1:
        return target
    return ''


def norm(text):
    text = str(text or '').replace('Y es', 'Yes')
    text = re.sub(r'-\s*\n\s*', '', text)
    text = unicodedata.normalize('NFKD', text).encode('ascii', 'ignore').decode()
    return re.sub(r'[^a-z0-9]+', ' ', text.lower()).strip()


def sheet_records(workbook, name):
    rows = workbook[name].iter_rows(values_only=True)
    header = next(rows)
    return [{k: ('' if v is None else str(v)) for k, v in zip(header, row)} for row in rows]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    workbook = load_workbook(XLSX, read_only=True, data_only=True)
    metadata = sheet_records(workbook, 'Study variable list')
    options = sheet_records(workbook, 'Field options')
    workbook.close()
    grouped_options = defaultdict(dict)
    for row in options:
        group, value, label = row['Option group name'], row['Option value'], row['Option name']
        if value in grouped_options[group] and grouped_options[group][value] != label:
            raise ValueError(f'Conflicting explicit options in {group}')
        grouped_options[group][value] = label
    by_name = defaultdict(list)
    for row in metadata:
        by_name[row['Variable name']].append(row)
    old_list = pd.read_csv(TXT, sep='\t', encoding='utf-16', dtype=str, keep_default_na=False)
    old_names = set(old_list['Variable name'])
    current = pd.read_csv(SOURCE, sep=';', dtype=str, keep_default_na=False)
    pages = json.loads((V2 / 'crf_pages.json').read_text())
    page_text = [p['text'] for p in pages]
    page_norm = [norm(p) for p in page_text]
    form_starts = {}
    for m in metadata:
        key = (m['Form number'], m['Form name'])
        if key in form_starts:
            continue
        needle = norm(m['Form number']+' '+m['Form name'])
        matches = [i for i in range(5, len(page_norm)) if needle and needle in page_norm[i]]
        if matches:
            form_starts[key] = matches[0]
    historical = json.loads((V2 / 'historical_value_labels.json').read_text())
    sav = historical['OPERATE_algoritme200924.sav']

    # Reviewed prefix renames are additionally checked against the same SAV
    # question label and section; a matching prefix alone never accepts an alias.
    aliases = {}
    for source_name, item in sav.items():
        if source_name not in current or source_name in by_name:
            continue
        question = item['label'].split('|', 1)[-1]
        number = item['label'].split('|', 1)[0]
        form = number.split('.')[0]
        candidates = [r for r in metadata
                      if norm(r['Field label']) == norm(question)
                      and r['Form number'] == form]
        # Several expanded checkboxes may share labels: never collapse these.
        if len(candidates) == 1:
            aliases[source_name] = (candidates[0], 'historical_SAV_same_question_label_and_form_number')

    rows, mapping, conflicts = [], {}, []
    for name in current.columns:
        candidates = by_name.get(name, [])
        exact = bool(candidates)
        evidence = 'exact_2024_metadata_variable_name' if exact else ''
        if len(candidates) > 1:
            # Duplicate definitions must be fully identical to remain usable.
            unique = {json.dumps(r, sort_keys=True) for r in candidates}
            if len(unique) != 1:
                conflicts.append(dict(source_name=name, issue='conflicting_duplicate_metadata_names'))
                candidates = []
        renamed = reviewed_rename(name, by_name)
        # POLI intake was renamed to intakeglad; current intake now denotes the
        # separate pain/intake question, so exact old-name matching is unsafe.
        force_alias = name in ('POLI_A_intake', 'POLI_AC_intake')
        if renamed and (not candidates or force_alias):
            m = by_name[renamed][0]
            evidence = 'explicit_reviewed_first_visit_rename_current_CRF_form_question'
        elif candidates:
            m = candidates[0]
        elif name in aliases:
            m, evidence = aliases[name]
        else:
            m = {}
        question = m.get('Field label', '')
        group = m.get('Optiongroup name', '')
        typ = m.get('Exported field type', 'unresolved')
        original = m.get('Original field type', '')
        raw_codes = dict(grouped_options.get(group, {}))
        option_label = name.split('#', 1)[1] if '#' in name else ''
        is_checkbox = typ == 'checkbox' and original == 'checkbox' and bool(option_label)
        option_codes = [code for code, label in raw_codes.items() if norm(label) == norm(option_label)]
        codes = {'0': 'Unselected', '1': 'Selected'} if is_checkbox and len(option_codes) == 1 else raw_codes
        categorical = typ in ('radio', 'dropdown', 'checkbox')
        values = current[name].str.strip()
        nonempty = values[values != '']
        observed_set = set(nonempty.unique()) if categorical else set()
        unexpected = observed_set - set(codes) - set(MISSING) if categorical else set()
        # Export only numerical tokens/counts; never clinical narrative.
        safe_unexpected = sorted(v for v in unexpected if re.fullmatch(r'-?\d+(?:\.\d+)?', v))
        nonnumeric_n = int(sum(not re.fullmatch(r'-?\d+(?:\.\d+)?', v) for v in unexpected))
        if nonnumeric_n:
            safe_unexpected.append('[NONNUMERIC_VALUES_REDACTED]')
        coverage = bool(m) and not unexpected
        form_start = form_starts.get((m.get('Form number',''), m.get('Form name','')))
        later_starts = sorted(v for v in form_starts.values() if form_start is not None and v > form_start)
        form_end = later_starts[0] if later_starts else len(page_norm)
        search_pages = range(form_start, form_end) if form_start is not None else range(len(page_norm))
        question_pages = [i+1 for i in search_pages if norm(question) and
                          re.search(r'\b'+re.escape(norm(question))+r'\b', page_norm[i])]
        unit = m.get('Field units','')
        unit_confirmed = not unit or any(norm(unit) in page_norm[i-1] for i in question_pages)
        if len(norm(question)) < 5 and (form_start is None or not unit_confirmed):
            question_pages = []
        # Explicit option-group appendix gives the complete current option set.
        # Group blocks end at the next unbulleted group title in the appendix.
        current_group_context = ''
        group_pages = []
        if group:
            for i in range(848, len(page_text)):
                lines = page_text[i].splitlines()
                for j, line in enumerate(lines):
                    if norm(line) == norm(group):
                        block = []
                        for line2 in lines[j+1:]:
                            if line2.startswith('Exported at'):
                                break
                            if line2.startswith('•') or (block and not line2.strip()):
                                block.append(line2)
                            elif block and not line2.startswith('•'):
                                # Wrapped option labels are retained until a known group title.
                                if norm(line2) in {norm(g) for g in grouped_options}:
                                    break
                                block.append(line2)
                            else:
                                break
                        current_group_context = '\n'.join(block)
                        group_pages.append(i+1)
        # Inline options checked on question-containing pages; workbook links
        # question to this exact group, rather than deriving codes from the CRF.
        context_norm = norm(current_group_context) if current_group_context else ' '.join(page_norm[p-1] for p in question_pages)
        options_confirmed = bool(raw_codes) and all(norm(label) in context_norm for label in raw_codes.values())
        if is_checkbox:
            options_confirmed = len(option_codes) == 1 and norm(option_label) in context_norm
        question_confirmed = bool(question_pages)
        verified = bool(m and coverage and question_confirmed and (not categorical or options_confirmed))
        state = 'verified_versioned_crosswalk' if verified else 'provisional' if m else 'unresolved'
        parent = m.get('Dependency parent name', '')
        parent_rows = [r for r in metadata if parent and r['Field label'] == parent
                       and r['Visit name'] == m.get('Visit name')]
        same_form = [r for r in parent_rows if r['Form number'] == m.get('Form number')]
        parent_candidates = [r['Variable name'] for r in (same_form or parent_rows)]
        parent_candidates = sorted(set(parent_candidates))
        record = dict(source_name=name, metadata_variable_name=m.get('Variable name', ''),
                      question=question, label=question, type=typ,
                      original_type=original, exported_type=typ,
                      visit=m.get('Visit name', ''), form=m.get('Form name', ''),
                      form_number=m.get('Form number', ''), field_number=m.get('Field number', ''),
                      option_group=group, codes=codes, codes_json=json.dumps(codes, ensure_ascii=False),
                      option_label=option_label, selected_option_group_code=option_codes[0] if len(option_codes)==1 else '',
                      missing_reason_codes=MISSING,
                      units=m.get('Field units', ''), required=m.get('Field required', ''),
                      minimum=m.get('Field min', ''), maximum=m.get('Field max', ''),
                      info=m.get('Field info', ''), calculation_template=m.get('Calculation/summary template', ''),
                      hidden=m.get('Field hidden', ''), dependency_parent=parent,
                      dependency_operator=m.get('Field dependency operator', ''),
                      dependency_value=m.get('Field dependency value', ''),
                      dependency_parent_metadata_candidates=parent_candidates,
                      verification_state=state, verified_crosswalk=verified,
                      current_code_coverage=coverage, unexpected_codes=safe_unexpected,
                      nonnumeric_unexpected_distinct_count=nonnumeric_n,
                      current_question_confirmed=question_confirmed,
                      current_units_confirmed=unit_confirmed,
                      current_options_confirmed=options_confirmed,
                      current_crf_pages=question_pages, current_crf_option_group_pages=group_pages,
                      current_crf_option_group_context=current_group_context,
                      evidence=evidence,
                      metadata_source_version='Castor Excel export 2024-06-10; crosschecked current CSV2026-06-22 / CRF2026-08-11',
                      name_present_in_2022_list=name in old_names,
                      metadata_name_present_in_2022_list=m.get('Variable name','') in old_names,
                      checkbox_missing_policy='0/1 encode option selection only; blank and untriggered remain missing/structural' if is_checkbox else '')
        rows.append(record)
        mapping[name] = record
    metadata_to_source = defaultdict(list)
    for row in rows:
        if row['metadata_variable_name']:
            metadata_to_source[row['metadata_variable_name']].append(row['source_name'])
    for row in rows:
        sources = sorted(set(source for old in row['dependency_parent_metadata_candidates']
                             for source in metadata_to_source.get(old, [])))
        row['dependency_parent_source_candidates'] = sources
        row['dependency_parent_source'] = sources[0] if len(sources) == 1 else ''
    critical_checks = {
        'FU_einde_werkdiagnose': {'6':'Cholecystolithiasis','7':'Cholecystitis'},
        'FU_einde_werkdiagnose_1_DEF': {'5':'Cholecystolithiasis','6':'Cholecystitis'},
        'FU_einde_diagnose_definitief_ordinaal': {'1':'Symptomatisch galsteenlijden','3':'Cholecystitis'},
        'FU_einde_definitieve_diagnose_1_DEF_1': {'5':'Cholecystolithiasis','6':'Cholecystitis'},
        'FU_einde_diagnose_definitief_ordinaal_2': {'5':'Cholecystolithiasis','6':'Cholecystitis'},
        'FU_einde_definitieve_diagnose_3_ordinaal': {'5':'Cholecystolithiasis','6':'Cholecystitis'},
        'Her_op_soort': {'2':'Laparoscopische cholecystectomie','3':'Laparoscopische cholecystolithiotomie'},
        'Her_op_soort_1': {'2':'Laparoscopische cholecystectomie','3':'Laparoscopische cholecystolithiotomie'},
        'Her_op_soort_2': {'2':'Laparoscopische cholecystectomie','3':'Laparoscopische cholecystolithiotomie'},
        'VG_cholecystectomie': {'1':'Nee','2':'Ja (manier onbekend/onzeker)','3':'Ja laparoscopisch','4':'Ja klassiek'},
    }
    for number in range(1,6):
        critical_checks[f'Echo{number}_cholecystitis'] = {'1':'Yes','0':'No'}
    for number in range(3):
        suffix = '' if number == 0 else '_' + str(number)
        critical_checks['Her_op'+suffix] = {'1':'Yes','0':'No'}
        diagnosis = ['Diagnose_heroperatie_definitief_ordinaal',
                     'Diagnose_heroperatie_2_definitief_ordinaal',
                     'Diagnose_heroperatie_3_definitief_ordinaal'][number]
        critical_checks[diagnosis] = {'6':'Cholecystitis'}
    for name, expected in critical_checks.items():
        if not mapping[name]['verified_crosswalk'] or any(mapping[name]['codes'].get(k) != v for k,v in expected.items()):
            raise ValueError(f'Critical diagnosis/surgery code validation failed: {name}')
    evidence_fields = []
    for number in range(1,6):
        evidence_fields.extend((f'Echo{number}_datum',f'Echo{number}_cholecystitis'))
    for number in range(3):
        suffix = '' if number == 0 else '_' + str(number)
        evidence_fields.extend(['Her_op'+suffix,'Her_op_datum'+suffix,'Her_op_soort'+suffix,
                                'Her_Op_overig'+suffix,
                                ['Diagnose_heroperatie_definitief_ordinaal','Diagnose_heroperatie_2_definitief_ordinaal',
                                 'Diagnose_heroperatie_3_definitief_ordinaal'][number],
                                ['Diagnose_heroperatie_definitief','Diagnose_heroperatie_definitief_1',
                                 'Diagnose_heroperatie_definitief_2'][number]])
    if not all(mapping[name]['verified_crosswalk'] for name in evidence_fields):
        raise ValueError('Unverified field in ultrasound/surgery outcome evidence')
    (OUT / 'mappings.json').write_text(json.dumps(mapping, ensure_ascii=False, indent=2))
    csv_rows = [{k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict,list)) else v
                 for k,v in row.items()} for row in rows]
    pd.DataFrame(csv_rows).to_csv(OUT / 'field_codebook.csv', index=False)
    pd.DataFrame([{k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v
                   for k,v in mapping[name].items()} for name in evidence_fields]).to_csv(
        OUT / 'outcome_evidence_metadata_links.csv', index=False)
    pd.DataFrame(metadata).to_csv(OUT / 'study_variable_list_2024.csv', index=False)
    pd.DataFrame(options).to_csv(OUT / 'field_options_2024.csv', index=False)
    old_list.to_csv(OUT / 'variable_list_2022.csv', index=False)
    pd.DataFrame([r for r in csv_rows if r['evidence'] and r['source_name'] != r['metadata_variable_name']]).to_csv(OUT / 'documented_alias_crosswalk.csv', index=False)
    with zipfile.ZipFile(DOCX) as archive:
        doc = ET.fromstring(archive.read('word/document.xml'))
        doc_text = '\n'.join(''.join(p.itertext()) for p in doc.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}p'))
    (OUT / 'diagnosis_document_text.txt').write_text(doc_text)
    summary = dict(current_source_shape=list(current.shape), metadata_fields=len(metadata),
                   exact_name_matches=sum(r['evidence']=='exact_2024_metadata_variable_name' for r in rows),
                   documented_aliases=sum(bool(r['evidence']) and r['source_name'] != r['metadata_variable_name'] for r in rows),
                   states=pd.Series([r['verification_state'] for r in rows]).value_counts().to_dict(),
                   checkbox_fields=sum(r['type']=='checkbox' for r in rows),
                   metadata_conflicts=conflicts,
                   critical_diagnosis_surgery_mapping_assertions_passed=list(critical_checks),
                   source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [XLSX, TXT, DOCX, SOURCE]},
                   limitations=['Metadata source version is 2024; 2026 observations retained unchanged.',
                                'No Study results participant rows read.',
                                'Blank checkbox values are never filled as unselected.',
                                'Dependency parents are labels; multiple candidate variables require contextual resolution.'])
    (OUT / 'validation_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('source_sha256','metadata_conflicts')}))


if __name__ == '__main__':
    main()
