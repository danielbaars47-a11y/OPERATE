#!/usr/bin/env python3
"""Audit historical SPSS dictionaries against the current CSV and printed CRF.

Reads dictionary records only, never SPSS participant records. Historical numeric
codes remain provisional because the current printed CRF does not expose codes.
No labels are inferred from option order, prevalence, or patient values. Outputs
contain only aggregate code counts for labelled categorical fields, never IDs.
"""
from pathlib import Path
import hashlib
import json
import re
import struct
import unicodedata
import pandas as pd
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/first_visit_2026-10-27/codebook'
SOURCE = ROOT / 'OPERATE data/Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
PDF = ROOT / 'OPERATE data/empty_crf.pdf'
SAVS = [ROOT / 'Predictiemodel interne herniatie' / n for n in
        ('OPERATE_algoritme270724.sav', 'OPERATE_algoritme200924.sav')]

# Field-specific current CRF blocks manually located from the complete pages.
# Boundaries isolate this question's options, rather than testing a whole page.
REVIEWED_BLOCKS = {
    'Geslacht': (6, 'Geslacht', '1.2', 'unconditional baseline'),
    'episodes': (6, 'Is dit de eerste episode van pa-', 'If 1.3', 'unconditional episode flag'),
    'VG_cholecystectomie': (17, 'Cholecystectomie', '3.25', 'unconditional baseline history; history date still required'),
    'Roken': (45, 'Rookt patient', 'If 12.1', 'unconditional baseline smoking question'),
    'SEH_1stecontact': (59, 'Was het eerste studie contact via', 'If 15.1', 'unconditional first-setting flag; No alone does not prove POLI'),
    'SEH_A_buikpijn': (60, 'Buikpijn?', 'If 15.1', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_COA': (62, 'Continu of in aanvallen?', 'If 16.11', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_vervoerspijn': (63, 'Vervoerspijn:', 'Exported at', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_bewegingsdrang': (64, 'Bewegingsdrang', 'If 15.1', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_misselijk': (64, 'Misselijk', 'If 15.1', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_braken_direct': (64, 'Komt het eten direct terug:', 'If 15.1', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_passageklachten': (64, 'Passageklachten', 'If 15.1', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_maagzuur': (64, 'Maagzuur?', 'Exported at', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_ontlasting': (65, 'Ontlasting afwijkend', 'If 16.23', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_mictie': (65, 'Mictie afwijkend?', 'If 16.24', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'SEH_A_koorts': (66, 'Koorts gehad?', 'If 16.25', 'SEH_1stecontact == 1; no symptom-positive parent'),
    'FU_einde_diagnose': (730, 'Was er aan het einde van de', 'If 134.5', 'end follow-up question134.5; outcome, never predictor'),
}


def sav_dictionary(path):
    """Read documented SAV dictionary record types; stop before case data."""
    blob = path.read_bytes()
    if blob[:4] not in (b'$FL2', b'$FL3'):
        raise ValueError(f'Unsupported SAV header: {path.name}')
    endian = '<' if struct.unpack('<i', blob[64:68])[0] in (2, 3) else '>'
    pos = 176
    variables, value_labels, longnames, pending = [], {}, {}, []

    def take(n):
        nonlocal pos
        if n < 0 or pos + n > len(blob):
            raise ValueError('Invalid dictionary record length')
        result = blob[pos:pos+n]
        pos += n
        return result

    def ints(n=1):
        return struct.unpack(endian + 'i' * n, take(4*n))

    def decode(value):
        try:
            return value.decode('utf-8').rstrip('\x00 ')
        except UnicodeDecodeError:
            return value.decode('cp1252').rstrip('\x00 ')

    while True:
        record = ints()[0]
        if record == 2:
            kind, has_label, missing_n, _, _ = ints(5)
            name = decode(take(8))
            label = ''
            if has_label:
                n = ints()[0]
                label = decode(take((n+3)//4*4)[:n])
            missing_raw = take(abs(missing_n)*8) if missing_n else b''
            variables.append(dict(short_name=name, label=label,
                                  type='numeric' if kind == 0 else 'string',
                                  continuation=kind == -1,
                                  declared_missing_count=missing_n))
        elif record == 3:
            pending = []
            for _ in range(ints()[0]):
                value = struct.unpack(endian+'d', take(8))[0]
                n = take(1)[0]
                label = decode(take(((n+1+7)//8)*8-1)[:n])
                pending.append((str(int(value)) if value.is_integer() else str(value), label))
        elif record == 4:
            for index in ints(ints()[0]):
                value_labels[index] = dict(pending)
        elif record == 6:
            take(ints()[0]*80)
        elif record == 7:
            subtype, size, count = ints(3)
            payload = take(size*count)
            if subtype == 13:
                for pair in decode(payload).split('\t'):
                    if '=' in pair:
                        key, value = pair.split('=', 1)
                        longnames[key] = value
        elif record == 999:
            ints()  # dictionary terminator filler
            break
        else:
            raise ValueError(f'Unknown SAV dictionary record {record} at {pos}')
    result = {}
    for i, variable in enumerate(variables, 1):
        if variable.pop('continuation'):
            continue
        name = longnames.get(variable['short_name'], variable['short_name'])
        result[name] = dict(variable, codes=value_labels.get(i, {}))
    return result


def norm(value):
    value = value.replace('Y es', 'Yes')
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode()
    value = re.sub(r'-\s*\n\s*', '', value)
    return re.sub(r'[^a-z0-9]+', ' ', value.lower()).strip()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dictionaries = [sav_dictionary(path) for path in SAVS]
    pdf_pages = [page.extract_text() or '' for page in PdfReader(PDF).pages]
    normalized = [norm(page) for page in pdf_pages]
    current = pd.read_csv(SOURCE, sep=';', dtype=str, keep_default_na=False)
    rows, contradictions = [], []
    for name in current.columns:
        a, b = [dictionary.get(name, {}) for dictionary in dictionaries]
        metadata = b or a
        codes = metadata.get('codes', {})
        concordant = bool(a and b and a.get('codes') == b.get('codes')
                          and a.get('label') == b.get('label') and a.get('type') == b.get('type'))
        label = metadata.get('label', '')
        question = label.split('|', 1)[-1]
        needle = norm(question)
        page_matches = [i for i, text in enumerate(normalized)
                        if len(needle) >= 6 and needle in text]
        # Require a unique question context; common repeated questions are ambiguous.
        context = pdf_pages[page_matches[0]] if len(page_matches) == 1 else ''
        text_codes = {k: v for k, v in codes.items() if not k.startswith('-')}
        confirmed_options = [k for k, v in text_codes.items()
                             if len(norm(v)) >= 2 and norm(v) in norm(context)]
        counts = current[name].str.strip().value_counts() if codes else pd.Series(dtype='int64')
        observed = {str(k): int(v) for k, v in counts.items()
                    if re.fullmatch(r'-?\d+(?:\.\d+)?', str(k))} if codes else {}
        nonnumeric_n = int(sum(v for k, v in counts.items()
                              if k and not re.fullmatch(r'-?\d+(?:\.\d+)?', str(k)))) if codes else 0
        outside = sorted(set(observed) - set(codes)) if codes else []
        if nonnumeric_n:
            outside.append('[NONNUMERIC_VALUES_REDACTED]')
        state = 'provisional' if codes and concordant and not outside else 'unresolved'
        evidence = []
        if metadata:
            evidence.append('Historical SPSS variable dictionary only; no case values read')
        if concordant:
            evidence.append('Exact name, question label, type and code-label concordance across both historical SAVs')
        if context:
            evidence.append('Unique exact normalized current CRF question text match')
        if codes:
            evidence.append('Current CRF exposes option text but no numerical code mapping; current-release code identity unverified')
        if outside:
            contradictions.append(dict(source_name=name, unexpected_codes=outside))
        dependencies = re.findall(r'If .*?answer this question:', context, re.S) if context else []
        reviewed = REVIEWED_BLOCKS.get(name)
        specific_context, applicability, reviewed_page = '', '', None
        if reviewed:
            reviewed_page, start, end, applicability = reviewed
            page = pdf_pages[reviewed_page-1]
            start_pos = page.find(start)
            end_pos = page.find(end, start_pos+len(start))
            if start_pos < 0 or end_pos < 0:
                raise ValueError(f'CRF question block boundary missing: {name}')
            specific_context = page[start_pos:end_pos]
        all_options_match = bool(text_codes and specific_context) and all(
            norm(v) in norm(specific_context) for v in text_codes.values())
        checked_mapping = bool(concordant and codes and not outside and all_options_match)
        rows.append(dict(source_name=name, label=label, type=metadata.get('type', 'unresolved'),
                         codes_json=json.dumps(codes, ensure_ascii=False), verification_state=state,
                         evidence='; '.join(evidence), current_crf_context=context,
                         current_crf_pdf_pages_json=json.dumps([i+1 for i in page_matches]),
                         dependencies_json=json.dumps(dependencies, ensure_ascii=False),
                         dependencies_scope='page_context_only_not_field_specific' if dependencies else 'not_extracted',
                         historical_concordant=concordant,
                         historical_mapping_current_question_and_options_concordant=checked_mapping,
                         reviewed_question_context=specific_context,
                         reviewed_question_pdf_page=reviewed_page,
                         reviewed_applicability=applicability,
                         option_text_confirmed_codes_json=json.dumps(confirmed_options),
                         observed_codes_json=json.dumps(observed), unexpected_codes_json=json.dumps(outside),
                         nonnumeric_observation_count=nonnumeric_n,
                         observed_code_coverage_valid=bool(codes) and not outside))
    pd.DataFrame(rows).to_csv(OUT / 'field_codebook.csv', index=False)
    (OUT / 'historical_value_labels.json').write_text(json.dumps(
        {p.name: d for p, d in zip(SAVS, dictionaries)}, ensure_ascii=False, indent=2))
    (OUT / 'crf_pages.json').write_text(json.dumps(
        [{'pdf_page': i+1, 'text': text} for i, text in enumerate(pdf_pages)], ensure_ascii=False, indent=2))
    (OUT / 'reviewed_mappings.json').write_text(json.dumps(
        [r for r in rows if r['reviewed_question_context']], ensure_ascii=False, indent=2))
    summary = dict(source_shape=list(current.shape), fields=len(rows),
                   question_and_options_checked_mapping_count=sum(
                       r['historical_mapping_current_question_and_options_concordant'] for r in rows),
                   status_counts=pd.Series([r['verification_state'] for r in rows]).value_counts().to_dict(),
                   historical_labelled_counts={p.name: sum(bool(v['codes']) for v in d.values())
                                               for p, d in zip(SAVS, dictionaries)},
                   contradictions=contradictions,
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [SOURCE, PDF, *SAVS]},
                   limitations=['Printed CRF contains neither variable names nor numerical codes.',
                                'No current-release numeric category mapping is verified by these sources.',
                                'Automatic question text matches/dependencies are page-level evidence; reviewed_question_context isolates manually located question blocks.',
                                'Unlabelled observations are not exported to avoid patient information leakage.'])
    (OUT / 'validation_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k not in ('source_sha256', 'contradictions')}))


if __name__ == '__main__':
    main()
