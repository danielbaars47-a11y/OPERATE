#!/usr/bin/env python3
"""Audit Castor record links, without assigning patient IDs or altering raw data.

Run with python3 OPERATE/analysis/04_audit_operate_linkage.py. Uses standard
library only; outputs an exact-link ledger and aggregate summary in
OPERATE/results/cleaned_2026-09-30. Repeated targets are not deduplicated.
"""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'OPERATE data' / 'Buikpijn_na_bariatrische_chirurg_export_20260622.csv'
OUTPUT = ROOT / 'results' / 'cleaned_2026-09-30'
FIELDS = ['Episode_nummer_1ste episode_Nummer', 'Episode_nummer_2de episode_Nummer',
          'Episode_nummer_3de episode_Nummer', 'Episode_nummer_4de episode_Nummer']
SENTINELS = {'-95', '-96', '-97', '-98', '-99'}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identifier_pattern(value):
    # Syntax only: no clinical interpretation of the prefixes.
    if value.isdigit():
        return 'numeric'
    if value.casefold().startswith('cc'):
        return 'cc_prefix_case_insensitive'
    if value.casefold().startswith('c'):
        return 'c_prefix_case_insensitive'
    return 'other'


def comparable_equal(left, right):
    if not left or not right or left in SENTINELS or right in SENTINELS:
        return 'not_comparable'
    return 'yes' if left == right else 'no'


def first_visit(row):
    seh, poli = row['SEH_1stecontact'] == '1', row['POLI_1stecontact'] == '1'
    if seh and not poli:
        return 'SEH', row['SEH_datum']
    if poli and not seh:
        return 'POLI', row['POLI_datum']
    return ('both' if seh and poli else 'unassigned'), ''


def audit_candidate_duplicates(records, columns):
    """Flag shared complete demographic triplets; do not infer identity.

    Equality of sex, birth date and bariatric surgery date is a candidate key,
    not a validated person identifier. Administrative variation or different
    visit dates do not establish whether a pair represents one or two people.
    """
    key_fields = ['Geslacht', 'Geboortedatum', 'B_operatie_datum']
    groups = defaultdict(list)
    for row in records:
        if row['episodes'] == '1':
            key = tuple(row[field] for field in key_fields)
            if all(value and value not in SENTINELS and not value.endswith('2999')
                   and not value.endswith('2995') and not value.endswith('2996')
                   and not value.endswith('2997') and not value.endswith('2998') for value in key):
                groups[key].append(row)
    pairs = []
    administrative = {'Participant Id', 'Participant Creation Date', 'Participant Status', 'Site Abbreviation'}
    for group_number, group in enumerate((v for v in groups.values() if len(v) > 1), 1):
        for left, right in combinations(group, 2):
            different = [field for field in columns if left[field] != right[field]]
            both_values_different = [field for field in different if left[field] and right[field]]
            clinical_different = [field for field in different if field not in administrative]
            shared_long_text = [field for field in columns if len(left[field]) >= 30 and left[field] == right[field]]
            left_setting, left_date = first_visit(left)
            right_setting, right_date = first_visit(right)
            direct_link = any(left[field] == right['Participant Id'] or right[field] == left['Participant Id'] for field in FIELDS)
            left_incoming = sum(row[field] == left['Participant Id'] for row in records for field in FIELDS)
            right_incoming = sum(row[field] == right['Participant Id'] for row in records for field in FIELDS)
            exact_duplicate = not clinical_different
            pairs.append({
                'candidate_group_id': f'first_episode_triplet_{group_number:03d}',
                'candidate_group_records': len(group),
                'participant_id_a': left['Participant Id'],
                'participant_id_b': right['Participant Id'],
                'setting_a': left_setting, 'setting_b': right_setting,
                'index_date_a_raw': left_date, 'index_date_b_raw': right_date,
                'same_observed_index_date': int(bool(left_date) and left_date == right_date),
                'complete_matching_key_fields': '|'.join(key_fields),
                'differing_source_column_count': len(different),
                'differing_nonadministrative_column_count': len(clinical_different),
                'differing_both_populated_column_count': len(both_values_different),
                'identical_long_text_field_count': len(shared_long_text),
                'identical_long_text_field_names': '|'.join(shared_long_text),
                'direct_episode_link_between_pair': int(direct_link),
                'incoming_episode_links_a': left_incoming,
                'incoming_episode_links_b': right_incoming,
                'exact_duplicate_nonadministrative_record': int(exact_duplicate),
                'classification': 'confirmed_duplicate_record' if exact_duplicate else 'possible_same_patient_requires_verification',
                'identity_confirmed': int(exact_duplicate),
                'relation_evidence': 'Complete demographic triplet matches; ' + ('explicit episode link also present; ' if direct_link else 'no explicit episode link between pair; ') + ('all nonadministrative source values identical.' if exact_duplicate else 'nonadministrative values differ; demographic match alone does not confirm identity.'),
                'recommended_action': 'Review source identity before selecting one observation per patient; no automated record removal.',
                'differing_source_column_names': '|'.join(different),
            })
    fields = list(pairs[0]) if pairs else ['candidate_group_id', 'participant_id_a', 'participant_id_b', 'classification']
    with (OUTPUT / 'candidate_duplicate_first_episodes.csv').open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(pairs)
    return {
        'screening_population': 'episodes code1; complete non-sentinel sex, birth date, surgery date triplet',
        'candidate_key_fields': key_fields,
        'candidate_groups': len({row['candidate_group_id'] for row in pairs}),
        'candidate_pairs': len(pairs),
        'candidate_records': len({row[field] for row in pairs for field in ['participant_id_a', 'participant_id_b']}),
        'confirmed_duplicate_records_pairs': sum(row['identity_confirmed'] for row in pairs),
        'possible_same_patient_pairs': sum(not row['identity_confirmed'] for row in pairs),
        'records_removed': 0,
        'identity_rule': 'Matching demographic triplets identify candidates only, not confirmed persons. No record merge or deletion.',
        'audit_file': 'candidate_duplicate_first_episodes.csv',
    }


def main():
    source_hash = sha256(SOURCE)
    with SOURCE.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter=';')
        columns = reader.fieldnames
        records = list(reader)
    required = ['Participant Id', 'episodes', 'Episode_hoeveel', 'Geslacht', 'Geboortedatum', *FIELDS]
    assert all(field in columns for field in required), 'Required column missing'
    ids = [row['Participant Id'] for row in records]
    assert len(ids) == len(set(ids)) and all(ids), 'Missing or duplicated record key'
    by_id = {row['Participant Id']: row for row in records}
    pair_counts = Counter((field, row[field]) for row in records for field in FIELDS if row[field])
    target_counts = Counter(row[field] for row in records for field in FIELDS if row[field])
    audit = []
    for source_number, source in enumerate(records, start=1):
        for target_episode_number, field in enumerate(FIELDS, start=1):
            raw_target = source[field]
            if not raw_target:
                continue
            target = by_id.get(raw_target)
            self_link = raw_target == source['Participant Id']
            expected_flag = '1' if target_episode_number == 1 else '0'
            # CRF: Episode_hoeveel is 1=second, 2=third, etc., not a count.
            expected_order = '' if target_episode_number == 1 else str(target_episode_number - 1)
            role_consistent = target is not None and target['episodes'] == expected_flag and (
                target_episode_number == 1 or target['Episode_hoeveel'] == expected_order)
            status = ('self_link' if self_link else 'unresolved_target' if target is None else
                      'target_episode_role_inconsistent' if not role_consistent else 'matched_role_consistent')
            audit.append({
                'source_record_number': source_number,
                'source_participant_id': source['Participant Id'],
                'source_episodes_code': source['episodes'],
                'source_episode_hoeveel_code': source['Episode_hoeveel'],
                'link_field': field,
                'target_episode_number': target_episode_number,
                'target_participant_id_raw': raw_target,
                'target_value_has_outer_whitespace': int(raw_target != raw_target.strip()),
                'exact_target_match': int(target is not None),
                'trimmed_target_match_only': int(target is None and raw_target.strip() in by_id),
                'self_link': int(self_link),
                'target_episodes_code': target['episodes'] if target else '',
                'target_episode_hoeveel_code': target['Episode_hoeveel'] if target else '',
                'target_role_matches_requested_episode': int(role_consistent),
                'source_is_repeat_episode': int(source['episodes'] == '0'),
                'target_reference_count_same_link_field': pair_counts[(field, raw_target)],
                'target_reference_count_all_link_fields': target_counts[raw_target],
                'sex_codes_equal': comparable_equal(source['Geslacht'], target['Geslacht']) if target else 'not_comparable',
                'birth_dates_equal': comparable_equal(source['Geboortedatum'], target['Geboortedatum']) if target else 'not_comparable',
                'audit_status': status,
            })
    fields_summary = {}
    for field in FIELDS:
        subset = [row for row in audit if row['link_field'] == field]
        referenced = Counter(row['target_participant_id_raw'] for row in subset)
        n_match = sum(row['exact_target_match'] for row in subset)
        fields_summary[field] = {
            'populated_links': len(subset), 'exact_matches': n_match,
            'exact_match_percent': round(100 * n_match / len(subset), 3) if subset else None,
            'unresolved_targets': len(subset) - n_match,
            'self_links': sum(row['self_link'] for row in subset),
            'unique_referenced_targets': len(referenced),
            'targets_referenced_more_than_once': sum(n > 1 for n in referenced.values()),
            'excess_references_to_repeated_targets': sum(n - 1 for n in referenced.values()),
            'matched_target_episodes_code_counts': dict(Counter(row['target_episodes_code'] for row in subset if row['exact_target_match'])),
            'audit_status_counts': dict(Counter(row['audit_status'] for row in subset)),
            'sex_code_disagreements': sum(row['sex_codes_equal'] == 'no' for row in subset),
            'birth_date_disagreements': sum(row['birth_dates_equal'] == 'no' for row in subset),
        }
    patterns = sorted(set(identifier_pattern(value) for value in ids))
    summary = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'source_relative_path': str(SOURCE.relative_to(ROOT)), 'source_sha256': source_hash,
        'source_rows': len(records), 'source_columns': len(columns),
        'record_identifier': 'Participant Id', 'record_ids_unique_nonmissing': True,
        'record_ids_have_outer_whitespace': sum(value != value.strip() for value in ids),
        'link_match_rule': 'Exact string match to Participant Id; no correction, patient ID assignment, merge or deduplication.',
        'source_episodes_code_counts': dict(Counter(row['episodes'] for row in records)),
        'source_episode_hoeveel_code_counts': dict(Counter(row['Episode_hoeveel'] for row in records)),
        'identifier_pattern_groups': {pattern: {
            'records': sum(identifier_pattern(row['Participant Id']) == pattern for row in records),
            'episodes_code_counts': dict(Counter(row['episodes'] for row in records if identifier_pattern(row['Participant Id']) == pattern)),
        } for pattern in patterns},
        'link_fields': fields_summary,
        'total_populated_links': len(audit),
        'total_exact_matches': sum(row['exact_target_match'] for row in audit),
        'total_unresolved_links': sum(not row['exact_target_match'] for row in audit),
        'total_self_links': sum(row['self_link'] for row in audit),
        'unique_source_records_with_populated_links': len({row['source_participant_id'] for row in audit}),
        'repeat_episode_records_without_first_episode_link': sum(row['episodes'] == '0' and not row[FIELDS[0]] for row in records),
        'review_notes': [
            'Participant Id identifies an export record, not necessarily a unique person.',
            'Repeated targets may legitimately reflect several episodes of one person; no repeated link was deleted.',
            'ID prefix groups describe syntax only; no clinical group is assigned from a prefix.',
            'CRF PDF page6: first-episode question Yes/No; Episode_hoeveel categories are second/third/fourth/fifth episode. Numeric role codes follow the existing source conventions and require official Castor dictionary confirmation.',
            'Role-inconsistent, unresolved, self-linked or demographic-discordant records require source verification before longitudinal analysis.',
            'The ledger preserves identifiers locally; this summary contains no individual identifiers.',
        ],
        'source_modified': False,
    }
    assert sha256(SOURCE) == source_hash, 'Source changed during audit'
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary['candidate_duplicate_first_episodes'] = audit_candidate_duplicates(records, columns)
    with (OUTPUT / 'record_linkage_audit.csv').open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit[0]))
        writer.writeheader()
        writer.writerows(audit)
    (OUTPUT / 'linkage_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({key: summary[key] for key in [
        'source_rows', 'source_columns', 'total_populated_links', 'total_exact_matches',
        'total_unresolved_links', 'total_self_links', 'repeat_episode_records_without_first_episode_link'
    ]}, indent=2))
    print('Wrote record_linkage_audit.csv, candidate_duplicate_first_episodes.csv and linkage_summary.json')
    print(json.dumps(summary['candidate_duplicate_first_episodes'], indent=2))


if __name__ == '__main__':
    main()
