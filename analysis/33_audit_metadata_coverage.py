"""Read-only audit of current-source and historical form-metadata coverage.

No participant values, cleaning pipeline, imputations, or fitted models are read
or regenerated. The codebook describes the current source headers and verified
2024-to-current crosswalk, so absence of a verified mapping is kept unresolved.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import json
import hashlib

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/first_visit_2026-10-27_v3"
OUT = ROOT / "outputs/statistician_handover_2026-11-02/review"
OUT.mkdir(parents=True, exist_ok=True)
VISITS = {"Basis gegevens", "SEH bezoek", "Polikliniek bezoek", "Vragenlijsten Basis", "Aanvullend onderzoek"}


def read(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write(name, records):
    with (OUT / name).open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


codebook_path = SOURCE / "codebook/field_codebook.csv"
metadata_path = SOURCE / "codebook/study_variable_list_2024.csv"
codebook = read(codebook_path)
metadata = read(metadata_path)
refs = {}
inputs = [codebook_path, metadata_path]
for setting in ["POLI", "SEH"]:
    path = SOURCE / f"variable_dictionary_source_{setting}.csv"
    inputs.append(path)
    refs[setting] = {name for row in read(path) for name in json.loads(row["source_columns"] or "[]")}
all_refs = set.union(*refs.values())

current_gaps = []
current_summary = []
for visit in sorted(VISITS):
    current = [row for row in codebook if row["visit"] == visit]
    missing = [row for row in current if row["source_name"] not in all_refs]
    current_summary.append(dict(visit=visit, current_mapped_source_fields=len(current),
                                referenced_in_release=len(current)-len(missing),
                                not_referenced_in_release=len(missing)))
    for row in missing:
        name = row["source_name"]
        colon = name.startswith("AO_overig_colon_")
        reason = ("No retained same-day colonoscopy values; source extraction materializes result columns only for same-day matches. Frozen primary cohort has zero same-day colonoscopy matches in both settings."
                  if colon else
                  "General other-investigation source field is not extracted into the corrected first-visit release. It may contain investigations from any episode time; baseline availability is unverified.")
        current_gaps.append(dict(source_name=name, metadata_variable_name=row["metadata_variable_name"],
            visit=visit, form=row["form"], question=row["question"], field_type=row["type"], unit=row["units"],
            verification_state=row["verification_state"], disposition="outside_corrected_release",
            reason=reason, available_at_presentation="not established", observed_n="not assessed",
            followup_action="Keep inventoried; source/timing QC required before potential use; do not infer or impute values"))

by_name = defaultdict(list)
for row in codebook:
    if row["metadata_variable_name"]:
        by_name[row["metadata_variable_name"]].append(row)
    if row["source_name"] != row["metadata_variable_name"]:
        by_name[row["source_name"]].append(row)
historical_inventory = []
for row in metadata:
    if row["Visit name"] not in VISITS:
        continue
    name = row["Variable name"]
    candidates = sorted({r["source_name"] for r in by_name[name]})
    represented = sorted(set(candidates) & all_refs)
    if name in all_refs and name not in represented:
        represented.append(name)
    disposition = ("represented_in_corrected_release" if represented else
                   "current_source_exists_outside_corrected_release" if candidates else
                   "historical_metadata_without_verified_current_source_mapping")
    historical_inventory.append(dict(metadata_variable_name=name, visit=row["Visit name"],
        form=row["Form name"], question=row["Field label"],
        field_type=row["Exported field type"], disposition=disposition,
        current_source_names=json.dumps(candidates, ensure_ascii=False),
        represented_source_names=json.dumps(represented, ensure_ascii=False),
        note="Historical metadata is not proof that a variable was measured at presentation or that a changed name maps to a current source field."))

write("current_source_fields_outside_release.csv", current_gaps)
write("historical_firstvisit_metadata_coverage.csv", historical_inventory)
summary = dict(current_source_by_visit=current_summary,
    current_source_gaps=len(current_gaps),
    current_source_gap_forms=dict(Counter(r["form"] for r in current_gaps)),
    historical_metadata_fields=len(historical_inventory),
    historical_dispositions=dict(Counter(r["disposition"] for r in historical_inventory)),
    scope="Complete inventory of corrected v3 physical fields plus separate current-source/historical-metadata gaps; not a claim of exhaustive available baseline measurements in the original export.",
    limitations=["No participant raw-source values read; availability of the 40 general AO fields is unassessed.",
        "Diagnostic forms can refer to later workup, so existence in the raw export does not establish baseline timing.",
        "Unmapped historical variables are not assumed absent or equivalent to similarly named current variables."],
    input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
(OUT / "metadata_coverage_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False))
