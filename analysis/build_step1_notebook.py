#!/usr/bin/env python3
"""Build and execute the reader-facing step-1 audit notebook.

The bundled environment does not include nbformat/nbclient, so this small builder
uses the public .ipynb JSON schema and executes each Python cell in one shared
namespace before writing the notebook with captured outputs.
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import traceback
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "01_gallstone_step1_data_quality.ipynb"


def markdown(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(keepends=True)}


def code(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


cells = [
    markdown(
        """## tl;dr

Stap 1 bevat 789 eerste SEH-presentaties en 776 eerste poliklinische presentaties. Volgens de aangeleverde tabel met definities zijn er 96 bevestigde SEH-events en 95 bevestigde poli-events. Dertien dossiers blijven outcome-missing tot klinische adjudicatie (11 SEH, 2 poli). Na uitsluiting van negatieve controles met een eerdere cholecystectomie en van records met onzekere eligibility zijn in beide cohorten 578 records primair analyseerbaar. Er is nog geen model gefit.
"""
    ),
    markdown(
        """## Context & Methods

Doel: een reproduceerbare eerste stap voor een diagnostisch predictiemodel voor definitief galsteen-gerelateerd lijden na bariatrische chirurgie.

### Key Assumptions

- Alleen `episodes = 1` wordt gebruikt.
- SEH en polikliniek blijven afzonderlijke cohorten.
- Alle drie definitieve diagnoses aan het einde van follow-up dragen bij aan de outcome.
- Cholecystolithiasis vereist een passend definitief diagnoselabel plus laparoscopische cholecystectomie, laparoscopische cholecystolithiotomie of steen-positieve ERCP op of na het indexmoment.
- Cholecystitis vereist een passend definitief diagnoselabel plus positieve echo, CT of relevante chirurgische bevestiging op of na het indexmoment.
- Een procedure zonder passend definitief diagnoselabel is geen event.
- Een galsteenlabel zonder het vereiste tabelbewijs en andere mogelijke galsteen-gerelateerde diagnoses, zoals biliaire pancreatitis of cholangitis, blijven te adjudiceren en worden niet als controles gebruikt.
- DEF-codes 5 en 6 identificeren de kandidaatlabels; de oude eerste-diagnosecodes 1 en 3 worden alleen gebruikt als de nieuwe DEF-code ontbreekt.
- Negatieve controles met een eerdere cholecystectomie worden uitgesloten. Positieve cases met eerdere cholecystectomie blijven behouden.
- Onderzoek wordt alleen als indexinformatie gekoppeld bij dezelfde datum als de eerste presentatie.
- Outcome, structureel niet verricht onderzoek en onbekende controle-eligibility worden niet geïmputeerd.
- Omdat de aparte Castor-variabelenlijst ontbreekt, blijven niet-ondubbelzinnige numerieke categorische codes ongewijzigd; hun frequenties staan in de categorie-inventaris voor latere labelkoppeling.
"""
    ),
    markdown("## Data\n\nDe bron is de lokale Castor CSV-export. `Presumed and definitive table picture diagnosis.jpeg` bepaalt de definitieve-diagnosecriteria; het lege CRF is gebruikt om de numerieke episode-, setting-, diagnose-, operatie- en antwoordcodes te vertalen."),
    code(
        """from pathlib import Path
import os
import json
import subprocess
import sys
import pandas as pd

ROOT = Path.cwd()
PIPELINE = ROOT / "analysis" / "01_prepare_gallstone_step1.py"
OUTPUT = ROOT / "results" / "step1_gallstone"
env = dict(os.environ)
env["PYTHONWARNINGS"] = "ignore"
run = subprocess.run([sys.executable, str(PIPELINE)], cwd=ROOT, env=env, capture_output=True, text=True, check=True)
print("Pipeline executed successfully")
print("Output directory:", OUTPUT)
"""
    ),
    markdown("## Results\n\n### Cohort flow"),
    code(
        """flow = pd.read_csv(OUTPUT / "cohort_flow.csv")
print(flow.pivot(index="stage", columns="setting", values="n").fillna(0).astype(int).to_string())
"""
    ),
    markdown("### Outcome en indexonderzoek"),
    code(
        """rows = []
for setting in ["SEH", "POLI"]:
    cohort = pd.read_csv(OUTPUT / f"cohort_{setting}_first_presentations.csv", low_memory=False)
    rows.append({
        "cohort": setting,
        "first_presentations": len(cohort),
        "events": int(cohort["gallstone_related_disease"].eq(1).sum()),
        "adjudication_required": int(cohort["outcome_status"].str.startswith("adjudication_required", na=False).sum()),
        "eligible_controls": int(cohort["gallstone_related_disease"].eq(0).sum()),
        "primary_analysis_ready": int(cohort["primary_analysis_include"].eq(1).sum()),
        "same_day_lab": int(pd.to_numeric(cohort["index_lab_same_day_record_count"], errors="coerce").fillna(0).gt(0).sum()),
        "same_day_CT": int(pd.to_numeric(cohort["index_ct_same_day_record_count"], errors="coerce").fillna(0).gt(0).sum()),
        "same_day_echo": int(pd.to_numeric(cohort["index_echo_same_day_record_count"], errors="coerce").fillna(0).gt(0).sum()),
    })
print(pd.DataFrame(rows).to_string(index=False))
"""
    ),
    markdown("### Toepassing van de definitieve-diagnosetabel"),
    code(
        """definition_rows = []
for setting in ["SEH", "POLI"]:
    cohort = pd.read_csv(OUTPUT / f"cohort_{setting}_first_presentations.csv", low_memory=False)
    definition_rows.append({
        "cohort": setting,
        "confirmed_events": int(cohort["gallstone_related_disease"].eq(1).sum()),
        "cholecystolithiasis_definition_met": int(cohort["outcome_definition_met_cholecystolithiasis"].eq(1).sum()),
        "cholecystitis_definition_met": int(cohort["outcome_definition_met_cholecystitis"].eq(1).sum()),
        "both_definitions_met": int((cohort["outcome_definition_met_cholecystolithiasis"].eq(1) & cohort["outcome_definition_met_cholecystitis"].eq(1)).sum()),
        "label_without_required_evidence": int(cohort["outcome_status"].eq("adjudication_required_label_without_table_evidence").sum()),
        "other_possible_gallstone_diagnosis": int(cohort["outcome_status"].eq("adjudication_required_other_gallstone_related_diagnosis").sum()),
        "changed_from_previous_provisional_rule": int(cohort["outcome_changed_from_previous_provisional"].eq(1).sum()),
    })
print(pd.DataFrame(definition_rows).to_string(index=False))
"""
    ),
    markdown("### Reproduceerbare validatie van de eventcodering"),
    code(
        """seh = pd.read_csv(OUTPUT / "cohort_SEH_first_presentations.csv", low_memory=False)
poli = pd.read_csv(OUTPUT / "cohort_POLI_first_presentations.csv", low_memory=False)
combined = pd.concat([seh.assign(_setting="SEH"), poli.assign(_setting="POLI")], ignore_index=True)
flow = pd.read_csv(OUTPUT / "cohort_flow.csv")
adjudication = pd.read_csv(OUTPUT / "outcome_adjudication_review.csv")
changes = pd.read_csv(OUTPUT / "outcome_definition_change_review.csv")
metadata = json.loads((OUTPUT / "analysis_metadata.json").read_text(encoding="utf-8"))

checks = {}
checks["one_record_per_patient_SEH"] = seh["Participant Id"].is_unique
checks["one_record_per_patient_POLI"] = poli["Participant Id"].is_unique
checks["no_setting_overlap"] = set(seh["Participant Id"]).isdisjoint(set(poli["Participant Id"]))
checks["reference_image_recorded"] = metadata["outcome_definition"]["reference_image"] == "Presumed and definitive table picture diagnosis.jpeg"
checks["reference_image_exists"] = (ROOT / metadata["outcome_definition"]["reference_image"]).exists()

for setting, cohort, expected_events, expected_adjudication, expected_ready in [
    ("SEH", seh, 96, 11, 578),
    ("POLI", poli, 95, 2, 578),
]:
    prefix = setting.lower()
    positive = cohort["gallstone_related_disease"].eq(1)
    negative = cohort["gallstone_related_disease"].eq(0)
    adjudicate = cohort["outcome_status"].str.startswith("adjudication_required", na=False)
    lith_met = cohort["outcome_definition_met_cholecystolithiasis"].eq(1)
    itis_met = cohort["outcome_definition_met_cholecystitis"].eq(1)
    lith_evidence = cohort[[
        "outcome_evidence_laparoscopic_cholecystectomy",
        "outcome_evidence_laparoscopic_cholecystolithotomy",
        "outcome_evidence_positive_ercp",
    ]].eq(1).any(axis=1)
    itis_evidence = cohort[[
        "outcome_evidence_ultrasound_cholecystitis",
        "outcome_evidence_ct_cholecystitis",
        "outcome_evidence_surgery_cholecystitis",
    ]].eq(1).any(axis=1)
    checks[f"{prefix}_event_count"] = int(positive.sum()) == expected_events
    checks[f"{prefix}_adjudication_count"] = int(adjudicate.sum()) == expected_adjudication
    checks[f"{prefix}_analysis_ready_count"] = int(cohort["primary_analysis_include"].eq(1).sum()) == expected_ready
    checks[f"{prefix}_positive_iff_table_definition_met"] = positive.equals(lith_met | itis_met)
    checks[f"{prefix}_lithiasis_requires_label"] = bool((~lith_met | cohort["outcome_candidate_cholecystolithiasis"].eq(1)).all())
    checks[f"{prefix}_lithiasis_requires_table_evidence"] = bool((~lith_met | lith_evidence).all())
    checks[f"{prefix}_cholecystitis_requires_label"] = bool((~itis_met | cohort["outcome_candidate_cholecystitis"].eq(1)).all())
    checks[f"{prefix}_cholecystitis_requires_table_evidence"] = bool((~itis_met | itis_evidence).all())
    checks[f"{prefix}_procedure_alone_not_event"] = bool((~positive | cohort[["outcome_candidate_cholecystolithiasis", "outcome_candidate_cholecystitis"]].eq(1).any(axis=1)).all())
    checks[f"{prefix}_adjudication_has_no_binary_outcome"] = bool(cohort.loc[adjudicate, "gallstone_related_disease"].isna().all())
    checks[f"{prefix}_adjudication_not_analysis_ready"] = bool(cohort.loc[adjudicate, "primary_analysis_include"].eq(0).all())
    checks[f"{prefix}_eligible_controls_no_prior_cholecystectomy"] = bool(cohort.loc[negative, "prior_cholecystectomy_status"].eq("no_prior_cholecystectomy").all())
    checks[f"{prefix}_status_partitions_all_rows"] = int(positive.sum() + negative.sum() + adjudicate.sum() + cohort["outcome_status"].eq("negative_excluded_prior_cholecystectomy").sum() + cohort["outcome_status"].eq("negative_control_eligibility_unknown").sum() + cohort["outcome_status"].eq("outcome_missing").sum()) == len(cohort)

checks["adjudication_review_reconciles"] = len(adjudication) == int(combined["outcome_status"].str.startswith("adjudication_required", na=False).sum())
checks["definition_change_review_reconciles"] = len(changes) == int(combined["outcome_changed_from_previous_provisional"].eq(1).sum()) == 11
checks["flow_event_totals_reconcile"] = all(
    int(flow.loc[(flow["setting"] == setting) & (flow["stage"] == "gallstone_positive_total"), "n"].iloc[0])
    == int(cohort["gallstone_related_disease"].eq(1).sum())
    for setting, cohort in [("SEH", seh), ("POLI", poli)]
)

failed = [name for name, passed in checks.items() if not bool(passed)]
assert not failed, f"Failed validation checks: {failed}"
print(f"All {len(checks)} event-definition validation checks passed.")
"""
    ),
    markdown("### Controle op waardebereik"),
    code(
        """flags = pd.read_csv(OUTPUT / "value_quality_flags.csv")
print("Gemarkeerde waarden:", len(flags))
if len(flags):
    print(flags.groupby(["setting", "variable"]).size().sort_values(ascending=False).head(20).to_string())
"""
    ),
    markdown(
        """## Takeaways

- De tabeldefinitie bevestigt 96 SEH-events en 95 poli-events. De som van de twee subtypen kan hoger zijn dan het totaal doordat twee SEH-patiënten aan beide definities voldoen.
- Elf SEH-dossiers en twee poli-dossiers blijven buiten event- én controlegroep tot klinische adjudicatie.
- Same-day laboratorium is vrijwel standaard op de SEH (743/789), maar zeldzaam op de poli (26/776); een primair polimodel met routinematige labwaarden is daarom niet passend.
- Beeldvorming op de indexdatum is veel beperkter beschikbaar, vooral op de poli, en hoort uitsluitend in een secundair extended model.
- De kalenderdatum geeft geen binnen-dagvolgorde. Ook bestaan enkele conflicten tussen onderzoeksdatum en de CRF-vlag; gebruik daarom de afzonderlijke linkage-status bij de definitieve predictordefinitie.
- Predictorselectie, transformaties, imputatie en modellering zijn bewust nog niet uitgevoerd.
- De volgende stap is adjudicatie van de dertien open dossiers en daarna parameterbudget/kandidaatpredictorselectie met de statisticus.
"""
    ),
]


namespace: dict[str, object] = {}
execution_count = 0
for cell in cells:
    if cell["cell_type"] != "code":
        continue
    execution_count += 1
    stream = io.StringIO()
    try:
        with contextlib.redirect_stdout(stream):
            exec(compile("".join(cell["source"]), f"cell-{execution_count}", "exec"), namespace)
    except Exception as exc:
        cell["outputs"] = [
            {
                "output_type": "error",
                "ename": type(exc).__name__,
                "evalue": str(exc),
                "traceback": traceback.format_exc().splitlines(),
            }
        ]
        raise
    cell["execution_count"] = execution_count
    output = stream.getvalue()
    if output:
        cell["outputs"] = [{"output_type": "stream", "name": "stdout", "text": output.splitlines(keepends=True)}]

notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": sys.version.split()[0]},
        "codex_validation": {
            "status": "executed_top_to_bottom",
            "engine": "stdlib_shared_namespace_fallback",
            "reason": "nbformat and nbclient were not installed in the bundled runtime",
        },
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

NOTEBOOK.parent.mkdir(parents=True, exist_ok=True)
NOTEBOOK.write_text(json.dumps(notebook, indent=1, ensure_ascii=False), encoding="utf-8")

with NOTEBOOK.open(encoding="utf-8") as handle:
    validated = json.load(handle)
assert validated["nbformat"] == 4
assert all(cell.get("execution_count") is not None for cell in validated["cells"] if cell["cell_type"] == "code")
assert not any(
    output.get("output_type") == "error"
    for cell in validated["cells"]
    for output in cell.get("outputs", [])
)
print(f"Built and executed {NOTEBOOK}")
