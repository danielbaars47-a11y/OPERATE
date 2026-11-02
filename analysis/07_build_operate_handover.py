#!/usr/bin/env python3
"""Build and execute a compact notebook and HTML handover for OPERATE cleaning."""
from pathlib import Path
import contextlib
import html
import io
import json
import sys
import nbformat
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/first_visit_2026-09-30"
MASTER = ROOT / "results/cleaned_2026-09-30"
metadata = json.loads((OUT / "validation_summary.json").read_text())
ascertainment = json.loads((OUT / "outcome_ascertainment_summary.json").read_text())
linkage = json.loads((MASTER / "linkage_summary.json").read_text())
rows = []
for setting, m in metadata["cohorts"].items():
    rows.append({"Setting": setting, "First presentations": m["n_first_visits"],
                 "Provisionally positive": m["n_existing_positive_events"],
                 "Provisional controls": m["n_existing_primary_nonevents"],
                 "Clean first-visit columns": m["n_analysis_columns"]})
table = pd.DataFrame(rows)
scope = """These are the complete first-visit cleaning outputs requested in the staged protocol.
No clinical predictor selection, category pooling, imputation or model fitting has been performed.
Outcomes are stored separately and retain the previous project definition. All outcome coding remains
provisional until the negative-outcome definition, ascertainment and flagged cases are resolved.
The source export and earlier analysis files have not been overwritten."""
limitations = [
"First episodes are selected from the explicit CRF flag. Record IDs are unique, but four possible duplicate pairs (eight records) need source identity review before a one-person analysis freeze. None was automatically deleted.",
"Twenty-eight first episodes have no resolved first-presentation setting. They remain in a separate file and in the full master.",
"Of the current controls, 318 SEH and 410 POLI records say no definitive diagnosis was established. Agree the operational negative-outcome definition and follow-up adequacy; this is not an instruction to automatically exclude all such records.",
"Thirteen previously identified gallstone outcome cases require clinical adjudication (11 SEH and two POLI). Prior-cholecystectomy negatives and unresolved cholecystectomy history remain excluded by the existing rule.",
"Numeric category labels and conditional checkbox semantics still need the official machine-readable Castor data dictionary. Unresolved codes are preserved, not given guessed labels or pooled.",
"Histories require a date before index; incomplete 30-day postoperative windows, undated context, and known retrospective imaging fields are withheld from the analysis values and retained in the source snapshot.",
"Same-day diagnostic linkage does not establish availability before the clinician's assessment. Linkage/performance conflicts remain explicit for later review.",
"Numerical plausibility flags do not prove transcription error. Source values are retained; screened numeric companions are missing until source review. Unspecified or ambiguous laboratory units remain marked as such."
]
files = [
("SEH cleaned first-visit data", "cohort_SEH_cleaned.csv.gz"),
("POLI cleaned first-visit data", "cohort_POLI_cleaned.csv.gz"),
("SEH missingness for every variable", "missingness_SEH_all_first_visits.csv"),
("POLI missingness for every variable", "missingness_POLI_all_first_visits.csv"),
("SEH dictionary and provenance", "variable_dictionary_SEH.csv"),
("POLI dictionary and provenance", "variable_dictionary_POLI.csv"),
("SEH provisional outcomes", "outcomes_SEH_existing_definition.csv"),
("POLI provisional outcomes", "outcomes_POLI_existing_definition.csv"),
("Outcome ascertainment review", "outcome_ascertainment_review.csv"),
("Numerical/source review flags", "quality_flags.csv"),
("Timing restriction summary", "timing_summary.csv"),
("Full cleaned master (all source records and variables)", "../cleaned_2026-09-30/operate_cleaned_master.csv.gz"),
("Possible duplicate first episodes", "../cleaned_2026-09-30/candidate_duplicate_first_episodes.csv"),
("Executed audit notebook", "../../notebooks/02_operate_first_visit_cleaning.ipynb")
]
readme = """OPERATE first-visit data package - 30 September 2026

""" + scope + "\n\n" + table.to_string(index=False) + """

START HERE
Use cohort_SEH_cleaned.csv.gz and cohort_POLI_cleaned.csv.gz for first-visit values.
Use outcomes_{setting}_existing_definition.csv for provisional outcome labels.
Join on participant_id only after checking uniqueness; both files preserve the same
record order. Do not use outcome/follow-up fields as predictors.

The cleaned files retain every first-visit variable, with timing-restricted cells
masked and new numeric companions appended. No clinical candidate-variable selection
has occurred. Matching analysis_missing_status files explain cell states; matching
source_missing_status files distinguish original absence from later timing masking.
before_timing_filter snapshots preserve source values for review.
*_numeric columns contain exact, screened numeric values; *_bound and *_qualifier
retain measurement/censoring information. Never substitute a detection-limit bound
for an exact observation without an analysis-specific decision.
Source measurement strings remain alongside numeric companions. Value-labelled
categorical coding remains pending the official Castor codebook.

READING THE DATA
Python: pd.read_csv(path, dtype=str, keep_default_na=False)
The empty string is an absent analytical value. Literal NA is preserved as source
text and must not be silently reinterpreted by default CSV NA parsing.
Convert numeric companion fields using the dictionary; use the status matrix to
distinguish observed data from structural or unresolved absence.
R: read.csv(gzfile(path), colClasses="character", na.strings="",
            check.names=FALSE, stringsAsFactors=FALSE)
Do not fit a model against every column automatically.

REVIEW BEFORE IMPUTATION
""" + "\n".join("- " + t for t in limitations) + """

REPRODUCING
The scripts use project-relative paths and never overwrite source data.
From the PHD workspace, use the same Python environment for:
  python OPERATE/analysis/02_clean_operate_master.py
  python OPERATE/analysis/04_audit_operate_linkage.py
  python OPERATE/analysis/03_prepare_operate_first_visit.py
  python OPERATE/analysis/05_audit_outcome_ascertainment.py
  python OPERATE/analysis/06_validate_operate_master.py
  python OPERATE/analysis/07_build_operate_handover.py
The last command requires nbformat. Core cleaning requires pandas and numpy.
The isolated environment at OPERATE/analysis/.venv supplies the notebook validator.
Exact runtime versions are recorded in validation metadata.

IMPUTATION HANDOVER
Agree outcome/control definitions and resolve identity/coding/value issues, then
choose predictors and parameters with the statistician. Genuine sporadic missingness
may be multiply imputed separately within SEH and POLI, conditioning on the relevant
predictors and outcome. Never impute outcome/setting or structurally unperformed
investigations. Integrate imputation with internal resampling for prediction
validation. No imputed output is delivered at this stage, as requested by the
detailed staged protocol.

Method references:
https://amices.org/mice/reference/mice.html
https://amices.org/mice/reference/pool.html
"""
(OUT / "README.txt").write_text(readme, encoding="utf-8")
md = nbformat.v4.new_markdown_cell
cc = nbformat.v4.new_code_cell
cells = [
md("# OPERATE: first-visit cleaning audit\n\n" + scope),
md("## Context and methods\n\n" + "\n\n".join(limitations[:3]) + "\n\nSee analysis/first_visit_protocol_2026-09-30.txt and the scripts for the complete workflow. This notebook validates saved outputs; rerun the listed scripts to regenerate them."),
cc("""from pathlib import Path
import hashlib
import json
import pandas as pd
candidates = [Path.cwd(), *Path.cwd().parents]
ROOT = next((p if (p / "OPERATE data").is_dir() else p / "OPERATE")
            for p in candidates if (p / "OPERATE data").is_dir() or (p / "OPERATE/OPERATE data").is_dir())
OUT = ROOT / "results/first_visit_2026-09-30"
MASTER = ROOT / "results/cleaned_2026-09-30"
meta = json.loads((OUT / "validation_summary.json").read_text())
source = ROOT / "OPERATE data/Buikpijn_na_bariatrische_chirurg_export_20260622.csv"
raw = pd.read_csv(source, sep=";", dtype=str, keep_default_na=False, encoding="utf-8-sig")
print(f"Source: {len(raw):,} records and {raw.shape[1]:,} columns")
print("Raw SHA256:", hashlib.sha256(source.read_bytes()).hexdigest())
"""),
md("## Cohort counts\n\nEvent and control numbers below are the existing provisional classification, not a frozen final model cohort."),
cc("""summary = []
cohorts = {}
for setting in ["SEH", "POLI"]:
    d = pd.read_csv(OUT / f"cohort_{setting}_cleaned.csv.gz", dtype=str, keep_default_na=False)
    outcome = pd.read_csv(OUT / f"outcomes_{setting}_existing_definition.csv", dtype=str, keep_default_na=False)
    status = pd.read_csv(OUT / f"cohort_{setting}_analysis_missing_status.csv.gz", dtype=str, keep_default_na=False)
    expected = raw.loc[raw.episodes.eq("1") & raw[f"{setting}_1stecontact"].eq("1"), "Participant Id"]
    assert d.participant_id.is_unique and d.participant_id.tolist() == expected.tolist()
    assert outcome.participant_id.tolist() == d.participant_id.tolist()
    assert status.participant_id.tolist() == d.participant_id.tolist()
    assert d.columns.tolist() == status.columns.tolist()
    counts = pd.to_numeric(outcome.gallstone_related_disease, errors="coerce")
    summary.append({"setting": setting, "first_presentations": len(d),
                    "provisional_positive": int(counts.eq(1).sum()),
                    "provisional_control": int(counts.eq(0).sum()),
                    "binary_outcome_unassigned": int(counts.isna().sum())})
    cohorts[setting] = d
print(pd.DataFrame(summary).to_string(index=False))
assert set(cohorts["SEH"].participant_id).isdisjoint(cohorts["POLI"].participant_id)
print("Verified exact source membership, unique record keys, aligned statuses and outcome joins.")
"""),
md("## Variable availability\n\nComplete per-variable tables are saved in missingness_SEH_all_first_visits.csv and missingness_POLI_all_first_visits.csv, with additional tables for the provisional primary subset. Missingness reasons remain separate; blank does not automatically mean imputable."),
cc("""for setting in ["SEH", "POLI"]:
    availability = pd.read_csv(OUT / f"missingness_{setting}_all_first_visits.csv")
    dictionary = pd.read_csv(OUT / f"variable_dictionary_{setting}.csv")
    print(setting, "availability rows:", len(availability), "dictionary rows:", len(dictionary))
    print("Availability columns:", ", ".join(availability.columns[:9]))
    print(availability.head(5).iloc[:, :8].to_string(index=False))
"""),
md("## Outstanding source review\n\nReview flags preserve the existing coding; they are not automatic exclusions. No source identity is resolved using demographic similarity alone."),
cc("""a = json.loads((OUT / "outcome_ascertainment_summary.json").read_text())
for setting, item in a["settings"].items():
    print(setting, "ascertainment/date/code review records:", item["records_flagged_for_review_n"])
    print("Controls with no definitive diagnosis:", item["review_reason_counts_controls"]["control_no_definitive_diagnosis_flag_0"])
pairs = pd.read_csv(MASTER / "candidate_duplicate_first_episodes.csv")
print("Possible duplicate pairs:", len(pairs))
flags = pd.read_csv(OUT / "quality_flags.csv")
print("Quality flag rows:", len(flags))
print(flags.groupby("setting").size().to_string())
"""),
md("## Verification\n\nThe full master retains all records and original variables. The analytical cohorts retain source snapshots, provide cell-level status matrices, and separate outcome information."),
cc("""master_check = json.loads((MASTER / "independent_validation.json").read_text())
assert master_check["independent_master_validation"] == "passed"
assert meta["imputation_performed"] is False
assert meta["predictor_selection_performed"] is False
for setting, values in meta["cohorts"].items():
    assert all(values["validation"].values())
    print(setting, values["validation"])
print("Master validation:", master_check["independent_master_validation"])
print("Imputation performed:", meta["imputation_performed"])
print("Clinical predictor selection performed:", meta["predictor_selection_performed"])
"""),
md("## Next decisions\n\n" + "\n\n".join(limitations[3:]) + "\n\nChoose and verify candidate variables jointly with the statistician after this stage. Multiple imputation must use the relevant predictors and outcome, preserve structural absence, and be integrated with internal validation. References: [MICE](https://amices.org/mice/reference/mice.html), [pooling estimates](https://amices.org/mice/reference/pool.html).")
]
namespace = {}
count = 0
for cell in cells:
    if cell.cell_type != "code":
        continue
    count += 1
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        exec(compile(cell.source, f"notebook-cell-{count}", "exec"), namespace)
    cell.execution_count = count
    cell.outputs = [nbformat.v4.new_output("stream", name="stdout", text=buffer.getvalue())]
notebook = nbformat.v4.new_notebook(cells=cells)
notebook.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": sys.version.split()[0]},
    "validation": {"execution": "all_code_cells_executed_in_order", "engine": "Python shared namespace",
                   "schema": "nbformat.validate", "jupyter_kernel_execution": False}
}
nbformat.validate(notebook)
nbpath = ROOT / "notebooks/02_operate_first_visit_cleaning.ipynb"
nbformat.write(notebook, nbpath)
nbformat.validate(nbformat.read(nbpath, as_version=4))
links = "".join('<li><a href="' + html.escape(path, quote=True) + '">' + html.escape(label) + '</a></li>' for label, path in files)
bullets = "".join("<li>" + html.escape(t) + "</li>" for t in limitations)
doc = """<!doctype html><html lang="en"><head><meta charset="utf-8"><title>OPERATE first-visit data</title>
<style>body{font:16px/1.55 Arial,sans-serif;color:#172b35;max-width:1050px;margin:48px auto;padding:0 28px}
h1{font-size:30px;line-height:1.2}h2{font-size:22px;margin-top:32px}table{border-collapse:collapse;width:100%;font-size:15px}
th,td{text-align:left;padding:10px;border-bottom:1px solid #ccd6da}th{background:#eef4f5}li{margin:10px 0}
a{color:#155d81}code{background:#eef4f5;padding:2px 5px}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f2f5f6;padding:18px}
.small{color:#52636b;font-size:14px}</style></head><body><h1>OPERATE first-visit data</h1>
<p class="small">30 September 2026 · Cleaning and missingness stage</p><p>""" + html.escape(scope) + "</p>" + table.to_html(index=False, border=0) + "<h2>Files</h2><ul>" + links + "</ul><h2>What needs review before imputation</h2><ul>" + bullets + """</ul><h2>How to use the package</h2>
<p>The two cleaned first-visit files contain the complete variable inventory for each setting.
Outcomes are separate and join by <code>participant_id</code>. The matching status matrices preserve the reason a cell is absent.
Keep the dictionary and source snapshots with the data. Read CSV strings without automatically converting the literal <code>NA</code> token.</p>
<p>No multiple imputation is delivered yet. The detailed protocol puts the variable/missingness review and joint predictor decision first.</p>
<h2>Validation and reproducibility</h2><p>The executed notebook checks source membership, unique record keys, join alignment, outcome counts and the saved validation results.
Its code cells executed top to bottom in Python and the file passed nbformat schema validation. The scripts and source hash are retained locally.</p>
<p>See <a href="README.txt">README.txt</a> for commands and import examples, and
<a href="../../analysis/first_visit_protocol_2026-09-30.txt">the protocol decisions</a> for scope.</p></body></html>"""
(OUT / "START_HERE.html").write_text(doc, encoding="utf-8")
# Text-only HTML companion for visually inspecting the saved notebook outputs.
preview = ["<!doctype html><meta charset='utf-8'><title>OPERATE notebook preview</title><style>body{font:16px/1.5 Arial;max-width:1100px;margin:36px auto}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f2f5f6;padding:18px}article{margin:28px 0}</style>"]
for cell in notebook.cells:
    if cell.cell_type == "markdown":
        preview.append("<article><pre>" + html.escape(cell.source) + "</pre></article>")
    else:
        preview.append("<details><summary>Code</summary><pre>" + html.escape(cell.source) + "</pre></details>")
        preview.extend("<pre>" + html.escape(o.get("text", "")) + "</pre>" for o in cell.outputs)
(OUT / "notebook_preview.html").write_text("\n".join(preview), encoding="utf-8")
print("Created and executed", nbpath)
print("Created", OUT / "START_HERE.html")
