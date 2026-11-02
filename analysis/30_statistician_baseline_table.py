"""Exhaustive observed-data first-presentation table for statistician handover.

Inputs: corrected first_visit_2026-10-27_v3 clean/status/dictionaries and the
existing frozen_analysis_cohort.csv. Eligibility/outcome are never re-derived.
One row per frozen record within setting; no imputation, fitting, hypothesis
tests, predictor selection, new exclusions, or code-label inference.
Outputs are aggregate, local handover files. Raw data/releases are read only.
Run with bundled Python: python analysis/30_statistician_baseline_table.py
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import html
import json
import re
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/first_visit_2026-10-27_v3"
FROZEN = ROOT / "results/analysis_2026-10-27_no_adjudication/frozen_analysis_cohort.csv"
OUT = ROOT / "outputs/statistician_handover_2026-11-02"
SUPPORT = OUT / "table_support"
SETTINGS = ("POLI", "SEH")
ROLE_LABELS = {
    "derived_index_measurement": "Derived presentation measurements",
    "baseline_or_first_visit_source_field": "History and presentation",
    "historical_information_with_timing_gate": "Historical information (timing gate)",
    "same_day_laboratory": "Same-day laboratory (availability order unproven)",
    "same_day_diagnostic_information": "Same-day investigations (availability order unproven)",
    "same_day_questionnaire": "Same-day questionnaire",
    "assessment_or_care_process_review": "Assessment or treatment (suitability review)",
    "undated_context_timing_unconfirmed": "Undated context (timing unconfirmed)",
    "reference_standard_timing_review": "Reference standard (not a baseline predictor)",
    "administrative_or_selection": "Administration or cohort selection",
    "administrative_assessor_metadata": "Assessor metadata",
    "linkage_audit": "Linkage audit",
}
ROLE_ORDER = {r: i for i, r in enumerate(ROLE_LABELS)}
DERIVED = {
    "age_years_at_index": ("Age at presentation", "years"),
    "years_since_bariatric_surgery": ("Time since bariatric surgery", "years"),
    "baseline_bmi_calculated": ("BMI before bariatric surgery (calculated)", "kg/m2"),
}


def read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)


def clean_label(value):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]*>", " ", value))).strip()


def json_value(value, default):
    return json.loads(value) if value else default


def fmt_number(value):
    return f"{value:.2f}".rstrip("0").rstrip(".")


def frac(n, den):
    return f"{n}/{den} ({100*n/den:.1f}%)" if den else "0/0 (no observations)"


def write_csv(name, records):
    pd.DataFrame(records).to_csv(OUT / name, index=False)


def main():
    SUPPORT.mkdir(parents=True, exist_ok=True)
    frozen = read(FROZEN)
    assert not frozen.duplicated(["setting", "participant_id"]).any()
    assert set(frozen.primaryeligibility) == {"0", "1"}
    eligible = frozen.loc[frozen.primaryeligibility.eq("1")]
    assert eligible.finalY.isin(["0", "1"]).all()
    assert not eligible.participant_id.duplicated().any(), "Eligible person repeats across settings"
    counts = {s: {"all": int(frozen.setting.eq(s).sum()),
                  "eligible": int(eligible.setting.eq(s).sum()),
                  "cases": int((eligible.setting.eq(s) & eligible.finalY.eq("1")).sum()),
                  "controls": int((eligible.setting.eq(s) & eligible.finalY.eq("0")).sum())}
              for s in SETTINGS}
    assert counts == {"POLI": {"all": 776, "eligible": 563, "cases": 88, "controls": 475},
                      "SEH": {"all": 789, "eligible": 567, "cases": 88, "controls": 479}}

    # 1. Exact keyed alignment, with a complete physical-column inventory.
    data, statuses, masks, fields, inventory = {}, {}, {}, {}, []
    inputs = [FROZEN]
    physical_coverage = {}
    for setting in SETTINGS:
        paths = [SOURCE / f"cohort_{setting}_cleaned.csv.gz",
                 SOURCE / f"cohort_{setting}_analysis_missing_status.csv.gz",
                 SOURCE / f"variable_dictionary_source_{setting}.csv",
                 SOURCE / f"variable_dictionary_analysis_{setting}.csv"]
        inputs.extend(paths)
        clean, status, source_dict, analytic_dict = map(read, paths)
        assert clean.participant_id.is_unique and status.participant_id.is_unique
        assert list(clean.columns) == list(status.columns)
        assert clean.participant_id.equals(status.participant_id)
        f = frozen.loc[frozen.setting.eq(setting)].set_index("participant_id")
        assert set(clean.participant_id) == set(f.index)
        f = f.loc[clean.participant_id].reset_index()
        data[setting], statuses[setting] = clean, status
        masks[setting] = {"all": pd.Series(True, index=clean.index),
                          "eligible": f.primaryeligibility.eq("1"),
                          "cases": f.primaryeligibility.eq("1") & f.finalY.eq("1"),
                          "controls": f.primaryeligibility.eq("1") & f.finalY.eq("0")}
        sd = source_dict.set_index("variable").to_dict("index")
        ad = analytic_dict.set_index("variable").to_dict("index")
        assert source_dict.variable.is_unique and analytic_dict.variable.is_unique
        analytic_for_source = {r["source_variable"]: v for v, r in ad.items()
                               if r["type"] != "label"}
        fields[setting] = {}
        for v in clean.columns:
            source = sd.get(v, {})
            analytic = ad.get(v, {})
            source_v = analytic.get("source_variable", source.get("companion_of") or v)
            companion = analytic_for_source.get(source_v, "")
            source_names = json_value(analytic.get("source_names", source.get("source_columns", "")), [])
            role = analytic.get("role", source.get("role", "administrative_or_selection"))
            label = clean_label(analytic.get("question", "") or source.get("original_or_derived_name", "") or v)
            typ = analytic.get("type", source.get("metadata_type", "unresolved"))
            note = analytic.get("note", source.get("note", ""))
            quarantined = any(term in note.lower() for term in
                              ["quarantin", "unresolved measurement scale/unit", "no target or auxiliary use"])
            codes = json_value(analytic.get("option_codes_json", ""), {})
            unit = analytic.get("unit", source.get("unit", ""))
            if v in DERIVED:
                label, unit = DERIVED[v]
                role, typ = "derived_index_measurement", "numeric"
                disposition, why = "descriptive_summary", "Existing observed derived measurement"
            elif analytic and typ != "label" and not quarantined:
                if typ in {"numeric", "slider"} or (typ in {"radio", "checkbox"} and codes):
                    disposition, why = "descriptive_summary", "Verified analytical companion; all fields retained"
                else:
                    disposition, why = "availability_only", "Verified summary type or labels unavailable"
            elif analytic and typ == "label":
                disposition, why = "duplicate_label_inventory", "Human-label mirror; analytic companion summarized once"
            elif analytic and quarantined:
                disposition, why = "quarantined_availability", note
            elif companion or source.get("companion_of"):
                disposition, why = "source_companion_inventory", "Source/qualifier copy retained; verified analytical companion owns summary"
            elif role in {"administrative_or_selection", "administrative_assessor_metadata", "linkage_audit"}:
                disposition, why = "administrative_inventory", "Administrative, identifying, selection, or linkage field; availability only"
            elif role == "reference_standard_timing_review":
                disposition, why = "reference_standard_inventory", "Reference-standard field, not a presentation predictor; availability only"
            else:
                disposition = "availability_only"
                why = "Source field retained without verified analytical companion; no values or guessed code meanings displayed"
            if typ == "checkbox" and source_names:
                # The option text comes from the verified source name, not code order.
                option = source_names[0].split("#", 1)
                if len(option) > 1:
                    label += ": " + clean_label(option[1])
            info = dict(variable=v, label=label, unit=unit, type=typ, role=role,
                        source_names=source_names, source_variable=source_v,
                        companion=companion, disposition=disposition, reason=why,
                        codes=codes, note=note, form=analytic.get("form", ""),
                        dependency_parent=analytic.get("dependency_parent", ""),
                        dependency_operator=analytic.get("dependency_operator", ""),
                        dependency_value=analytic.get("dependency_value", ""),
                        metadata_status=analytic.get("metadata_status", source.get("metadata_status", "")))
            fields[setting][v] = info
            primary = masks[setting]["eligible"]
            nonblank = clean[v].ne("")
            inventory.append(dict(setting=setting, variable=v, question=label, unit=unit,
                field_role=role, metadata_type=typ, disposition=disposition,
                summarized_companion=companion, disposition_reason=why,
                primary_n=int(primary.sum()), primary_nonblank_n=int((primary & nonblank).sum()),
                primary_blank_n=int((primary & ~nonblank).sum()), all_n=len(clean),
                all_nonblank_n=int(nonblank.sum()), all_blank_n=int((~nonblank).sum()),
                source_names=json.dumps(source_names, ensure_ascii=False),
                metadata_status=info["metadata_status"], dependency_parent=info["dependency_parent"],
                dependency_operator=info["dependency_operator"], dependency_value=info["dependency_value"],
                caution="Nonblank source counts do not establish usable analytical observations" if disposition != "descriptive_summary" else note))
        physical_coverage[setting] = dict(cleaned_columns=len(clean.columns),
            inventoried_columns=len(fields[setting]),
            source_dictionary_fields=len(sd), analysis_dictionary_fields=len(ad),
            omitted_columns=sorted(set(clean.columns) - set(fields[setting])),
            unexpected_dictionary_fields=sorted((set(sd) | set(ad)) - set(clean.columns)),
            disposition_counts=dict(Counter(x["disposition"] for x in fields[setting].values())))
        assert not physical_coverage[setting]["omitted_columns"]
        assert not physical_coverage[setting]["unexpected_dictionary_fields"]

    # 2. One display item per analytical field; unresolved source fields remain visible.
    # Pool only identically named, identically defined fields. Setting-specific
    # questions retain separate items; their other-setting cells are structural N/A.
    items = defaultdict(dict)
    for setting in SETTINGS:
        for v, info in fields[setting].items():
            if info["disposition"] not in {"descriptive_summary", "availability_only", "quarantined_availability"}:
                continue
            signature = (v, info["label"], info["unit"], info["type"], info["role"],
                         json.dumps(info["codes"], sort_keys=True), info["disposition"])
            items[signature][setting] = info
    ordered = sorted(items.items(), key=lambda x: (ROLE_ORDER.get(x[0][4], 99),
                     0 if x[0][6] == "descriptive_summary" else 1,
                     0 if x[0][0] in DERIVED else 1, x[0][0]))
    stats, missing_reasons, availability = [], [], []
    table_rows, within_rows, all_rows = [], [], []
    table_index = []
    for item_no, (signature, by_setting) in enumerate(ordered, 1):
        v, label, unit, typ, role, codes_json, disposition = signature
        settings = list(by_setting)
        item_id = f"F{item_no:04d}"
        codes = json.loads(codes_json)
        is_summary = disposition == "descriptive_summary"
        # Do not interpret unknown original source values. Count nonblank only.
        numeric = is_summary and typ in {"numeric", "slider"}
        categories = list(codes) if is_summary and not numeric else []
        summaries = {}
        for setting in settings:
            z = data[setting][v]
            for population in ("eligible", "cases", "controls", "all"):
                mask = masks[setting][population]
                observed = mask & z.ne("")
                n, obs = int(mask.sum()), int(observed.sum())
                reason_counts = statuses[setting].loc[mask & ~z.ne(""), v].value_counts().to_dict()
                rec = dict(item_id=item_id, setting=setting, population=population,
                    variable=v, question=label, unit=unit, role=role,
                    disposition=disposition, total_n=n, observed_n=obs,
                    missing_n=n-obs, median=None, q1=None, q3=None, categories={})
                if numeric:
                    values = pd.to_numeric(z.loc[observed], errors="raise")
                    assert np.isfinite(values).all(), v
                    if obs:
                        q = values.quantile([.25, .5, .75], interpolation="linear")
                        rec.update(median=float(q.loc[.5]), q1=float(q.loc[.25]), q3=float(q.loc[.75]))
                elif is_summary:
                    observed_codes = z.loc[observed].value_counts().to_dict()
                    assert not set(observed_codes) - set(codes), (v, set(observed_codes) - set(codes))
                    rec["categories"] = {c: int(observed_codes.get(c, 0)) for c in categories}
                    assert sum(rec["categories"].values()) == obs
                summaries[(setting, population)] = rec
                availability.append({k: val for k, val in rec.items() if k not in {"categories", "median", "q1", "q3"}})
                for reason, count in reason_counts.items():
                    missing_reasons.append(dict(item_id=item_id, setting=setting, population=population,
                        variable=v, missing_reason=reason or "blank_status", n=int(count)))
                assert sum(reason_counts.values()) == n - obs
                if is_summary:
                    if numeric:
                        stats.append({k: val for k, val in rec.items() if k != "categories"})
                    else:
                        for c in categories:
                            stats.append({**{k: val for k, val in rec.items() if k != "categories"},
                                "category_code": c, "category_label": codes[c],
                                "category_n": rec["categories"][c],
                                "percent_among_observed": 100*rec["categories"][c]/obs if obs else None})

        def combine(group):
            pop = group if group in {"cases", "controls"} else "eligible"
            ss = [group] if group in SETTINGS else settings
            if not set(ss).intersection(settings):
                return None
            selected = [summaries[(s, pop)] for s in ss if s in settings]
            result = dict(total_n=sum(r["total_n"] for r in selected),
                observed_n=sum(r["observed_n"] for r in selected),
                missing_n=sum(r["missing_n"] for r in selected), categories={})
            if numeric:
                values = pd.concat([pd.to_numeric(data[s].loc[masks[s][pop] & data[s][v].ne(""), v], errors="raise")
                                    for s in ss if s in settings], ignore_index=True)
                q = values.quantile([.25, .5, .75], interpolation="linear") if len(values) else None
                result.update(median=float(q.loc[.5]) if q is not None else None,
                              q1=float(q.loc[.25]) if q is not None else None,
                              q3=float(q.loc[.75]) if q is not None else None)
            else:
                result["categories"] = {c: sum(r["categories"].get(c, 0) for r in selected) for c in categories}
            return result

        def display(rec, kind, category=None):
            if rec is None:
                return "N/A (other setting)"
            if kind == "numeric":
                return (f"{fmt_number(rec['median'])} [{fmt_number(rec['q1'])}, {fmt_number(rec['q3'])}]"
                        if rec["observed_n"] else "No observations")
            if kind == "category":
                return frac(rec["categories"].get(category, 0), rec["observed_n"])
            if kind in {"observed", "nonblank"}:
                return frac(rec["observed_n"], rec["total_n"])
            return frac(rec["missing_n"], rec["total_n"])

        lines = []
        if numeric:
            lines += [("Median [Q1, Q3]", "numeric", None)]
        elif is_summary:
            lines += [(f"{codes[c]} [code {c}]", "category", c) for c in categories]
        lines += [("Observed / total" if is_summary else "Nonblank source / total (availability only)",
                   "observed" if is_summary else "nonblank", None),
                  ("Missing / total" if is_summary else "Blank source / total", "missing", None)]
        applicable = ", ".join(settings)
        caution = ("Only originating setting contributes to case/control denominator. " if len(settings) == 1 else "")
        if not is_summary:
            caution += "Source nonblank counts are not usable observed values; values withheld."
            if disposition == "quarantined_availability":
                caution += " " + next(iter(by_setting.values()))["note"]
        elif role.startswith("same_day"):
            caution += "Same date does not prove results were available before assessment."
        table_index.append(dict(item_id=item_id, variable=v, question=label, applicable_settings=applicable,
            role=role, disposition=disposition, primary_excel_first_row=len(table_rows)+8,
            primary_excel_last_row=len(table_rows)+len(lines)+7))
        for metric, kind, cat in lines:
            meta = [applicable, v, ROLE_LABELS.get(role, role), unit, item_id, caution]
            table_rows.append([label, metric] + [display(combine(g), kind, cat) for g in ["POLI", "SEH", "cases", "controls"]] + meta)
            within_rows.append([label, metric] + [display(summaries.get((s, p)), kind, cat)
                for s, p in [("POLI", "cases"), ("POLI", "controls"), ("SEH", "cases"), ("SEH", "controls")]] + meta)
            all_rows.append([label, metric] + [display(summaries.get((s, "all")), kind, cat) for s in SETTINGS] + meta)

    # 3. Availability/status coverage includes every physical field, even mirrors,
    # identifiers and fields outside the predictor window. No source values copied.
    full_availability, full_reasons = [], []
    for setting in SETTINGS:
        for v in data[setting].columns:
            z = data[setting][v]
            for pop in ["eligible", "all"]:
                mask = masks[setting][pop]
                nonblank = mask & z.ne("")
                full_availability.append(dict(setting=setting, population=pop, variable=v,
                    total_n=int(mask.sum()), nonblank_n=int(nonblank.sum()), blank_n=int((mask & z.eq("")).sum()),
                    disposition=fields[setting][v]["disposition"]))
                for reason, n in statuses[setting].loc[mask, v].value_counts().items():
                    full_reasons.append(dict(setting=setting, population=pop, variable=v,
                        status=reason or "blank_status", n=int(n)))
    assert len(inventory) == 4311
    assert len({(r["setting"], r["variable"]) for r in inventory}) == 4311
    assert all(r["total_n"] == r["nonblank_n"] + r["blank_n"] for r in full_availability)
    cohort = [{"setting": s, **counts[s]} for s in SETTINGS]
    cohort += [dict(setting="Pooled eligible", all=1565, eligible=1130, cases=176, controls=954),
               dict(setting="Unassigned (not described)", all=28, eligible=0, cases=0, controls=0)]
    write_csv("baseline_descriptive_statistics.csv", stats)
    write_csv("baseline_field_inventory.csv", inventory)
    write_csv("baseline_field_availability.csv", full_availability)
    write_csv("baseline_field_missing_status.csv", full_reasons)
    write_csv("baseline_display_field_index.csv", table_index)
    write_csv("baseline_cohort_counts.csv", cohort)
    headings = ["Variable / question", "Statistic / category", "POLI (n=563)", "SEH (n=567)",
                "Cases pooled (n=176)", "Controls pooled (n=954)", "Applicable settings", "Variable code",
                "Field role", "Unit", "Field ID", "Caution"]
    for name, rows, headers in [
        ("baseline_primary_table.csv", table_rows, headings),
        ("baseline_within_setting_table.csv", within_rows, headings[:2] + ["POLI cases (n=88)", "POLI controls (n=475)", "SEH cases (n=88)", "SEH controls (n=479)"] + headings[6:]),
        ("baseline_all_presentations_table.csv", all_rows, headings[:2] + ["POLI (n=776)", "SEH (n=789)"] + headings[6:])]:
        with (OUT/name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f); writer.writerow(headers); writer.writerows(rows)
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    coverage = dict(status="passed", source_release="first_visit_2026-10-27_v3", cohort_counts=counts,
        primary_n=1130, cases=176, controls=954, unassigned_first_presentations=28,
        columns=physical_coverage, total_physical_fields=4311, omitted_fields=0,
        displayed_items=len(table_index), displayed_verified_summary_items=sum(x[0][6] == "descriptive_summary" for x in ordered),
        displayed_availability_items=sum(x[0][6] != "descriptive_summary" for x in ordered),
        descriptive_rows=len(table_rows), input_sha256=hashes,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        quantile_method="pandas linear interpolation, Q1/median/Q3 of observed exact numeric values",
        code_labels="Only existing verified analysis dictionary labels; unknown values not reproduced",
        join="One-to-one participant_id within frozen setting; exact full set equality, no inflation",
        imputation=False, models_refitted=False, p_value_screening=False,
        restrictions=["Frozen eligibility and finalY unchanged", "No text narratives, IDs, individual dates or source values exported",
          "Pooled case/control groups overlap setting groups and are not an independent four-group comparison",
          "Setting-specific items use only originating-setting denominators in pooled case/control columns",
          "Source nonblank counts include preserved unavailable/unverified encodings and do not establish analytic availability",
          "Missing status reasons remain separated; blank never recoded as no",
          "Same-day investigations are not automatically eligible predictors; precise availability order unproven"])
    (OUT/"baseline_coverage_validation.json").write_text(json.dumps(coverage, indent=2, ensure_ascii=False))
    gap_path = OUT / "review/current_source_fields_outside_release.csv"
    metadata_path = OUT / "review/metadata_coverage_summary.json"
    source_gaps = read(gap_path).to_dict("records") if gap_path.exists() else []
    metadata_coverage = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    if source_gaps:
        assert len(source_gaps) == 56
        coverage["source_metadata_coverage"] = metadata_coverage
        (OUT/"baseline_coverage_validation.json").write_text(json.dumps(coverage, indent=2, ensure_ascii=False))
    payload = dict(primary_headers=headings, primary_rows=table_rows, within_rows=within_rows,
        all_rows=all_rows, inventory=inventory, field_index=table_index, cohort=cohort, coverage=coverage,
        source_gaps=source_gaps, metadata_coverage=metadata_coverage,
        missing_reason_totals=dict(Counter(r["status"] for r in full_reasons)))
    (SUPPORT/"baseline_workbook_data.json").write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False))
    print(json.dumps({k: coverage[k] for k in ["status", "primary_n", "cases", "controls", "total_physical_fields", "omitted_fields", "displayed_items", "displayed_verified_summary_items", "displayed_availability_items", "descriptive_rows"]}))


if __name__ == "__main__":
    main()
