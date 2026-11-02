# OPERATE

OPERATE studies patients presenting with abdominal pain after bariatric surgery. This repository contains the code for Castor data cleaning, missing-data assessment, multiple imputation and diagnostic prediction of recorded definitive gallstone-related disease during the index episode.

**For the statistician:** start with the comprehensive observed-data baseline workbook and the one-page QC, imputation and results summary supplied separately through the approved research channel. They are not patient datasets and do not limit the candidate-variable inventory to the nine variables in the previous laboratory model. Clinical datasets, generated workbooks, model objects and research documents remain outside this public repository.

## What the project is trying to predict

- The unit is the first study episode per patient. The first presentation determines the setting: emergency department (**SEH**) or bariatric outpatient clinic (**POLI**).
- Follow-up covers the whole index episode until recorded closure after three months without contact. It is not a fixed 90-day window after presentation.
- The target is recorded definitive cholecystolithiasis or cholecystitis, considering all three definitive-diagnosis fields and the required supporting evidence. A procedure alone does not define a case. Follow-up information establishes the outcome, not baseline predictors.
- Primary controls have an assessable outcome without qualifying definitive gallstone disease; an unclear final diagnosis does not itself imply a missing outcome. Missing/uninterpretable outcome registration remains missing.
- Prior cholecystectomy excludes a patient from both cases and controls. Unresolved identity, eligibility and outcome issues remain marked in the full cleaned files and excluded from the relevant analysis population. No outcomes or setting labels are imputed.

This is a diagnostic model in a symptomatic population, distinct from the BARIA project on incident gallstones after surgery.

## Which information is retained

All baseline and first-presentation fields retained in the corrected first-visit release are inventoried: demographic and bariatric history, operation details, weight and weight loss, comorbidity, previous surgery, medication, substance use, history, symptoms, examination, vital signs, laboratory tests, imaging, other investigations and initial assessment or management. There is no univariable p-value screening or preferred-variable shortlist for this inventory. Coverage is reconciled to every physical field in that release. An additional source-coverage audit found all currently mapped baseline, SEH and POLI form fields represented, but 56 mapped other-investigation fields absent from the extracted release. Some lacked index-day matches; others were not extracted by the earlier diagnostic-block code. The handover inventories these gaps rather than inventing baseline values or silently selecting them out. Historical metadata fields without a verified current mapping cannot be assumed to be dropped current variables.

Retaining a field does not establish its suitability for a final prediction model. Administrative fields and identifiers are not predictors; exact dates support derived ages/intervals. Free text requires separate coding. Unverified answer codes and quarantined measurements remain unresolved. Initial working diagnoses and treatment plans need separate assessment of timing and leakage. Imaging is an optional extension because it can contribute to the diagnostic reference standard.

## How to read the comprehensive baseline table

The primary table uses the frozen outcome-eligible population. Its four requested columns are **POLI**, **SEH**, **cases** and **controls**. The first pair splits patients by setting; the second pair pools settings and splits the same patients by outcome. These are overlapping views, not four independent groups. Within-setting case/control tables are supplied alongside them because case mix and diagnostic availability differ substantially.

- Descriptions use observed cleaned values before imputation, so imputation cannot conceal missingness.
- Continuous measurements are described by median and interquartile range. Categorical responses use counts with explicit observed denominators and percentages.
- For a setting-specific item, pooled case/control summaries use the originating setting's case/control denominators, not the larger pooled header totals.
- Not performed, not applicable, unknown/not described, genuinely missing and unresolved coding/timing are distinguished where the CRF supports the distinction. A blank is not automatically a negative response or proof of a test being ordered.
- The field inventory accounts for raw fields, analytical companions, labels and technical metadata. Validated companions avoid displaying duplicate versions as separate clinical characteristics; unresolved fields retain availability and their disposition.
- The full first-presentation overview also retains people outside the outcome-eligible population. Model-specific availability restrictions are additional and do not define the primary baseline table.

No hypothesis tests or association-based variable selection are needed to describe this inventory.

## Data releases and quality control

| Stage | Meaning | Main code |
|---|---|---|
| Original Castor export | Immutable patient source; versioned metadata checked against the current export and CRF | `13_import_castor_metadata.py` |
| Corrected first-visit v3 | First-episode/setting linkage, clean values, missingness reasons, endpoint evidence and exclusions | `14_prepare_first_visit_v3.py` |
| Preparatory blockwise imputation | Separate SEH/POLI MICE and release checks | `10_impute_first_visit_v2.R`, `15_validate_v3.py`, `15_validate_saved_mice_v3.R`, `16_package_first_visit_v3.py` |
| Frozen analysis cohort | Eligibility and model specifications under the decision to proceed without further adjudication | `23_freeze_analysis.py` |
| Model-specific inputs and estimation | Observed restricted populations, fresh imputation within validation, ridge models | `24_model_prepare_inputs.py`, `24_model_pipeline.R` |
| Statistician handover | Exhaustive observed baseline descriptions and field disposition inventory | `30_statistician_baseline_table.py`, `31_build_statistician_baseline.mjs` |

The active cleaned release is `results/first_visit_2026-10-27_v3`; saved model results are in `results/analysis_2026-10-27_no_adjudication`. These locations are private and are not bundled here. The initially withdrawn v3 output and earlier September/v2 releases remain historical only. File names can retain an older version suffix while operating on an explicitly supplied newer release: inspect arguments and configuration.

QC checks identities, linkage cardinality, record counts, code mappings, dependency gates, dates, numeric parsing, units, plausibility, missingness and index-visit timing. Implausible source coding is not repaired by guessing. ASA, mixed-scale urine specific gravity and unresolved-unit eGFR remain quarantined; out-of-range urine pH is also withheld analytically. Original source values are preserved privately. Same calendar date does not prove that a result was known before the initial assessment, and tests crossing midnight need an explicit visit link.

## Two different imputation stages

**Preparatory cleaning release.** Fifty imputations per setting use clinically defined form/applicability blocks, 30 iterations extended to 100 when needed, numeric predictive mean matching and supported categorical models. Known outcomes can be auxiliary variables in their stratum; people with missing outcomes form separate strata. Only permitted missing predictor cells can change. Unsupported or poorly converged cells remain empty in the released datasets. BMI is recalculated from weight and height rather than independently imputed. This blockwise release does not preserve every cross-block association and is not a universal completed dataset for any future model.

**Existing prediction models.** Model-specific imputation is rebuilt inside bootstrap training samples. The final models use 50 imputations; each bootstrap uses five. Only genuinely missing time since surgery is imputed in the existing model sets; the other selected measurements are observed by cohort restriction. Outcomes are excluded from this imputation, and evaluation rows cannot estimate the imputer or supply training donors. This differs from the broader preparatory imputation and requires review when specifying a new inferential model.

Structural nonmeasurement, such as outpatient laboratory tests that were never requested, is never turned into ordinary missing laboratory values.

## Existing analyses and how to interpret them

The previous SEH laboratory model contains age, sex, time since bariatric surgery, the documented RYGB code group, right upper abdominal pain, bilirubin, alkaline phosphatase, ALAT and CRP. Laboratories use **ln(1 + value)**. RYGB is a recorded-code grouping; the comparator includes an unspecified operation code and is not proven non-RYGB anatomy.

The adjusted model uses fixed ridge logistic regression with training-only standardisation and penalty `lambda = 1/n_training`, without automatic backward selection. Internal validation uses 200 bootstrap samples and reports optimism-corrected discrimination, Brier score and calibration. Penalised odds ratios have bootstrap percentile intervals. These models are existing provisional analyses, not a final decision on the exhaustive candidate-variable set.

The later univariable analyses fit each of the same nine predictors separately using ordinary logistic regression and Rubin pooling of 50 saved imputations. Their odds ratios are unadjusted and their confidence intervals use a different method from the ridge model. A change between crude and adjusted odds ratios is not automatically a coding error or a causal finding. Univariable results do not select predictors.

Imaging extensions are exploratory: the verified variables concern inflammatory findings, not a complete gallstone/duct imaging set. Inclusion in the reference standard creates incorporation bias. Availability selection, unresolved source adjudication, within-day timing and absence of external validation limit clinical interpretation.

## Script guide

| Scripts | Purpose |
|---|---|
| 01–07 | Initial cleaning, linkage/outcome audits and early handover |
| 08–12 | Historical v2 dictionary, preparation, imputation and checks |
| 13–18 | Metadata reconciliation, corrected v3, release/validation and source-review queues |
| 20–22 | Remaining decisions and per-patient clinical review workbook |
| 23–24 | Frozen cohorts, model inputs, imputation, fitting, bootstrap and saved-result processing |
| 25–26 | Earlier observed summaries and analytical report; the earlier core table is not the exhaustive inventory |
| 27 | Statistical checks of saved models and outputs |
| 28 | Forestplot regenerated from saved multivariable estimates |
| 29 | Exploratory univariable regressions |
| 30–31 | Comprehensive statistician baseline table and workbook |
| 32 | One-page document authoring from separately supplied reviewed prose |
| 33 | Read-only current-source and historical-metadata coverage audit |

## Running locally

This repository is a reviewable source snapshot, not a self-contained dataset or a single-command pipeline. Obtain the approved source export, metadata, frozen cohorts and saved results separately. Keep them in access-controlled research storage. Generated outputs must remain ignored by Git.

Clone into a folder named `OPERATE`; many scripts expect the working directory to be its parent, while others accept explicit project/release arguments. Read headers and configuration before execution. Some scripts execute on import/source. Do not run every script in numeric order or use historical files as current inputs.

Examples requiring private inputs:

```sh
python3 OPERATE/analysis/23_freeze_analysis.py
python3 OPERATE/analysis/24_model_prepare_inputs.py
Rscript OPERATE/analysis/24_model_pipeline.R OPERATE/results/analysis_2026-10-27_no_adjudication
Rscript OPERATE/analysis/29_SEH_lab_univariable.R OPERATE
```

Python cleaning dependencies are in `analysis/requirements-cleaning.txt`; metadata/report/workbook utilities also require the packages named in their imports. R uses `mice`, `glmnet`, `jsonlite` and dependencies. The recorded univariable run used R 4.6.1. Bundled library folders are excluded and this is not a fully locked environment. The workbook renderer uses `@oai/artifact-tool`, available in the originating workspace runtime. For the handover, run the read-only metadata audit (33) before descriptive preparation (30), then render the workbook (31) using the runtime instructions in its header. Script 32 requires the reviewed handover prose supplied privately. The shared copy of script 33 only adapts the original review helper's relative paths; its metadata audit logic is unchanged.

## Provenance and review

`docs/source_manifest.json` maps shared source files to SHA-256 hashes. The original code snapshot passed Python and R syntax parsing. The statistician handover adds descriptive calculations and a separate before/after statistical review; it does not refit the existing prediction models. A code/privacy review does not independently validate the clinical source data or missing-at-random assumptions.

The workbook and one-page document are supplied separately for senior statistical review. They describe observed availability and the existing results, preserving null findings and unresolved decisions. No clinical readiness or external validation is claimed.
