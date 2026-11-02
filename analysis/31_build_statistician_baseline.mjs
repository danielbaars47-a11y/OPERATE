// Aggregate Table 1 workbook. Run from OPERATE using the bundled runtime.
// Place a symlink to this builder in the output directory and invoke Node with
// --preserve-symlinks-main so @oai/artifact-tool resolves from that directory's
// node_modules symlink to the bundled dependency tree. No project dependencies.
import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root = process.env.OPERATE_ROOT || process.cwd();
const out = process.env.OPERATE_HANDOVER_OUT || path.join(root, 'outputs/statistician_handover_2026-11-02');
const support = path.join(out, 'table_support');
const data = JSON.parse(await fs.readFile(path.join(support, 'baseline_workbook_data.json'), 'utf8'));
await fs.mkdir(support, {recursive: true});
const wb = Workbook.create();
const navy = '#23435F', text = '#243442', pale = '#EAF0F5', gray = '#687784';
function col(n) { let s=''; for(n++;n;n=Math.floor((n-1)/26)) s=String.fromCharCode(65+(n-1)%26)+s; return s; }
function safe(v) { return typeof v==='string' && /^[=+@]/.test(v) ? "'"+v : v; }
function base(name, title, context, note, headers, rows, widths, tableName, groupCol=null) {
  const sheet=wb.worksheets.add(name), last=rows.length+7, end=col(headers.length-1);
  sheet.showGridLines=false;
  const area=sheet.getRange(`A1:${end}${last}`);
  area.format.font={name:'Arial',size:10,color:text};
  area.format.verticalAlignment='center';
  sheet.getRange('A2').values=[[title]];
  sheet.getRange('A2').format.font={name:'Arial',size:15,bold:true,color:navy};
  sheet.getRange(`A2:${end}2`).format.rowHeight=26;
  sheet.getRange('A3').values=[[context]];
  sheet.getRange('A4').values=[[note]];
  sheet.getRange('A5').values=[['Observed data only. Median [Q1, Q3]; categories n/N observed (%). Missing / total shown separately.']];
  sheet.getRange('A3:A5').format.font={name:'Arial',size:10,color:gray};
  sheet.getRange('A3:A5').format.wrapText=false;
  sheet.getRange(`A6:${end}6`).format.borders={bottom:{style:'thin',color:navy}};
  sheet.getRange(`A7:${end}7`).values=[headers];
  sheet.getRange(`A8:${end}${last}`).values=rows.map(r=>r.map(safe));
  const table=sheet.tables.add(`A7:${end}${last}`,true,tableName);
  table.style='TableStyleMedium2'; table.showFilterButton=true;
  sheet.getRange(`A7:${end}7`).format={fill:navy,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},
    wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:40};
  sheet.getRange(`A8:${end}${last}`).format.wrapText=true;
  widths.forEach((w,i)=>sheet.getRange(`${col(i)}1:${col(i)}${last}`).format.columnWidth=w);
  for (let i=0;i<rows.length;i++) {
    let lines=1;
    for(let j=0;j<rows[i].length;j++) {
      const val=String(rows[i][j]??'');
      const length=val.split('\n').reduce((n,line)=>n+Math.max(1,Math.ceil(line.length/Math.max(7,widths[j]*0.95))),0);
      lines=Math.max(lines,length);
    }
    sheet.getRange(`A${i+8}:${end}${i+8}`).format.rowHeight=Math.max(24,lines*13+8);
    if(groupCol!==null && (i===0 || rows[i][groupCol]!==rows[i-1][groupCol])) {
      sheet.getRange(`A${i+8}:${end}${i+8}`).format.borders={top:{style:'thin',color:'#BAC8D4'}};
      sheet.getRange(`A${i+8}:B${i+8}`).format.font={name:'Arial',size:10,bold:true,color:text};
    }
  }
  sheet.freezePanes.freezeRows(7); sheet.freezePanes.freezeColumns(2);
  sheet.tabColor=navy;
  return sheet;
}

const widths=[62,40,25,25,25,25,17,47,42,20,10,64];
const context='Corrected first-visit release, 27 October 2026. Frozen eligible cohort: 1,130 records, 176 cases and 954 controls.';
const note='POLI/SEH and pooled case/control columns overlap. Setting-specific rows use only the originating setting in pooled denominators.';
const primary=base('Baseline characteristics','OPERATE baseline characteristics',context,note,
  data.primary_headers,data.primary_rows,widths,'BaselineCharacteristics',10);
primary.getRange(`C8:F${data.primary_rows.length+7}`).format.horizontalAlignment='right';
const withinHeaders=data.primary_headers.slice();
withinHeaders.splice(2,4,'POLI cases (n=88)','POLI controls (n=475)','SEH cases (n=88)','SEH controls (n=479)');
const within=base('Within setting','Case and control characteristics within setting',context,
  'Cases and controls are separated within POLI and SEH. N/A denotes a field not included for the other setting.',
  withinHeaders,data.within_rows,widths,'WithinSettingCharacteristics',10);
within.getRange(`C8:F${data.within_rows.length+7}`).format.horizontalAlignment='right';
const allHeaders=data.primary_headers.slice(0,2).concat(['POLI (n=776)','SEH (n=789)'],data.primary_headers.slice(6));
const allWidths=widths.slice(0,4).concat(widths.slice(6));
const all=base('All first presentations','All assigned first presentations',
  'POLI 776 and SEH 789 records, including records outside the frozen eligible analysis cohort.',
  'Availability audit. No case/control classification is assigned to ineligible records. The 28 unassigned records are not tabulated.',
  allHeaders,data.all_rows,allWidths,'AllFirstPresentationCharacteristics',8);
all.getRange(`C8:D${data.all_rows.length+7}`).format.horizontalAlignment='right';

const ikeys=['setting','variable','question','unit','field_role','metadata_type','disposition','summarized_companion',
  'disposition_reason','primary_n','primary_nonblank_n','primary_blank_n','all_n','all_nonblank_n','all_blank_n',
  'source_names','metadata_status','dependency_parent','dependency_operator','dependency_value','caution'];
const iheaders=['Setting','Physical field','Question / source label','Unit','Field role','Metadata type','Disposition',
  'Analytical companion','Disposition reason','Eligible total','Eligible nonblank','Eligible blank','All total','All nonblank','All blank',
  'Source names','Metadata status','Dependency parent','Dependency operator','Dependency value','Caution'];
const irows=data.inventory.map(r=>ikeys.map(k=>r[k]??''));
const inventory=base('Field inventory','Complete corrected-release field inventory',
  '4,311 setting-fields: POLI 2,008 and SEH 2,303. Every physical cleaned column has exactly one disposition.',
  'Nonblank source values may contain unavailable or unverified encodings. These counts do not establish usable analytical observations.',
  iheaders,irows,[13,48,62,18,46,19,38,48,85,19,21,19,17,20,18,90,43,43,23,30,90],'CorrectedReleaseFieldInventory');
inventory.getRange('A5').values=[['Mirrors, raw values, numeric bounds and qualifiers remain here. Verified analytical companions own the descriptive summary once.']];
inventory.getRange(`J8:O${irows.length+7}`).setNumberFormat('#,##0');
inventory.getRange(`J8:O${irows.length+7}`).format.horizontalAlignment='right';
inventory.tabColor='#708597';

const gapKeys=['source_name','metadata_variable_name','visit','form','question','field_type','unit',
  'verification_state','disposition','reason','available_at_presentation','observed_n','followup_action'];
const gapHeaders=['Current source field','Metadata field','Visit','Form','Question','Field type','Unit',
  'Metadata verification','Release disposition','Reason','Presentation availability','Observed count','Required review'];
const gapRows=(data.source_gaps||[]).map(r=>gapKeys.map(k=>r[k]??''));
if(gapRows.length) {
  const gaps=base('Source coverage','Current source fields outside the corrected release',
    '56 current mapped diagnostic fields are not physically retained in the corrected first-visit release.',
    '40 general-investigation fields require source/timing QC; 16 colonoscopy fields were not materialized without same-day matches.',
    gapHeaders,gapRows,[50,50,30,43,65,20,20,45,40,100,38,23,110],'CurrentSourceCoverageGaps');
  gaps.getRange('A5').values=[['The separate metadata audit records 1,072 historical fields without a verified current source mapping. Their presence, equivalence and availability are not inferred.']];
  gaps.tabColor='#708597';
}

const methodsRows=[
  ['Population','Frozen eligible first-presentation records only for the main table. POLI: 563, including 88 cases and 475 controls. SEH: 567, including 88 cases and 479 controls.'],
  ['Outcome and exclusions','Existing frozen primaryeligibility and finalY are authoritative. No outcomes, exclusions, identity resolutions, clinical adjudications or cohort membership were changed.'],
  ['Unit and linkage','One frozen record per eligible person; one-to-one participant_id alignment within setting. Full identifier-set equality with clean/status data was checked. No identifiers or participant-level observations are included in this workbook.'],
  ['Overlapping columns','POLI and SEH partition the eligible cohort. Cases and controls pool the two settings. The four columns are overlapping descriptions and do not define four independent comparison groups.'],
  ['Setting-specific fields','A question specific to one setting remains a separate row. The other setting is structural N/A. Case/control denominators for these rows use the originating setting only, even though column headers give the full pooled group size.'],
  ['Field scope','Every physical field in the corrected 2026-10-27_v3 first-visit release is inventoried. This is release completeness, rather than an assertion that every field in every Castor CRF was exported or is an eligible predictor.'],
  ['Source coverage gaps','Source coverage lists 56 current mapped diagnostic fields outside the corrected release. The metadata audit also records 1,072 historical fields lacking a verified current source mapping. No source match, baseline timing or availability is inferred for these historical fields. See historical_firstvisit_metadata_coverage.csv and metadata_coverage_summary.json in the review evidence.'],
  ['All fields',`There is no preferred predictor subset, rarity filter, missingness filter, univariable p-value screening or significance testing. All ${data.coverage.displayed_verified_summary_items} distinct verified summary items and ${data.coverage.displayed_availability_items} availability-only items are displayed.`],
  ['Duplicates','Raw, numeric-bound, numeric-qualifier and label companions remain in Field inventory. A verified analytical companion owns the descriptive distribution once. Administration, linkage and reference-standard fields receive inventory availability only.'],
  ['Numerical variables','Median [Q1, Q3] uses observed exact analytical values. Quartiles use pandas linear interpolation. Observed / total and missing / total are shown separately. Censored, qualified, quarantined and unresolved source values are not silently treated as exact measurements.'],
  ['Categorical variables','Verified dictionary option labels and explicit codes are used. Counts are n / N observed (%), where N observed excludes blank analytical values. All documented categories, including zero counts, are retained. Checkbox options are conditional on an answered group and are not mutually exclusive across options.'],
  ['Availability-only fields','Only source nonblank and source blank counts are displayed. They are not usable observed clinical measurements. Unknown encodings, source text, dates and treatment/assessment fields are retained for suitability review without copying their participant values.'],
  ['Missingness','Blank analytical data are not recoded as absent or no. Field-level status counts are retained in baseline_field_missing_status.csv. Dependency not triggered, unanswered, unasked, not applicable, unperformed, failure, unknown, unresolved coding, timing uncertainty and quarantine remain distinct source statuses.'],
  ['Imputation','The entire workbook describes original observed values. No imputed datasets were loaded and no outcomes, settings or structural missingness were imputed.'],
  ['Timing and leakage','Same-day diagnostic/laboratory or questionnaire data are optional descriptive information, with field roles retained. Matching date does not prove availability before clinical assessment. Working diagnoses, treatment and reference-standard information are not automatically promoted to presentation predictors.'],
  ['All first presentations','This audit includes 776 POLI and 789 SEH records, including ineligible records. Twenty-eight frozen unassigned records remain excluded from setting summaries. Case/control classification is only used for eligible known-outcome records.'],
  ['Source release','first_visit_2026-10-27_v3: cohort_POLI/SEH_cleaned.csv.gz, cohort_POLI/SEH_analysis_missing_status.csv.gz and variable_dictionary_source/analysis_POLI/SEH.csv.'],
  ['Frozen manifest','analysis_2026-10-27_no_adjudication/frozen_analysis_cohort.csv. Source hashes and the complete column reconciliation are recorded in baseline_coverage_validation.json.'],
  ['Regeneration','analysis/30_statistician_baseline_table.py creates aggregate tables and checks; analysis/31_build_statistician_baseline.mjs creates this workbook with the bundled Artifact Tool runtime. Existing models and historical outputs are unchanged.'],
];
const methods=base('Methods','Definitions, scope and source provenance','Statistician handover prepared 2 November 2026.',
  'Observed descriptive information and release-completeness audit. No model performance or statistical validity claim is added.',
  ['Topic','Definition'],methodsRows,[30,135],'BaselineMethods');
methods.getRange('A5').values=[['Use Field inventory to inspect dispositions, metadata uncertainty, dependencies and the retained source-to-analysis mapping.']];
methods.freezePanes.unfreeze();methods.tabColor='#89929A';

wb.recalculate();
const inspection=await wb.inspect({kind:'table',range:'Baseline characteristics!A7:F13',include:'values,formulas',tableMaxRows:7,tableMaxCols:6,maxChars:5500});
await fs.writeFile(path.join(support,'workbook_inspection.ndjson'),inspection.ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',
  options:{useRegex:true,maxResults:30},maxChars:2500,summary:'Final workbook error scan'});
await fs.writeFile(path.join(support,'formula_scan.ndjson'),errors.ndjson);
for(const [name,range,file] of [
  ['Baseline characteristics','A2:F17','baseline'],
  ['Within setting','A2:F17','within_setting'],
  ['All first presentations','A2:D17','all_presentations'],
  ['Field inventory','A2:G14','inventory'],
  ['Methods','A2:B12','methods'],
  ['Field inventory','H7:O13','inventory_dispositions'],
  ...(gapRows.length ? [['Source coverage','A2:F13','source_coverage']] : []),
]) {
  const preview=await wb.render({sheetName:name,range,scale:1,format:'png'});
  await fs.writeFile(path.join(support,`${file}.png`),new Uint8Array(await preview.arrayBuffer()));
}
// Check a long question and later rows as well as the opening view.
const longest=data.primary_rows.reduce((best,r,i)=>r[0].length>data.primary_rows[best][0].length?i:best,0)+8;
for(const [range,file] of [[`A${longest}:F${longest+2}`,'long_question'],[`A${data.primary_rows.length+4}:F${data.primary_rows.length+7}`,'last_rows']]) {
  const preview=await wb.render({sheetName:'Baseline characteristics',range,scale:1,format:'png'});
  await fs.writeFile(path.join(support,`${file}.png`),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(out,'OPERATE_baseline_characteristics.xlsx'));
console.log(JSON.stringify({workbook:'OPERATE_baseline_characteristics.xlsx',sheets:gapRows.length?6:5,primaryRows:data.primary_rows.length,physicalFields:data.inventory.length,sourceGapFields:gapRows.length,renderedRanges:gapRows.length?9:8}));
