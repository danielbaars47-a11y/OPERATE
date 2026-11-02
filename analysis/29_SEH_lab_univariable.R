#!/usr/bin/env Rscript
# Exploratory unadjusted associations in the frozen SEH_lab population.
# Reuse 50 saved imputations, unchanged nine predictors and outcome. No MI,
# multivariable fit, selection, bootstrap or new transformation is performed.
# Each model: binomial logit with intercept and one slope. Pool slope and
# variance using Rubin's rules and Barnard-Rubin finite df (mice::pool).
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(if (length(args)) args[1] else 'OPERATE', mustWork = TRUE)
.libPaths(c(file.path(project, 'analysis/R-library'), .libPaths()))
stopifnot(requireNamespace('mice', quietly = TRUE))
source_dir <- file.path(project, 'results/analysis_2026-10-27_no_adjudication')
out <- file.path(project, 'results/univariable_2026-10-28')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = '')

# 1. Frozen inputs and exact ID-based linkage; outcomes never taken from MI.
input_path <- file.path(source_dir, 'model_input_SEH_lab.csv')
original <- read.csv(input_path, colClasses = c(participant_id = 'character'))
variables <- c('age', 'male', 'years_since_surgery', 'RYGB', 'RUQ_pain', 'bilirubin', 'ALP', 'ALT', 'CRP')
stopifnot(nrow(original) == 455L, !anyNA(original$participant_id),
          !anyDuplicated(original$participant_id), all(original$outcome %in% 0:1),
          sum(original$outcome) == 79L,
          identical(names(original), c('participant_id', 'outcome', variables)),
          sum(is.na(original[variables])) == 3L,
          sum(is.na(original$years_since_surgery)) == 3L)
paths <- file.path(source_dir, sprintf('SEH_lab_completed_%02d.csv.gz', 1:50))
stopifnot(all(file.exists(paths)))
checks <- list()
draws <- lapply(seq_along(paths), function(i) {
  d <- read.csv(paths[i], colClasses = c(participant_id = 'character'))
  stopifnot(!anyNA(d$participant_id), !anyDuplicated(d$participant_id),
            nrow(d) == 455L, setequal(d$participant_id, original$participant_id),
            all(d$.imp == i), identical(names(d), c('participant_id', '.imp', variables)))
  d <- d[match(original$participant_id, d$participant_id), ]
  stopifnot(all(vapply(d[variables], is.numeric, logical(1))),
            all(is.finite(as.matrix(d[variables]))))
  for (v in variables) {
    observed <- !is.na(original[[v]])
    delta <- max(abs(d[[v]][observed] - original[[v]][observed]))
    stopifnot(delta < 1e-12)
    checks[[length(checks) + 1L]] <<- data.frame(imputation = i, predictor = v,
      observed_n = sum(observed), imputed_n = sum(!observed), max_observed_absolute_difference = delta)
  }
  stopifnot(all(d$male %in% 0:1), all(d$RYGB %in% 0:1), all(d$RUQ_pain %in% 0:1))
  d$outcome <- original$outcome
  d
})
write_csv(do.call(rbind, checks), 'input_checks.csv')

# 2. Ordinary maximum-likelihood single-predictor logistic regressions.
fit_checks <- list()
results <- lapply(variables, function(v) {
  fits <- lapply(seq_along(draws), function(i) {
    warnings <- character()
    fit <- withCallingHandlers(glm(reformulate(v, response = 'outcome'), data = draws[[i]],
      family = binomial(), na.action = na.fail), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
      })
    valid <- isTRUE(fit$converged) && fit$rank == 2L && !isTRUE(fit$boundary) &&
      all(is.finite(coef(fit))) && all(is.finite(vcov(fit))) && nobs(fit) == 455L
    fit_checks[[length(fit_checks) + 1L]] <<- data.frame(predictor = v, imputation = i,
      n = nobs(fit), events = sum(draws[[i]]$outcome), converged = fit$converged,
      boundary = fit$boundary, rank = fit$rank, iterations = fit$iter,
      min_fitted = min(fitted(fit)), max_fitted = max(fitted(fit)),
      warnings = paste(warnings, collapse = '; '), valid = valid)
    if (!valid) stop('Invalid fit: ', v, ', imputation ', i)
    fit
  })
  pooled <- mice::pool(fits, dfcom = 453L)
  z <- pooled$pooled[pooled$pooled$term == v, ]
  stopifnot(nrow(z) == 1L, is.finite(z$t), z$t > 0, is.finite(z$df), z$df > 0)
  se <- sqrt(z$t); critical <- qt(.975, df = z$df)
  data.frame(predictor = v, n = 455L, events = 79L, imputations = 50L,
    log_OR = z$estimate, pooled_SE = se, OR = exp(z$estimate),
    CI95_lower = exp(z$estimate - critical * se), CI95_upper = exp(z$estimate + critical * se),
    p_value = 2 * pt(-abs(z$estimate / se), df = z$df), df = z$df,
    within_variance = z$ubar, between_variance = z$b, total_variance = z$t,
    fraction_missing_information = z$fmi)
})
results <- do.call(rbind, results)
write_csv(results, 'SEH_lab_univariable.csv')
write_csv(do.call(rbind, fit_checks), 'fit_checks.csv')

# 3. Human-readable Dutch results, scales, provenance and limits.
labels <- c('Leeftijd (per jaar)', 'Man versus vrouw', 'Tijd sinds bariatrische operatie (per jaar)',
  'Gedocumenteerde RYGB-codegroep versus overige codes', 'Pijn rechter bovenbuik: ja versus nee',
  'Bilirubine (per eenheid ln[1+x], x in µmol/L)', 'AF (per eenheid ln[1+x], x in U/L)',
  'ALAT (per eenheid ln[1+x], x in U/L)', 'CRP (per eenheid ln[1+x], x in mg/L)')
fmt_p <- function(p) ifelse(p < .001, '<0,001', sub('\\.', ',', sprintf('%.3f', p)))
lines <- c('# Univariable logistische regressies: SEH met laboratoriumgegevens', '',
  'Dezelfde 455 deelnemers en 79 uitkomsten als het bestaande SEH_lab-model; 50 bestaande imputaties.', '',
  '| Predictor / vergelijking | OR | 95%-BI | p (ongecorrigeerd) |', '|---|---:|---:|---:|')
for (j in seq_len(nrow(results))) lines <- c(lines, sprintf('| %s | %.2f | %.2f–%.2f | %s |',
  labels[j], results$OR[j], results$CI95_lower[j], results$CI95_upper[j], fmt_p(results$p_value[j])))
lines <- c(lines, '', '## Methoden en controles', '',
  'Per predictor is afzonderlijk een gewone, onbestrafte logistische regressie met intercept geschat. OR’s zijn ongecorrigeerde associaties met de geregistreerde definitieve galsteengerelateerde uitkomst (1 versus 0), geen causale effecten. Alle negen vooraf gekozen predictoren zijn gerapporteerd; er is geen selectie op significantie uitgevoerd. P-waarden zijn tweezijdig en niet voor meervoudig toetsen gecorrigeerd; deze aanvullende analyses zijn exploratief.', '',
  'De 50 opgeslagen SEH_lab_completed-bestanden zijn gekoppeld via unieke deelnemer-ID aan de oorspronkelijke model_input_SEH_lab-uitkomst. Alle 455 ID’s, de 79 uitkomsten en alle oorspronkelijk waargenomen predictorwaarden zijn gecontroleerd; alleen drie ontbrekende operatie-intervallen waren al geïmputeerd. Geen nieuwe imputatie, transformatie, multivariabele fit of bootstrap is uitgevoerd.', '',
  'Laboratoriumwaarden zijn in zowel model_input als de opgeslagen imputaties al ln(1+x)-getransformeerd door analysis/24_model_prepare_inputs.py. Ze zijn zonder verdere transformatie gebruikt. Een eenheid op deze schaal vermenigvuldigt (1+x) met e; dit is geen OR per oorspronkelijke lab-eenheid. Leeftijd en operatie-interval blijven per jaar. Binaire vergelijkingen zijn code 1 versus 0. De overige operatiecodes bevatten ook de ongespecificeerde code 9: de referentie is daarom geen anatomisch bevestigde niet-RYGB-groep.', '',
  'Log-OR en variantie zijn over alle 50 imputaties gecombineerd met Rubin-regels via mice::pool; 95%-BI en p gebruiken de Barnard–Rubin-eindige vrijheidsgraden en de t-verdeling (complete-data df = 455 − 2 = 453). Alle 450 fits convergeerden, hadden rang 2 en eindige coëfficiënten/varianties, zonder grensoplossing. Modelwaarschuwingen en fitkenmerken staan in fit_checks.csv; inputcontroles in input_checks.csv.', '',
  'De bestaande imputaties bevatten bewust geen uitkomstinformatie (oorspronkelijke voorspelmodelstrategie). Deze keuze is behouden en vormt een beperking voor inferentiële associaties. De beschikbaarheidsselectie voor laboratoriumgegevens, onzekere volgorde binnen dezelfde dag en bestaande uitkomst-/brononzekerheden blijven gelden. Dit is geen externe validatie.', '',
  'Bron: results/analysis_2026-10-27_no_adjudication; bronspecificatie model_specification.json. Reproduceerbaar commando vanaf de PhD-map: Rscript OPERATE/analysis/29_SEH_lab_univariable.R OPERATE. Exacte numerieke resultaten en variantiecomponenten: SEH_lab_univariable.csv. Sessiedetails: sessionInfo.txt; broncontrolesommen: source_md5.csv.')
writeLines(lines, file.path(out, 'results.md'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
source_paths <- c(input_path, paths, file.path(source_dir, 'model_specification.json'),
  file.path(project, 'analysis/24_model_prepare_inputs.py'), file.path(project, 'analysis/29_SEH_lab_univariable.R'))
write_csv(data.frame(file = sub(paste0(project, '/'), '', source_paths, fixed = TRUE),
  md5 = unname(tools::md5sum(source_paths))), 'source_md5.csv')
print(results[, c('predictor', 'OR', 'CI95_lower', 'CI95_upper', 'p_value')], row.names = FALSE)
