#!/usr/local/bin/Rscript
# Figure-only rendering of stored SEH laboratory ridge odds ratios.
# Input: 2026-10-27 no-adjudication release, coefficients and percentile limits.
# Population: n=455; 79 recorded definitive gallstone events. No fitting occurs.
# Labs are log1p transformed; intervals are descriptive bootstrap intervals.
# Run from workspace root: /usr/local/bin/Rscript OPERATE/analysis/28_forestplot_SEH_lab.R

# 1. Locate the established project and validate the stored results.
args <- commandArgs(trailingOnly = FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)])
project <- normalizePath(file.path(dirname(script), '..'), mustWork = TRUE)
.libPaths(c(file.path(project, 'analysis', 'R-library'), .libPaths()))
output_dir <- file.path(project, 'results', 'analysis_2026-10-27_no_adjudication')
source_file <- file.path(output_dir, 'SEH_lab_coefficients.csv')
source_md5 <- unname(tools::md5sum(source_file))
d <- read.csv(source_file, stringsAsFactors = FALSE)
d <- d[d$term != '(Intercept)', ]
expected <- c('age', 'male', 'years_since_surgery', 'RYGB', 'RUQ_pain',
              'bilirubin', 'ALP', 'ALT', 'CRP')
stopifnot(identical(d$term, expected), nrow(d) == 9L)
values <- as.matrix(d[c('exp_coefficient', 'exp_lower', 'exp_upper')])
stopifnot(all(is.finite(values)), all(values > 0),
          all(d$exp_lower <= d$exp_coefficient),
          all(d$exp_coefficient <= d$exp_upper))
expected_values <- exp(as.matrix(d[c('MI_mean_log_coefficient',
                                     'bootstrap_lower', 'bootstrap_upper')]))
stopifnot(max(abs(values - expected_values)) < 1e-10)
labels <- c('Leeftijd (per jaar)', 'Man (versus vrouw)',
            'Tijd sinds operatie (per jaar)', 'RYGB (versus overige operatie)',
            'Pijn rechterbovenbuik (ja versus nee)',
            'Bilirubine (per 1 eenheid ln[1 + waarde])',
            'AF (per 1 eenheid ln[1 + waarde])',
            'ALAT (per 1 eenheid ln[1 + waarde])',
            'CRP (per 1 eenheid ln[1 + waarde])')
fmt <- function(x) formatC(x, format = 'f', digits = 2, decimal.mark = ',')
numeric_labels <- paste0(fmt(d$exp_coefficient), ' [', fmt(d$exp_lower),
                         '; ', fmt(d$exp_upper), ']')

# 2. Shared print-size theme: 220 x 160 mm, Helvetica, uniform blue marks.
draw <- function() {
  par(family = 'Helvetica', ps = 10, oma = c(5.5, 0.7, 4, 0.7))
  layout(matrix(1:3, nrow = 1), widths = c(0.48, 0.33, 0.19))
  y <- 9:1
  ylim <- c(0.3, 10)
  par(mar = c(3, 0.3, 0, 0.2))
  plot.new(); plot.window(xlim = c(0, 1), ylim = ylim, xaxs = 'i', yaxs = 'i')
  text(0, 9.9, 'Voorspeller', adj = c(0, 0.5), font = 2, cex = 0.95)
  text(0, y, labels, adj = c(0, 0.5), cex = 0.94)
  par(mar = c(3, 0.7, 0, 0.5))
  plot.new(); plot.window(xlim = c(0.2, 10), ylim = ylim, log = 'x',
                          xaxs = 'i', yaxs = 'i')
  abline(v = 1, lty = 2, col = 'grey45', lwd = 0.9)
  segments(d$exp_lower, y, d$exp_upper, y, col = '#0072B2', lwd = 1.4)
  segments(d$exp_lower, y - 0.10, d$exp_lower, y + 0.10, col = '#0072B2')
  segments(d$exp_upper, y - 0.10, d$exp_upper, y + 0.10, col = '#0072B2')
  points(d$exp_coefficient, y, pch = 15, col = '#0072B2', cex = 0.85)
  axis(1, at = c(0.2, 0.5, 1, 2, 5, 10),
       labels = c('0,2', '0,5', '1', '2', '5', '10'),
       tck = -0.025, mgp = c(2, 0.55, 0), cex.axis = 0.9)
  mtext('Oddsratio (logaritmische schaal)', side = 1, line = 1.8, cex = 0.86)
  par(mar = c(3, 0.2, 0, 0.1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = ylim, xaxs = 'i', yaxs = 'i')
  text(1, 9.9, 'OR [95%-BI]', adj = c(1, 0.5), font = 2, cex = 0.95)
  text(1, y, numeric_labels, adj = c(1, 0.5), cex = 0.94)
  mtext('SEH-model met laboratoriumwaarden', side = 3, outer = TRUE,
        line = 2, adj = 0, font = 2, cex = 1.2)
  mtext('Geregistreerde definitieve galsteenziekte | n = 455; 79 events',
        side = 3, outer = TRUE, line = 0.6, adj = 0, cex = 0.95)
  foot <- c('95%-bootstrap-percentielintervallen; 200 resamples; ridge-gepenaliseerd model.',
            'Laboratorium-OR: per 1 eenheid ln(1 + waarde), dus niet per oorspronkelijke meeteenheid.',
            'Oorspronkelijke eenheden: bilirubine µmol/L; AF en ALAT U/L; CRP mg/L.')
  for (i in seq_along(foot)) mtext(foot[i], side = 1, outer = TRUE,
                                  line = 1.3 + (i - 1) * 1.2, adj = 0, cex = 0.86)
}

# 3. Export editable vector and high-resolution raster versions; record checks.
prefix <- file.path(output_dir, 'SEH_lab_forestplot')
pdf(paste0(prefix, '.pdf'), width = 220/25.4, height = 160/25.4,
    family = 'Helvetica', useDingbats = FALSE)
draw(); invisible(dev.off())
png(paste0(prefix, '.png'), width = 220, height = 160, units = 'mm',
    res = 300, type = 'quartz', bg = 'white')
draw(); invisible(dev.off())
# SVG omitted: system Cairo libraries unavailable; PDF is the vector master.
stopifnot(all(file.exists(paste0(prefix, c('.pdf', '.png')))))
stopifnot(identical(source_md5, unname(tools::md5sum(source_file))))
writeLines(c(paste('Source:', source_file), paste('Source MD5:', source_md5),
             'Validated: nine original-order non-intercept terms; finite positive ORs and limits;',
             'ordered CI bounds; numeric agreement with exponentiated stored log coefficients <1e-10.',
             'Source checksum unchanged after figure rendering. No models refitted.',
             'PDF: vector; PNG: 300 dpi; dimensions: 220 x 160 mm.',
             'SVG omitted because system Cairo libraries are unavailable.',
             capture.output(sessionInfo())), paste0(prefix, '_render_log.txt'))
message('Rendered forest plot from unchanged stored results: ', prefix)
