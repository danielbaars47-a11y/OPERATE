#!/usr/bin/env Rscript
# Independent read-only checks of saved MICE objects; never calls mice().
.libPaths(c('OPERATE/analysis/R-library', .libPaths()))
root <- 'OPERATE/results/first_visit_2026-10-27_v3'
result <- list()
for (s in c('SEH','POLI')) {
  path <- file.path(root,paste0('mice_',s,'.rds'))
  z <- readRDS(path)
  input <- file.path(root,paste0('imputation_input_',s,'.csv.gz'))
  md5 <- identical(unname(z$input_md5),unname(tools::md5sum(input)))
  ms <- vapply(z$fits,function(f) as.integer(f$m),integer(1))
  iterations <- vapply(z$fits,function(f) as.integer(f$iteration),integer(1))
  no_outcome <- all(vapply(z$fits,function(f) !'gallstone_related_disease' %in% names(f$method)[nzchar(f$method)],logical(1)))
  stopifnot(md5,length(ms)>0L,all(ms==50L),all(iterations %in% c(30L,100L)),no_outcome)
  result[[s]] <- list(saved_fits=length(z$fits),all_m_50=all(ms==50L),iterations=as.list(table(iterations)),input_md5_match=md5,outcome_not_imputed=no_outcome)
}
jsonlite::write_json(list(review='Independent saved MICE provenance and configuration checks',passed=TRUE,settings=result,R_version=R.version.string),file.path(root,'rds_validation.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Independent saved RDS checks passed for both settings\n')
