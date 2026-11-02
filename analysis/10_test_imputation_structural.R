#!/usr/bin/env Rscript
# Synthetic regression: one imputed cell per outcome stratum, structural NA,
# and partly missing auxiliary must not block eligible cells.
.libPaths(c(file.path(getwd(),'OPERATE','analysis','R-library'),.libPaths()))
set.seed(20261028)
p<-tempfile('10_operate_structural_');dir.create(p)
jsonlite::write_json(list(outcome_mapping_status='validated'),file.path(p,'configuration.json'),auto_unbox=TRUE)
n<-100L
x<-data.frame(participant_id=sprintf('S%03d',1:n),age_years_at_index=runif(n,20,70),years_since_bariatric_surgery=runif(n),gallstone_related_disease=c(rep(0:1,25),rep(NA,50)),value=rnorm(n))
x$age_alias<-x$age_years_at_index
x$complete_aux<-rnorm(n)
x$years_since_bariatric_surgery[seq(1,n,2)]<-NA
x$value[c(10,60,100)]<-NA
w<-data.frame(participant_id=x$participant_id,value=seq_len(n) %in% c(10,60))
e<-data.frame(participant_id=x$participant_id,value=seq_len(n)!=100)
sp<-data.frame(variable=c('value','age_alias','complete_aux'),method=c('pmm','',''),block=c('clinical','clinical','other_form'),levels_json='',reason='synthetic',auxiliary_for_blocks=c('[]','[]','["clinical"]'),auxiliary_rank=c(1,2,3))
wr<-function(z,f){con<-gzfile(f,'wt');write.csv(z,con,row.names=FALSE,na='');close(con)}
for(st in c('SEH','POLI')) {wr(x,file.path(p,paste0('imputation_input_',st,'.csv.gz')));wr(w,file.path(p,paste0('where_',st,'.csv.gz')));wr(e,file.path(p,paste0('eligibility_',st,'.csv.gz')));write.csv(sp,file.path(p,paste0('imputation_spec_',st,'.csv')),row.names=FALSE,na='')}
status<-system2(file.path(R.home('bin'),'Rscript'),c('OPERATE/analysis/10_impute_first_visit_v2.R',shQuote(p)))
stopifnot(status==0L)
for(st in c('SEH','POLI')) {
 inv<-read.csv(file.path(p,paste0('imputation_inventory_',st,'.csv')))
 stopifnot(nrow(inv)==2L,all(inv$status=='imputed'))
 cv<-read.csv(file.path(p,paste0('imputation_convergence_',st,'.csv')))
 stopifnot(all(!cv$sd_applicable),all(is.na(cv$psrf_sd)),all(is.finite(cv$psrf)))
 fit<-readRDS(file.path(p,paste0('mice_',st,'.rds')))
 stopifnot(length(fit$fits)==2L)
 audit<-read.csv(file.path(p,paste0('imputation_predictor_audit_',st,'.csv')))
 stopifnot(all(audit$status[audit$predictor=='age_alias']=='redundant_or_rank_deficient_predictor'),all(audit$status[audit$predictor=='complete_aux']=='selected'))
 for(f in fit$fits) stopifnot(!'years_since_bariatric_surgery' %in% names(f$data))
 for(i in 1:50) {
  cc<-read.csv(gzfile(file.path(p,sprintf('imputed_%s_%02d.csv.gz',st,i))),na.strings='')
  stopifnot(identical(cc$participant_id,x$participant_id),identical(cc$gallstone_related_disease,x$gallstone_related_disease),is.na(cc$value[100]),all(is.finite(cc$value[c(10,60)])),identical(is.na(cc$years_since_bariatric_surgery),is.na(x$years_since_bariatric_surgery)))
  for(v in names(x)) {obs<-!is.na(x[[v]]);stopifnot(isTRUE(all.equal(cc[[v]][obs],x[[v]][obs],tolerance=1e-12)))}
 }
}
cat('Single-cell mean-only convergence, structural preservation, partial auxiliary exclusion, and outcome-unknown stratum passed for both settings.\n')
