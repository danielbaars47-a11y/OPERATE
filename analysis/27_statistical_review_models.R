#!/usr/bin/env Rscript
# Independent saved-output checks; never refits a clinical model or reruns study MI.
# A small synthetic PMM donor/ignore test extracts only the impute function AST.
.libPaths(c('OPERATE/analysis/R-library',.libPaths()))
stopifnot(requireNamespace('glmnet',quietly=TRUE),requireNamespace('mice',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE))
root <- 'OPERATE/results/analysis_2026-10-27_no_adjudication'
config <- jsonlite::fromJSON(file.path(root,'model_specification.json'),simplifyVector=FALSE)
checks <- list(); summaries <- list(); pending <- character()
ck <- function(name,x) checks[[length(checks)+1L]] <<- list(check=name,passed=isTRUE(x))
near <- function(x,y,tol=1e-8) length(x)==length(y)&&all(is.finite(x))&&all(is.finite(y))&&max(abs(x-y))<tol
read <- function(name) read.csv(file.path(root,name),check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(participant_id='character'))
independent_metrics <- function(y,p) {
  n1<-sum(y==1);n0<-sum(y==0);auc<-(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)
  lp<-qlogis(pmin(1-1e-8,pmax(1e-8,p)))
  intercept<-glm(y~1,offset=lp,family=binomial())
  slope<-glm(y~lp,family=binomial())
  stopifnot(intercept$converged,slope$converged)
  c(AUC=auc,Brier=mean((y-p)^2),calibration_intercept=unname(coef(intercept)[1]),calibration_slope=unname(coef(slope)[2]))
}
for (index in seq_along(config$models)) {
  model<-config$models[[index]];id<-model$id
  if(!file.exists(file.path(root,paste0(id,'_summary.json')))){pending<-c(pending,id);next}
  d<-read(paste0('model_input_',id,'.csv'));p<-read(paste0(id,'_predictions.csv'));perf<-read(paste0(id,'_performance.csv'));status<-read(paste0(id,'_bootstrap_status.csv'));coeff<-read(paste0(id,'_coefficients.csv'))
  z<-readRDS(file.path(root,paste0(id,'_final_models.rds')));boot<-readRDS(file.path(root,paste0(id,'_bootstrap.rds')))
  y<-d$outcome;x<-d[setdiff(names(d),c('participant_id','outcome'))];n<-nrow(d)
  ck(paste0(id,'_population'),n==model$expected_n&&sum(y)==model$expected_events&&!anyDuplicated(d$participant_id))
  ck(paste0(id,'_prediction_identity'),identical(d$participant_id,p$participant_id)&&identical(y,p$outcome))
  ck(paste0(id,'_fifty_saved_models'),length(z$fitted_models)==50&&length(z$imputation$data)==50)
  ck(paste0(id,'_predictor_contract'),identical(names(x),vapply(model$predictors,'[[',character(1),'name')))
  ck(paste0(id,'_no_outcome_in_imputation'),!any(grepl('outcome|gallstone|participant',names(z$imputation$fit$data),ignore.case=TRUE)))
  ck(paste0(id,'_final_imputer_no_evaluation_rows'),!any(z$imputation$fit$ignore))
  ck(paste0(id,'_final_convergence'),!isTRUE(z$imputation$convergence_concern))
  pp<-matrix(NA_real_,n,50)
  for(i in seq_len(50)) {
    dd<-z$imputation$data[[i]];fit<-z$fitted_models[[i]]
    ck(paste0(id,'_draw',i,'_ridge_penalty'),near(fit$lambda,1/n,1e-12))
    ck(paste0(id,'_draw',i,'_finite_complete'),all(is.finite(as.matrix(dd))))
    protected<-!is.na(as.matrix(x));ck(paste0(id,'_draw',i,'_observed_predictors_preserved'),near(as.matrix(dd)[protected],as.matrix(x)[protected]))
    ck(paste0(id,'_draw',i,'_PMM_training_donors'),all(dd$years_since_surgery[is.na(x$years_since_surgery)] %in% x$years_since_surgery[!is.na(x$years_since_surgery)]))
    pp[,i]<-as.numeric(predict(fit,as.matrix(dd),s=1/n,type='response'))
  }
  ck(paste0(id,'_saved_predictions_reproduced_without_refit'),near(rowMeans(pp),p$prediction_apparent))
  ap<-independent_metrics(y,p$prediction_apparent);oo<-independent_metrics(y,p$prediction_OOB)
  ck(paste0(id,'_apparent_metrics'),near(ap,perf$apparent))
  ck(paste0(id,'_OOB_metrics'),near(oo,perf$OOB_averaged_prediction_metric))
  valid<-vapply(boot,function(a)isTRUE(a$ok),logical(1));success<-boot[valid]
  ck(paste0(id,'_bootstrap_denominators'),length(boot)==200&&sum(valid)>=180&&sum(status$status=='completed')==sum(valid))
  ck(paste0(id,'_no_accepted_convergence_concerns'),!any(vapply(success,function(a)isTRUE(a$convergence_concern),logical(1))))
  train<-do.call(rbind,lapply(success,'[[','train'));test<-do.call(rbind,lapply(success,'[[','test'))
  ck(paste0(id,'_finite_bootstrap_metrics'),all(is.finite(train))&&all(is.finite(test)))
  optimism<-colMeans(train-test)
  ck(paste0(id,'_optimism_arithmetic'),near(optimism,perf$mean_optimism)&&near(ap-optimism,perf$optimism_corrected))
  bc<-do.call(cbind,lapply(success,'[[','coef'));limits<-t(apply(bc,1,quantile,c(.025,.975)))
  ck(paste0(id,'_coefficient_bootstrap_intervals'),near(rowMeans(z$coefficient_draws),coeff$MI_mean_log_coefficient)&&near(limits[,1],coeff$bootstrap_lower)&&near(limits[,2],coeff$bootstrap_upper))
  seed<-if(!is.null(model$population_seed))as.integer(model$population_seed) else if(!is.null(model$imputation_seed))as.integer(model$imputation_seed) else 27102026L+index*100000L
  set.seed(seed);samples<-replicate(200,sample.int(n,n,replace=TRUE),simplify=FALSE)
  for(b in which(valid))ck(paste0(id,'_bootstrap',b,'_true_OOB_membership'),identical(is.na(boot[[b]]$oob_prediction),seq_len(n)%in%unique(samples[[b]])))
  om<-do.call(cbind,lapply(success,'[[','oob_prediction'));on<-rowSums(!is.na(om));op<-rowMeans(om,na.rm=TRUE)
  ck(paste0(id,'_OOB_ensemble_reproduced'),near(op,p$prediction_OOB)&&identical(as.integer(on),as.integer(p$OOB_replicates))&&all(on>0))
  summaries[[id]]<-list(n=n,events=sum(y),slopes=ncol(x),valid_bootstrap=sum(valid),failed_bootstrap=sum(!valid),OOB_min=min(on),OOB_max=max(on),corrected_metrics=as.list(ap-optimism))
}
# Independently verify paired incremental comparisons on common valid resamples.
pair_path<-file.path(root,'same_population_model_comparisons.csv')
if(file.exists(pair_path)) {
  comparisons<-read.csv(pair_path,stringsAsFactors=FALSE)
  pairs<-unique(comparisons[c('base','extended')])
  for(j in seq_len(nrow(pairs))) {
    base<-pairs$base[j];extended<-pairs$extended[j];label<-paste0(base,'_vs_',extended)
    a<-read(paste0('model_input_',base,'.csv'));b<-read(paste0('model_input_',extended,'.csv'))
    ck(paste0(label,'_identical_population_outcomes'),identical(a$participant_id,b$participant_id)&&identical(a$outcome,b$outcome))
    ba<-readRDS(file.path(root,paste0(base,'_bootstrap.rds')));bb<-readRDS(file.path(root,paste0(extended,'_bootstrap.rds')))
    common<-which(vapply(ba,function(z)isTRUE(z$ok),logical(1))&vapply(bb,function(z)isTRUE(z$ok),logical(1)))
    ta<-do.call(rbind,lapply(ba[common],'[[','train'));ea<-do.call(rbind,lapply(ba[common],'[[','test'));tb<-do.call(rbind,lapply(bb[common],'[[','train'));eb<-do.call(rbind,lapply(bb[common],'[[','test'))
    pa<-read(paste0(base,'_performance.csv'));pb<-read(paste0(extended,'_performance.csv'))
    ad<-pb$apparent-pa$apparent;od<-(tb-eb)-(ta-ea);ed<-eb-ea;given<-comparisons[comparisons$base==base&comparisons$extended==extended,]
    ck(paste0(label,'_paired_denominator'),all(given$paired_valid_bootstraps==length(common)))
    ck(paste0(label,'_incremental_optimism_arithmetic'),near(ad,given$apparent_difference)&&near(colMeans(od),given$mean_paired_optimism_difference)&&near(ad-colMeans(od),given$paired_optimism_corrected_difference))
    ck(paste0(label,'_raw_distribution_quantiles'),near(apply(ed,2,quantile,.025),given$raw_original_evaluation_difference_p025)&&near(apply(ed,2,quantile,.975),given$raw_original_evaluation_difference_p975)&&near(apply(od,2,quantile,.025),given$paired_optimism_difference_p025)&&near(apply(od,2,quantile,.975),given$paired_optimism_difference_p975))
  }
} else pending<-c(pending,'paired_comparisons')
# Only extract function assignments, avoiding the study execution at script end.
expressions<-parse('OPERATE/analysis/24_model_pipeline.R');env<-new.env(parent=globalenv())
for(ex in expressions)if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&as.character(ex[[2]])%in%c('trace_check','impute'))eval(ex,env)
set.seed(5001);synthetic<-data.frame(age=runif(90,20,70),male=sample(0:1,90,TRUE),years_since_surgery=runif(90,1,12),RUQ=sample(0:1,90,TRUE));synthetic$years_since_surgery[c(2,8,72,76)]<-NA_real_;synthetic$years_since_surgery[80:90]<-10000
syn<-env$impute(synthetic,70,5,5002,10,10)
ck('synthetic_ignore_excludes_evaluation_fit_rows',identical(as.logical(syn$fit$ignore),seq_len(90)>70))
donors<-synthetic$years_since_surgery[seq_len(70)];donors<-donors[is.finite(donors)]
ck('synthetic_evaluation_values_never_PMM_donors',all(vapply(syn$data,function(dd)all(dd$years_since_surgery[is.na(synthetic$years_since_surgery)]%in%donors),logical(1))))
rejected<-tryCatch({env$impute(transform(synthetic,outcome=sample(0:1,90,TRUE)),70,5,5002,10,10);FALSE},error=function(e)TRUE)
ck('synthetic_outcome_column_rejected',rejected)
# Hash provenance without exposing patient values.
paths<-c('model_specification.json','frozen_analysis_cohort.csv','frozen_model_populations.csv',paste0('model_input_',vapply(config$models,'[[',character(1),'id'),'.csv'))
paths<-paths[file.exists(file.path(root,paths))]
provenance<-as.list(tools::md5sum(file.path(root,paths)));names(provenance)<-paths
failed<-Filter(function(a)!a$passed,checks)
result<-list(review='Independent post-execution saved-fit/output review plus synthetic leakage guards',complete=!length(pending),passed=!length(failed),pending=pending,checks_n=length(checks),failed_checks=failed,models=summaries,input_md5=provenance,limitations=c('Internal validation does not resolve ascertainment, timing or transportability limitations.','OOB-averaged ensemble metrics differ from optimism-corrected full-model metrics.','Bootstrap optimism percentiles are not confidence intervals for model performance.'))
jsonlite::write_json(result,file.path(root,'post_execution_statistical_review.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Complete:',result$complete,' Passed:',result$passed,' Checks:',length(checks),' Pending:',paste(pending,collapse=','),'\n')
if(length(failed)){print(failed);quit(status=1)}
