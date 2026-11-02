#!/usr/bin/env Rscript
# Recover/check recorded-seed bootstrap indices and export saved diagnostics.
# No model estimation, imputation, or performance recomputation.
.libPaths(c('OPERATE/analysis/R-library',.libPaths()))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))args[1] else 'OPERATE/results/analysis_2026-10-27_no_adjudication'
spec<-jsonlite::fromJSON(file.path(root,'model_specification.json'),simplifyVector=FALSE)
allcv<-list()
for(i in seq_along(spec$models)){
 model<-spec$models[[i]];id<-model$id
 d<-read.csv(file.path(root,paste0('model_input_',id,'.csv')),colClasses=c(participant_id='character'));n<-nrow(d)
 seed<-if(!is.null(model$imputation_seed))as.integer(model$imputation_seed)else 27102026L+i*100000L
 if(!is.null(model$population_seed))seed<-as.integer(model$population_seed)
 set.seed(seed);indices<-replicate(200,sample.int(n,n,replace=TRUE),simplify=FALSE)
 boots<-readRDS(file.path(root,paste0(id,'_bootstrap.rds')))
 stopifnot(length(boots)==length(indices))
 for(b in seq_along(boots))if(isTRUE(boots[[b]]$ok))stopifnot(identical(is.na(boots[[b]]$oob_prediction),seq_len(n)%in%unique(indices[[b]])))
 saveRDS(indices,file.path(root,paste0(id,'_bootstrap_indices.rds')))
 final<-readRDS(file.path(root,paste0(id,'_final_models.rds')))
 cv<-final$imputation$convergence
 if(nrow(cv)){cv$model<-id;cv$iterations<-final$imputation$iterations;allcv[[id]]<-cv}
 capture.output(final$sessionInfo,file=file.path(root,paste0(id,'_run_environment.txt')))
}
write.csv(do.call(rbind,allcv),file.path(root,'model_imputation_convergence.csv'),row.names=FALSE)
cat('Bootstrap seed reconstruction matches all saved OOB masks; saved indices and final MI diagnostics exported.\n')
