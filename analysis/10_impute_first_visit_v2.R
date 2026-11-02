#!/usr/bin/env Rscript
# OPERATE first-visit MI, 2026-10-27. Separate settings, no outcome imputation
# or prediction model fitting. Inputs are frozen canonical data/specification,
# cell-level where and eligibility matrices. Raw exports are never modified.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork=TRUE) else normalizePath('OPERATE/results/first_visit_2026-10-27', mustWork=TRUE)
.libPaths(c(file.path(dirname(dirname(root)), 'analysis', 'R-library'), .libPaths()))
.libPaths(c(Sys.getenv('OPERATE_R_LIBRARY', unset=file.path(getwd(),'OPERATE','analysis','R-library')), .libPaths()))
stopifnot(requireNamespace('mice', quietly=TRUE), requireNamespace('jsonlite', quietly=TRUE))
config_path <- file.path(root,'configuration.json')
if(!file.exists(config_path)) stop('Validated frozen outcome configuration.json is required')
config <- jsonlite::fromJSON(config_path)
if(!identical(config$outcome_mapping_status,'validated')) stop('Outcome mapping is provisional: clinical MI execution blocked')
read_csv <- function(p) read.csv(gzfile(p), check.names=FALSE, stringsAsFactors=FALSE, colClasses=c(participant_id='character'), na.strings='')
write_csv <- function(x,p) { con <- gzfile(p,'wt'); on.exit(close(con)); write.csv(x,con,row.names=FALSE,na='') }
as_bool <- function(x) { y <- tolower(as.character(x)); if(any(!is.na(y) & !y %in% c('true','false','1','0'))) stop('Invalid boolean matrix'); !is.na(y) & y %in% c('true','1') }
run_setting <- function(setting) {
  input_path <- file.path(root,paste0('imputation_input_',setting,'.csv.gz'))
  dat <- read_csv(input_path)
  spec <- read.csv(file.path(root,paste0('imputation_spec_',setting,'.csv')),check.names=FALSE,stringsAsFactors=FALSE,na.strings='')
  where_raw <- read_csv(file.path(root,paste0('where_',setting,'.csv.gz')))
  eligible_raw <- read_csv(file.path(root,paste0('eligibility_',setting,'.csv.gz')))
  stopifnot(!anyDuplicated(dat$participant_id), all(c('variable','method','block','levels_json','reason') %in% names(spec)))
  for (mat in list(where_raw,eligible_raw)) {
    stopifnot(nrow(mat)==nrow(dat))
    if ('participant_id' %in% names(mat)) stopifnot(identical(as.character(mat$participant_id),as.character(dat$participant_id)))
  }
  targets <- spec$variable[!is.na(spec$method) & nzchar(spec$method)]
  stopifnot(!anyDuplicated(spec$variable),all(targets %in% names(dat)), all(targets %in% names(where_raw)),all(targets %in% names(eligible_raw)))
  stopifnot(all(spec$method[match(targets,spec$variable)] %in% c('pmm','logreg','polyreg','polr')))
  stopifnot(!any(targets %in% c('participant_id','gallstone_related_disease','setting')))
  where <- matrix(FALSE,nrow(dat),ncol(dat),dimnames=list(NULL,names(dat)))
  eligibility <- where
  for(v in targets) {
    where[,v] <- as_bool(where_raw[[v]])
    eligibility[,v] <- as_bool(eligible_raw[[v]])
    stopifnot(all(!where[,v] | (is.na(dat[[v]]) & eligibility[,v])))
    s <- spec[match(v,spec$variable),]
    # read.csv infers an entirely empty numeric column as logical. Preserve
    # its empty values while restoring the declared numeric type.
    if(s$method=='pmm' && is.logical(dat[[v]]) && all(is.na(dat[[v]]))) dat[[v]]<-as.numeric(dat[[v]])
    if(s$method=='pmm' && (!is.numeric(dat[[v]]) || any(!is.finite(dat[[v]][!is.na(dat[[v]])])))) stop('PMM requires finite numeric observed values: ',v)
    if(s$method!='pmm' && (!'verification_state' %in% names(s) || is.na(s$verification_state) || !s$verification_state %in% c('historical_mapping_current_question_and_options_concordant','current_machine_dictionary_verified'))) stop('Categorical mapping not verified: ',v)
    if (!is.na(s$levels_json) && nzchar(s$levels_json)) {
      lev <- as.character(jsonlite::fromJSON(s$levels_json))
      stopifnot(all(is.na(dat[[v]]) | as.character(dat[[v]]) %in% lev))
      dat[[v]] <- factor(as.character(dat[[v]]),levels=lev,ordered=identical(s$method,'polr'))
    } else if (s$method!='pmm') stop('Categorical target lacks verified levels: ',v)
  }
  original <- dat
  completed <- replicate(50L,dat,simplify=FALSE)
  fits <- list(); inventory <- list(); diagnostics <- list(); predictor_audit <- list()
  # Targets share a chain only when their clinical block and eligibility agree.
  keys <- vapply(targets,function(v) paste(spec$block[match(v,spec$variable)],paste(which(eligibility[,v]),collapse=','),sep=':'),character(1))
  groups0 <- split(targets,keys)
  groups <- list(); group_rows <- list()
  for(vs0 in groups0) {
    rr <- which(eligibility[,vs0[[1]]])
    strata <- if('gallstone_related_disease' %in% names(dat)) split(rr,is.na(dat$gallstone_related_disease[rr])) else list(rr)
    for(rr0 in strata) if(length(rr0)) { groups[[length(groups)+1L]]<-vs0; group_rows[[length(group_rows)+1L]]<-rr0 }
  }
  for(g in seq_along(groups)) {
    vs <- groups[[g]]; rows <- group_rows[[g]]; active <- vs[colSums(where[rows,vs,drop=FALSE])>0]
    if(!length(active)) next
    # Prespecified priorities: fixed core, then CRF question order within form.
    # Candidate availability/rank/support checks never use outcome association.
    core <- intersect(c('age_years_at_index','geslacht','years_since_bariatric_surgery','bariatric_surgery_type','gallstone_related_disease'),names(dat))
    block_label <- spec$block[match(vs[[1]],spec$variable)]
    ranked <- spec$variable
    rr <- if('auxiliary_rank' %in% names(spec)) suppressWarnings(as.numeric(spec$auxiliary_rank)) else seq_len(nrow(spec))
    rr[!is.finite(rr)] <- Inf
    ranked <- ranked[order(rr,seq_len(nrow(spec)))]
    declared <- spec$variable[!is.na(spec$block) & spec$block==block_label]
    if('auxiliary_for_blocks' %in% names(spec)) for(j in seq_len(nrow(spec))) {
      st<-spec$auxiliary_for_blocks[j]
      if(!is.na(st) && nzchar(st) && block_label %in% jsonlite::fromJSON(st)) declared<-union(declared,spec$variable[j])
    }
    candidates <- unique(c(core,ranked[ranked %in% declared]))
    candidates <- intersect(candidates,setdiff(names(dat),c('participant_id','setting')))
    candidates <- candidates[!grepl('(^eligible|_eligible|_flag$|date|free.?text|comment)',candidates,ignore.case=TRUE)]
    # Apply verified categorical metadata also to complete auxiliary fields.
    for(v in setdiff(candidates,targets)) {
      j<-match(v,spec$variable)
      if(!is.na(j) && !is.na(spec$levels_json[j]) && nzchar(spec$levels_json[j])) {
        if(!'verification_state' %in% names(spec) || is.na(spec$verification_state[j]) || !spec$verification_state[j] %in% c('historical_mapping_current_question_and_options_concordant','current_machine_dictionary_verified')) next
        lev<-as.character(jsonlite::fromJSON(spec$levels_json[j]))
        if(any(!is.na(dat[[v]]) & !as.character(dat[[v]]) %in% lev)) stop('Auxiliary outside verified levels: ',v)
        dat[[v]]<-factor(as.character(dat[[v]]),levels=lev)
      }
    }
    usable <- active[vapply(active,function(v) { z<-dat[rows,v]; obs<-z[!is.na(z)]; if(spec$method[match(v,spec$variable)]=='pmm') length(obs)>=20L && length(unique(obs))>=5L else length(obs)>0L && all(table(z)>=5L) },logical(1))]
    rejected <- setdiff(active,usable)
    for(v in rejected) inventory[[length(inventory)+1L]] <- data.frame(variable=v,block=g,status='unreliable_not_imputed',reason='Insufficient observed support: continuous requires 20 observations/5 distinct; categorical requires 5 in every declared level',stringsAsFactors=FALSE)
    if(!length(usable)) next
    aux <- setdiff(candidates,usable)
    aux <- aux[vapply(dat[rows,aux,drop=FALSE],function(z) !anyNA(z) && (is.numeric(z)||is.factor(z)||is.logical(z)) && (!is.numeric(z)||all(is.finite(z))),logical(1))]
    cols <- unique(c(usable,aux)); x <- dat[rows,cols,drop=FALSE]; w <- where[rows,cols,drop=FALSE]
    meth <- setNames(rep('',length(cols)),cols); meth[usable] <- spec$method[match(usable,spec$variable)]
    pred <- matrix(0L,length(cols),length(cols),dimnames=list(cols,cols))
    for(v in usable) {
      chosen<-character(); observed<-!is.na(x[[v]])
      ordered <- unique(c(core,ranked[ranked %in% c(usable,candidates)]))
      ordered <- setdiff(ordered,v)
      for(pp in ordered) {
        reason<-''
        if(!pp %in% cols) reason<-'incomplete_unimputed_or_untyped_auxiliary'
        else if(any(is.na(x[[pp]]) & !w[,pp])) reason<-'structural_or_nonimputable_missing_predictor'
        else if(length(unique(x[[pp]][!is.na(x[[pp]])]))<2L) reason<-'constant_predictor'
        if(!nzchar(reason)) {
          dd<-x[observed,c(chosen,pp),drop=FALSE];dd<-dd[complete.cases(dd),,drop=FALSE]
          mm<-tryCatch(model.matrix(~.,dd),error=function(e) NULL)
          if(is.null(mm)||!nrow(mm)) reason<-'unavailable_observed_design'
          else if(ncol(mm)-1L>10L) reason<-'maximum_10_predictor_design_df'
          else if(sum(observed)<=10L*ncol(mm)) reason<-'observed_target_not_greater_than_10_times_design_columns'
          else if(qr(mm)$rank<ncol(mm)) reason<-'redundant_or_rank_deficient_predictor'
        }
        if(!nzchar(reason)) { chosen<-c(chosen,pp);pred[v,pp]<-1L;reason<-'selected' }
        predictor_audit[[length(predictor_audit)+1L]]<-data.frame(block=g,clinical_block=block_label,target=v,predictor=pp,status=reason,n_observed=sum(observed),stringsAsFactors=FALSE)
      }
    }
    # A target without any supported predictor is retained missing, never
    # filled by an intercept-only model. Other targets remain in the joint chain.
    repeat {
    unsupported<-usable[rowSums(pred[usable,,drop=FALSE])==0L]
    if(length(unsupported)) {
      for(v in unsupported) inventory[[length(inventory)+1L]]<-data.frame(variable=v,block=g,status='unreliable_not_imputed',reason='No supported informative predictor after prespecified technical guards')
      meth[unsupported]<-''; pred[,unsupported]<-0L;usable<-setdiff(usable,unsupported)
    } else break
    }
    if(!length(usable)) next
    seed <- 27102026L + match(setting,c('SEH','POLI'))*10000L + g
    fit <- tryCatch(mice::mice(x,m=50L,maxit=30L,method=meth,predictorMatrix=pred,where=w,seed=seed,printFlag=FALSE),error=function(e)e)
    if(inherits(fit,'error')) {
      for(v in usable) inventory[[length(inventory)+1L]] <- data.frame(variable=v,block=g,status='unreliable_not_imputed',reason=conditionMessage(fit))
      next
    }
    assess <- function(f) {
      # Split-chain potential scale reduction on chain means and SDs.
      # Uses saved MICE traces directly; avoids optional rstan installation.
      rhat <- function(z) {
        h<-floor(nrow(z)/2L); if(h<2L || any(!is.finite(z))) return(NA_real_)
        zz<-cbind(z[seq_len(h),,drop=FALSE],z[(nrow(z)-h+1L):nrow(z),,drop=FALSE])
        W<-mean(apply(zz,2,var)); B<-h*var(colMeans(zz))
        if(W==0) return(if(B==0) 1 else Inf)
        sqrt(((h-1)/h*W+B/h)/W)
      }
      cv<-do.call(rbind,lapply(usable,function(v) {
        z<-f$chainMean[v,,,drop=TRUE]; sdtrace<-sqrt(f$chainVar[v,,,drop=TRUE])
        # A single imputed cell has no within-iteration SD; only its mean
        # trace can be assessed. SD convergence is explicitly not applicable.
        sd_applicable<-sum(w[,v])>=2L
        data.frame(.it=f$iteration,vrb=v,psrf=rhat(z),psrf_sd=if(sd_applicable) rhat(sdtrace) else NA_real_,sd_applicable=sd_applicable)
      }))
      list(table=cv,poor=any(!is.finite(cv$psrf)|cv$psrf>1.1|(cv$sd_applicable & (!is.finite(cv$psrf_sd)|cv$psrf_sd>1.1))))
    }
    safe_assess <- function(f) tryCatch(assess(f),error=function(e) list(table=data.frame(),poor=TRUE,error=conditionMessage(e)))
    a <- safe_assess(fit); continuation_error <- NULL
    if(a$poor) {
      continued <- tryCatch(mice::mice.mids(fit,maxit=70L,printFlag=FALSE),error=function(e)e)
      if(inherits(continued,'error')) continuation_error <- conditionMessage(continued) else { fit<-continued; a<-safe_assess(fit) }
    }
    fits[[paste0('block_',g)]] <- fit
    diagnostics[[paste0('block_',g)]] <- list(convergence=a$table,poor=a$poor,continuation_error=continuation_error,assessment_error=a$error,loggedEvents=fit$loggedEvents,rows=rows,targets=usable,predictors=pred,seed=seed)
    unreliable <- usable[is.na(fit$method[usable]) | !nzchar(fit$method[usable])]
    if(!is.null(continuation_error)) unreliable<-usable
    draws <- tryCatch(lapply(seq_len(50L),function(i) mice::complete(fit,i)),error=function(e) NULL)
    if(is.null(draws)) unreliable<-usable else for(v in usable) {
      valid<-all(vapply(draws,function(cc) { z<-cc[[v]][w[,v]]; if(is.numeric(z)) all(is.finite(z)) else all(!is.na(z) & as.character(z) %in% levels(dat[[v]])) },logical(1)))
      if(!valid) unreliable<-union(unreliable,v)
    }
    if(a$poor && nrow(a$table)) { last<-a$table[a$table$.it==max(a$table$.it),,drop=FALSE]; unreliable<-union(unreliable,intersect(usable,unique(last$vrb[!is.finite(last$psrf)|last$psrf>1.1|(last$sd_applicable & (!is.finite(last$psrf_sd)|last$psrf_sd>1.1))]))) }
    if(a$poor && !nrow(a$table)) unreliable<-usable
    # Flag chains still unstable after 100; retain their MI draws for inspection,
    # but omit those draws from the released completed data.
    for(v in usable) inventory[[length(inventory)+1L]] <- data.frame(variable=v,block=g,status=if(v %in% unreliable)'unreliable_not_imputed' else 'imputed',reason=if(v %in% unreliable)'Postfit validity, method, continuation, or convergence check failed' else '',stringsAsFactors=FALSE)
    for(i in seq_len(50L)) {
      if(is.null(draws)) next
      cc <- draws[[i]]
      for(v in setdiff(usable,unreliable)) { idx <- rows[w[,v]]; completed[[i]][idx,v] <- cc[w[,v],v] }
    }
  }
  if(length(predictor_audit)) write.csv(do.call(rbind,predictor_audit),file.path(root,paste0('imputation_predictor_audit_',setting,'.csv')),row.names=FALSE,na='')
  if(length(diagnostics)) {
    dd<-do.call(rbind,lapply(names(diagnostics),function(nm) { z<-diagnostics[[nm]]$convergence; if(!nrow(z)) return(NULL); z$block<-nm; z$n_rows<-length(diagnostics[[nm]]$rows); z })); write.csv(dd,file.path(root,paste0('imputation_convergence_',setting,'.csv')),row.names=FALSE,na='')
  }
  saveRDS(list(setting=setting,fits=fits,diagnostics=diagnostics,spec=spec,input_md5=tools::md5sum(input_path),sessionInfo=sessionInfo()),file.path(root,paste0('mice_',setting,'.rds')))
  inv <- if(length(inventory)) do.call(rbind,inventory) else data.frame(variable=character(),block=integer(),status=character(),reason=character())
  event_tables<-lapply(names(diagnostics),function(nm) { ev<-diagnostics[[nm]]$loggedEvents; if(is.null(ev)||!nrow(ev)) return(NULL); ev$block<-nm;ev })
  events<-do.call(rbind,event_tables)
  if(!is.null(events)) write.csv(events,file.path(root,paste0('imputation_logged_events_',setting,'.csv')),row.names=FALSE,na='')
  write.csv(inv,file.path(root,paste0('imputation_inventory_',setting,'.csv')),row.names=FALSE,na='')
  for(i in seq_len(50L)) {
    cc<-completed[[i]]; changed<-where & FALSE
    for(v in names(dat)) {
      same <- (is.na(original[[v]]) & is.na(cc[[v]])) | (!is.na(original[[v]]) & !is.na(cc[[v]]) & as.character(original[[v]])==as.character(cc[[v]]))
      stopifnot(all(same[!where[,v]]))
      changed[,v] <- !same
      if(is.numeric(original[[v]]) && any(changed[,v])) stopifnot(all(cc[[v]][changed[,v]] %in% original[[v]][!is.na(original[[v]])]))
    }
    cc$.imp <- i
    write_csv(cc,file.path(root,sprintf('imputed_%s_%02d.csv.gz',setting,i)))
    ch <- data.frame(participant_id=dat$participant_id,changed[,setdiff(colnames(changed),'participant_id'),drop=FALSE],check.names=FALSE); ch$.imp<-i
    write_csv(ch,file.path(root,sprintf('imputed_changed_%s_%02d.csv.gz',setting,i)))
  }
  long <- do.call(rbind,lapply(seq_len(50L),function(i) {x<-completed[[i]];x$.imp<-i;x}))
  write_csv(long,file.path(root,paste0('imputed_',setting,'_long.csv.gz')))
  message(setting,': completed; ',length(unique(inv$variable[inv$status=='imputed'])),' unique targets released in ',sum(inv$status=='imputed'),' strata; ',length(unique(inv$variable[inv$status!='imputed'])),' unique targets unreliable in ',sum(inv$status!='imputed'),' strata.')
}
for(setting in c('SEH','POLI')) run_setting(setting)
