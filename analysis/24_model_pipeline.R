#!/usr/bin/env Rscript
# Vooraf gespecificeerde predictiemodellen, zonder aanvullende bronadjudicatie.
# Definitieve uitkomst uitsluitend voor fit/evaluatie; nooit in imputatie.
# Final50 MI;200 bootstrapresamples met opnieuw MI binnen elk trainingssample.
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)) normalizePath(args[1],mustWork=TRUE) else normalizePath('OPERATE/results/analysis_2026-10-27_no_adjudication',mustWork=TRUE)
.libPaths(c(file.path(dirname(dirname(root)),'analysis','R-library'),.libPaths()))
stopifnot(requireNamespace('mice',quietly=TRUE),requireNamespace('glmnet',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE))
config<-jsonlite::fromJSON(file.path(root,'model_specification.json'),simplifyVector=FALSE)
writecsv<-function(x,p)write.csv(x,p,row.names=FALSE,na='')
auc<-function(y,p){n1<-sum(y==1);n0<-sum(y==0);if(!n1||!n0)return(NA_real_);(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)}
metrics<-function(y,p){p<-pmin(1-1e-8,pmax(1e-8,p));lp<-qlogis(p);a<-tryCatch(suppressWarnings(glm(y~1+offset(lp),family=binomial())),error=function(e)NULL);s<-tryCatch(suppressWarnings(glm(y~lp,family=binomial())),error=function(e)NULL);if(is.null(a)||is.null(s)||!isTRUE(a$converged)||!isTRUE(s$converged)||any(!is.finite(c(coef(a),coef(s)))))stop('Calibration model failed convergence/finite coefficient check');c(AUC=auc(y,p),Brier=mean((y-p)^2),calibration_intercept=if(is.null(a))NA_real_ else unname(coef(a)[1]),calibration_slope=if(is.null(s))NA_real_ else unname(coef(s)[2]))}
# Split-Rhat for mean/SD traces; SD undefined for a single missing cell.
trace_check<-function(f,missingcols,w){rh<-function(z){h<-floor(nrow(z)/2);if(h<2||any(!is.finite(z)))return(NA_real_);z<-cbind(z[1:h,,drop=FALSE],z[(nrow(z)-h+1):nrow(z),,drop=FALSE]);W<-mean(apply(z,2,var));B<-h*var(colMeans(z));if(W==0)return(if(B==0)1 else Inf);sqrt(((h-1)*W/h+B/h)/W)};if(!length(missingcols))return(data.frame());do.call(rbind,lapply(missingcols,function(v)data.frame(variable=v,rhat_mean=rh(f$chainMean[v,,,drop=TRUE]),rhat_sd=if(sum(w[,v])>=2)rh(sqrt(f$chainVar[v,,,drop=TRUE])) else NA_real_,sd_applicable=sum(w[,v])>=2)))}
impute<-function(x,ntrain,m,seed,iterations,max_iterations){
 stopifnot(!any(grepl('outcome|finalY|gallstone|participant',names(x),ignore.case=TRUE)))
 missingcols<-names(x)[colSums(is.na(x))>0];if(!length(missingcols))return(list(data=replicate(m,x,simplify=FALSE),fit=NULL,convergence=data.frame(),iterations=0L))
 method<-setNames(rep('',ncol(x)),names(x));method[missingcols]<-'pmm'
 # Selected models have one numeric incomplete predictor: time since surgery.
 stopifnot(all(vapply(x[missingcols],is.numeric,logical(1))))
 pred<-matrix(1L,ncol(x),ncol(x),dimnames=list(names(x),names(x)));diag(pred)<-0;pred[setdiff(names(x),missingcols),]<-0
 set.seed(seed);initial<-x
 for(v in missingcols){donors<-x[[v]][seq_len(ntrain)];donors<-donors[is.finite(donors)];if(length(unique(donors))<5)stop('Insufficient training-only PMM support: ',v);initial[[v]][is.na(initial[[v]])]<-sample(donors,sum(is.na(initial[[v]])),replace=TRUE)}
 ignore<-seq_len(nrow(x))>ntrain;w<-is.na(x)
 f<-mice::mice(x,m=m,maxit=iterations,method=method,predictorMatrix=pred,ignore=ignore,data.init=initial,where=w,seed=seed,printFlag=FALSE,remove.constant=FALSE,remove.collinear=FALSE,remove_collinear=FALSE)
 check<-trace_check(f,missingcols,w)
 poor<-function(z)any(!is.finite(z$rhat_mean)|z$rhat_mean>1.1|(z$sd_applicable&(!is.finite(z$rhat_sd)|z$rhat_sd>1.1)))
 if(poor(check)&&max_iterations>iterations){f<-mice::mice.mids(f,maxit=max_iterations-iterations,printFlag=FALSE);check<-trace_check(f,missingcols,w)}
 if(any(!nzchar(f$method[missingcols])))stop('MICE dropped selected incomplete predictor')
 draws<-lapply(seq_len(m),function(i)mice::complete(f,i));if(!all(vapply(draws,function(z)all(is.finite(as.matrix(z))),logical(1))))stop('Nonfinite completed predictor')
 # Retain diagnostics, including remaining convergence concerns; no hidden fill.
 list(data=draws,fit=f,convergence=check,iterations=f$iteration,convergence_concern=poor(check))
}
fit_predict<-function(draws,y,ntrain,save_models=FALSE){
 pred<-matrix(NA_real_,nrow(draws[[1]]),length(draws));coefs<-NULL;fits<-list()
 for(i in seq_along(draws)){z<-as.matrix(draws[[i]]);storage.mode(z)<-'double';fit<-glmnet::glmnet(z[seq_len(ntrain),,drop=FALSE],y,alpha=0,lambda=1/ntrain,standardize=TRUE,intercept=TRUE,family='binomial');pred[,i]<-as.numeric(predict(fit,z,s=1/ntrain,type='response'));coefs<-cbind(coefs,as.matrix(coef(fit,s=1/ntrain))[,1]);if(save_models)fits[[i]]<-fit}
 list(prediction=rowMeans(pred),coef=rowMeans(coefs),coefs=coefs,fits=fits)
}
plot_results<-function(id,y,p,boot,apparent,corrected){
 pdf(file.path(root,paste0(id,'_discrimination.pdf')),width=3.5,height=3.3,family='Helvetica',bg='white',useDingbats=FALSE);dev.control(displaylist='enable');par(mar=c(4,4,1,1),las=1,mgp=c(2.4,.6,0),tcl=-.25,bg='white');ord<-order(p,decreasing=TRUE);yy<-y[ord];thresholds<-c(which(diff(p[ord])!=0),length(p));plot(c(0,(cumsum(yy==0)/sum(yy==0))[thresholds]),c(0,(cumsum(yy==1)/sum(yy==1))[thresholds]),type='l',lwd=1.6,col='#0072B2',xlim=c(0,1),ylim=c(0,1),xlab='False-positive rate',ylab='True-positive rate',axes=FALSE);axis(1);axis(2);abline(0,1,lty=2,col='grey50');legend('bottomright',legend=sprintf('OOB AUC %.3f\nCorrected AUC %.3f',auc(y,p),corrected['AUC']),bty='n',cex=.75);roc_record<-recordPlot();dev.off();png(file.path(root,paste0(id,'_discrimination.png')),width=3.5,height=3.3,units='in',res=300,bg='white');replayPlot(roc_record);dev.off()
 pdf(file.path(root,paste0(id,'_calibration.pdf')),width=3.5,height=3.3,family='Helvetica',bg='white',useDingbats=FALSE);dev.control(displaylist='enable');par(mar=c(4,4,1,1),las=1,mgp=c(2.4,.6,0),tcl=-.25,bg='white');g<-cut(rank(p,ties.method='first'),breaks=seq(0,length(p),length.out=9),include.lowest=TRUE);dd<-do.call(rbind,lapply(split(seq_along(y),g),function(ii){n<-length(ii);o<-mean(y[ii]);d<-1+1.96^2/n;center<-(o+1.96^2/(2*n))/d;half<-1.96*sqrt(o*(1-o)/n+1.96^2/(4*n^2))/d;c(p=mean(p[ii]),o=o,lo=center-half,hi=center+half,n=n)}));plot(dd[,'p'],dd[,'o'],xlim=c(0,1),ylim=c(0,1),pch=21,bg='#0072B2',xlab='Mean predicted probability',ylab='Observed event proportion',axes=FALSE);axis(1);axis(2);abline(0,1,lty=2,col='grey50');arrows(dd[,'p'],dd[,'lo'],dd[,'p'],dd[,'hi'],angle=90,code=3,length=.035);mtext('OOB calibration; Wilson 95% intervals',side=3,line=-.6,adj=0,cex=.65);cal_record<-recordPlot();dev.off();png(file.path(root,paste0(id,'_calibration.png')),width=3.5,height=3.3,units='in',res=300,bg='white');replayPlot(cal_record);dev.off();writecsv(as.data.frame(dd),file.path(root,paste0(id,'_calibration_bins.csv')))
}
run_model<-function(model,index){
 id<-model$id;message('Starting ',id)
 d<-read.csv(file.path(root,paste0('model_input_',id,'.csv')),check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(participant_id='character'))
 stopifnot(!anyDuplicated(d$participant_id),all(d$outcome %in% 0:1),sum(d$outcome==1)>0,sum(d$outcome==0)>0)
 y<-d$outcome;x<-d[,setdiff(names(d),c('participant_id','outcome')),drop=FALSE];stopifnot(ncol(x)>=2,all(vapply(x,is.numeric,logical(1))))
 n<-nrow(x);base_seed<-if(!is.null(model$imputation_seed))as.integer(model$imputation_seed) else 27102026L+index*100000L
 final<-impute(x,n,50L,base_seed,30L,100L)
 if(isTRUE(final$convergence_concern))stop('Final imputation convergence unresolved: ',id)
 f<-fit_predict(final$data,y,n,save_models=TRUE);apparent<-metrics(y,f$prediction)
 saveRDS(list(imputation=final,coefficient_draws=f$coefs,fitted_models=f$fits,specification=model,sessionInfo=sessionInfo()),file.path(root,paste0(id,'_final_models.rds')))
 writecsv(data.frame(participant_id=d$participant_id,outcome=y,prediction=f$prediction),file.path(root,paste0(id,'_predictions.csv')))
 for(i in 1:50){zz<-data.frame(participant_id=d$participant_id,.imp=i,final$data[[i]],check.names=FALSE);con<-gzfile(file.path(root,sprintf('%s_completed_%02d.csv.gz',id,i)),'wt');write.csv(zz,con,row.names=FALSE);close(con)}
 B<-200L;population_seed<-if(!is.null(model$population_seed))as.integer(model$population_seed) else base_seed
 set.seed(population_seed);samples<-replicate(B,sample.int(n,n,replace=TRUE),simplify=FALSE);saveRDS(samples,file.path(root,paste0(id,'_bootstrap_indices.rds')))
 oneboot<-function(b){tryCatch({ii<-samples[[b]];if(length(unique(y[ii]))<2)stop('Bootstrap has one outcome class');xx<-rbind(x[ii,,drop=FALSE],x);mi<-impute(xx,n,5L,base_seed+b,10L,100L);if(isTRUE(mi$convergence_concern))stop('Bootstrap MI convergence unresolved after100iterations');fp<-fit_predict(mi$data,y[ii],n);train<-metrics(y[ii],fp$prediction[seq_len(n)]);test<-metrics(y,fp$prediction[n+seq_len(n)]);list(ok=TRUE,train=train,test=test,coef=fp$coef,convergence_concern=isTRUE(mi$convergence_concern),iterations=mi$iterations,oob_prediction=replace(fp$prediction[n+seq_len(n)],unique(ii),NA_real_))},error=function(e)list(ok=FALSE,error=conditionMessage(e)))}
 cores<-as.integer(Sys.getenv('OPERATE_MODEL_CORES','2'));cores<-max(1,min(cores,4))
 boots<-parallel::mclapply(seq_len(B),oneboot,mc.cores=cores,mc.preschedule=TRUE,mc.set.seed=FALSE)
 saveRDS(boots,file.path(root,paste0(id,'_bootstrap.rds')))
 ok<-vapply(boots,function(z)isTRUE(z$ok),logical(1));success<-boots[ok]
 if(length(success)<180)stop('Fewer than 180 valid bootstrap replicates: ',id)
 train<-do.call(rbind,lapply(success,'[[','train'));test<-do.call(rbind,lapply(success,'[[','test'));stopifnot(all(is.finite(train)),all(is.finite(test)));optimism<-colMeans(train-test);corrected<-apparent-optimism
 oobmat<-do.call(cbind,lapply(success,'[[','oob_prediction'));oobn<-rowSums(!is.na(oobmat));oobp<-rowMeans(oobmat,na.rm=TRUE);stopifnot(all(oobn>0));oobmetric<-metrics(y,oobp)
 writecsv(data.frame(participant_id=d$participant_id,outcome=y,prediction_apparent=f$prediction,prediction_OOB=oobp,OOB_replicates=oobn),file.path(root,paste0(id,'_predictions.csv')))
 rawmetrics<-data.frame(replicate=which(ok),train,test,check.names=FALSE);names(rawmetrics)<-c('replicate',paste0('train_',colnames(train)),paste0('test_',colnames(test)));writecsv(rawmetrics,file.path(root,paste0(id,'_bootstrap_metrics.csv')))
 perf<-data.frame(metric=names(apparent),apparent=as.numeric(apparent),mean_optimism=as.numeric(optimism),optimism_corrected=as.numeric(corrected),OOB_averaged_prediction_metric=as.numeric(oobmetric),bootstrap_optimism_lower=apply(train-test,2,quantile,.025,na.rm=TRUE),bootstrap_optimism_upper=apply(train-test,2,quantile,.975,na.rm=TRUE))
 writecsv(perf,file.path(root,paste0(id,'_performance.csv')))
 bc<-do.call(cbind,lapply(success,'[[','coef'));ci<-t(apply(bc,1,quantile,c(.025,.975),na.rm=TRUE));coefsummary<-data.frame(term=names(f$coef),MI_mean_log_coefficient=as.numeric(f$coef),bootstrap_lower=ci[,1],bootstrap_upper=ci[,2],exp_coefficient=exp(f$coef),exp_lower=exp(ci[,1]),exp_upper=exp(ci[,2]))
 writecsv(data.frame(replicate=which(ok),t(bc),check.names=FALSE),file.path(root,paste0(id,'_bootstrap_coefficients.csv')))
 writecsv(coefsummary,file.path(root,paste0(id,'_coefficients.csv')))
 details<-do.call(rbind,lapply(seq_len(B),function(b){z<-boots[[b]];if(!isTRUE(z$ok))return(data.frame(replicate=b,status='failed',reason=z$error,convergence_concern=NA,iterations=NA));data.frame(replicate=b,status='completed',reason='',convergence_concern=z$convergence_concern,iterations=z$iterations)}));writecsv(details,file.path(root,paste0(id,'_bootstrap_status.csv')))
 plot_results(id,y,oobp,boots,apparent,corrected)
 result<-list(id=id,n=n,events=sum(y==1),predictor_df=ncol(x),m_final=50,m_bootstrap=5,bootstrap_requested=B,bootstrap_success=sum(ok),bootstrap_failed=sum(!ok),bootstrap_convergence_concerns=sum(vapply(success,function(z)isTRUE(z$convergence_concern),logical(1))),final_imputation_iterations=final$iterations,performance=perf,penalty='ridge alpha=0 lambda=1/n_training; training-only standardization',outcome_in_imputation=FALSE,OOB_min_replicates=min(oobn),OOB_max_replicates=max(oobn),assessment='Internally validated under recorded unresolved-source exclusions; no external validation or causal interpretation')
 jsonlite::write_json(result,file.path(root,paste0(id,'_summary.json')),pretty=TRUE,auto_unbox=TRUE,na='null');message('Completed ',id,': n=',n,' events=',sum(y==1),' bootstrap=',sum(ok),'/',B);result
}
results<-lapply(seq_along(config$models),function(i)run_model(config$models[[i]],i))
jsonlite::write_json(results,file.path(root,'model_results.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
