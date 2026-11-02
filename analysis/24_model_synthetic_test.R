#!/usr/bin/env Rscript
# Nonclinical regression test: ignored evaluation values cannot affect training
# imputations or fitted coefficients; outcomes are never imputer inputs.
.libPaths(c('OPERATE/analysis/R-library',.libPaths()))
e<-new.env(parent=globalenv())
needed<-c('auc','metrics','trace_check','impute','fit_predict')
for(z in parse('OPERATE/analysis/24_model_pipeline.R')) if(is.call(z)&&identical(z[[1]],as.name('<-'))&&is.symbol(z[[2]])&&as.character(z[[2]])%in%needed)eval(z,e)
set.seed(71);n<-150L;nt<-110L
x<-data.frame(age=rnorm(n,50,10),male=rep(0:1,length.out=n),time=runif(n,1,15),RUQ=sample(0:1,n,TRUE))
x$time[c(10,60,130)]<-NA;y<-rbinom(nt,1,plogis(-1+.02*x$age[1:nt]))
z<-x;z$age[(nt+1):n]<-z$age[(nt+1):n]+1000
one<-e$impute(x,nt,5L,101,10L,30L);two<-e$impute(z,nt,5L,101,10L,30L)
for(i in 1:5)stopifnot(isTRUE(all.equal(one$data[[i]][1:nt,],two$data[[i]][1:nt,],tolerance=0)))
f1<-e$fit_predict(one$data,y,nt);f2<-e$fit_predict(two$data,y,nt)
stopifnot(isTRUE(all.equal(f1$coef,f2$coef,tolerance=0)),all(is.finite(f1$prediction)),all(is.finite(e$metrics(y,f1$prediction[1:nt]))))
cat('Synthetic ignored-row invariance and finite prediction/metric checks passed.\n')
