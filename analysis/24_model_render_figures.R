#!/usr/bin/env Rscript
# Figures regenerated only from saved predictions/results; no fitting or MI.
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)) args[1] else 'OPERATE/results/analysis_2026-10-27_no_adjudication'
files<-list.files(root,pattern='_predictions[.]csv$',full.names=TRUE)
for(file in files){
 id<-sub('_predictions[.]csv$','',basename(file));d<-read.csv(file,check.names=FALSE)
 if(!'prediction_OOB'%in%names(d))next
 y<-d$outcome;p<-d$prediction_OOB;stopifnot(all(is.finite(p)),all(d$OOB_replicates>0))
 perf<-read.csv(file.path(root,paste0(id,'_performance.csv')))
 n1<-sum(y==1);n0<-sum(y==0);auc<-(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)
 corrected<-perf$optimism_corrected[perf$metric=='AUC']
 g<-cut(rank(p,ties.method='first'),breaks=seq(0,length(p),length.out=9),include.lowest=TRUE)
 dd<-do.call(rbind,lapply(split(seq_along(y),g),function(ii){n<-length(ii);o<-mean(y[ii]);den<-1+1.96^2/n;center<-(o+1.96^2/(2*n))/den;half<-1.96*sqrt(o*(1-o)/n+1.96^2/(4*n^2))/den;c(predicted=mean(p[ii]),observed=o,lower=center-half,upper=center+half,n=n,events=sum(y[ii]))}))
 write.csv(dd,file.path(root,paste0(id,'_calibration_bins.csv')),row.names=FALSE)
 for(kind in c('discrimination','calibration'))for(extension in c('pdf','png')){
  dest<-file.path(root,paste0(id,'_',kind,'.',extension))
  if(extension=='pdf')pdf(dest,width=3.5,height=3.3,family='Helvetica',pointsize=9,bg='white',useDingbats=FALSE) else png(dest,width=3.5,height=3.3,units='in',res=300,pointsize=9,bg='white')
  par(mar=c(4.2,4.2,1.6,.8),las=1,mgp=c(2.6,.6,0),tcl=-.25,bg='white',fg='black',col.axis='black',col.lab='black',family='sans',bty='l')
  if(kind=='discrimination'){
   orderp<-order(p,decreasing=TRUE);yy<-y[orderp];ends<-c(which(diff(p[orderp])!=0),length(p))
   plot(c(0,(cumsum(yy==0)/n0)[ends]),c(0,(cumsum(yy==1)/n1)[ends]),type='l',lwd=1.5,col='#0072B2',xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i',xlab='False-positive rate',ylab='True-positive rate',axes=FALSE)
   axis(1);axis(2);abline(0,1,lty=2,col='grey55',lwd=.7)
   legend('bottomright',legend=c(sprintf('OOB AUC %.3f',auc),sprintf('Optimism-corrected AUC %.3f',corrected)),bty='n',cex=.75)
   mtext(sprintf('OOB predictions; n = %d, events = %d',length(y),n1),side=3,line=.35,adj=0,cex=.7)
  }else{
   plot(dd[,'predicted'],dd[,'observed'],xlim=c(0,1),ylim=c(0,1),pch=21,bg='#0072B2',col='black',cex=.8,xlab='Mean predicted probability',ylab='Observed event proportion',axes=FALSE)
   axis(1);axis(2);abline(0,1,lty=2,col='grey55',lwd=.7)
   arrows(dd[,'predicted'],dd[,'lower'],dd[,'predicted'],dd[,'upper'],angle=90,code=3,length=.025,lwd=.7)
   mtext('OOB calibration; Wilson 95% intervals',side=3,line=.35,adj=0,cex=.7)
  }
  dev.off()
 }
}
cat('Saved white-background vector and 300 dpi raster figures from existing outputs only.\n')
