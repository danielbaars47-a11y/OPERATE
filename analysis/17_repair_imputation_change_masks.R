#!/usr/bin/env Rscript
# Repair change-mask schema only using existing frozen inputs/completed draws.
# No imputation, analysis fitting, or completed-data modification occurs.
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)) args[[1]] else 'OPERATE/results/first_visit_2026-10-27_v3'
read_data<-function(p) read.csv(gzfile(p),stringsAsFactors=FALSE,check.names=FALSE,na.strings='',colClasses=c(participant_id='character'))
for(setting in c('SEH','POLI')) {
 original<-read_data(file.path(root,paste0('imputation_input_',setting,'.csv.gz')))
 stopifnot(!anyDuplicated(names(original)),!anyDuplicated(original$participant_id))
 for(i in 1:50) {
  cc<-read_data(file.path(root,sprintf('imputed_%s_%02d.csv.gz',setting,i)))
  stopifnot(identical(original$participant_id,cc$participant_id))
  cols<-setdiff(names(original),'participant_id')
  flags<-setNames(lapply(cols,function(v) {
   a<-original[[v]];b<-cc[[v]]
   same<-(is.na(a)&is.na(b))|(!is.na(a)&!is.na(b)&as.character(a)==as.character(b))
   # CSV readers can choose integer/double classes independently. Compare
   # finite numeric values with a tight CSV round-trip tolerance.
   if(is.numeric(a)&&is.numeric(b)) {
    both<-!is.na(a)&!is.na(b)&is.finite(a)&is.finite(b)
    same[both]<-abs(a[both]-b[both])<=1e-12*pmax(1,abs(a[both]),abs(b[both]))
   }
   stopifnot(!anyNA(same));!same
  }),cols)
  result<-data.frame(participant_id=original$participant_id,flags,check.names=FALSE)
  result$.imp<-i
  stopifnot(!anyDuplicated(names(result)),sum(names(result)=='participant_id')==1L)
  dest<-file.path(root,sprintf('imputed_changed_%s_%02d.csv.gz',setting,i))
  con<-gzfile(dest,'wt');write.csv(result,con,row.names=FALSE,na='');close(con)
 }
 message(setting,': repaired 50 change masks; exactly one participant_id column each')
}
