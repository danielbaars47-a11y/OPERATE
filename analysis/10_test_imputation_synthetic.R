set.seed(20261027)
p<-tempfile('10_operate_synthetic_');dir.create(p)
.libPaths(c(file.path(getwd(),'OPERATE','analysis','R-library'),.libPaths()))
jsonlite::write_json(list(outcome_mapping_status='validated'),file.path(p,'configuration.json'),auto_unbox=TRUE)
n<-220L
x<-data.frame(participant_id=sprintf('P%04d',1:n),age_years_at_index=runif(n,20,70),gallstone_related_disease=rep(0:1,length.out=n),continuous=rnorm(n),binary=sample(c('a','b'),n,TRUE),nominal=sample(c('a','b','c'),n,TRUE),ordinal=sample(c('a','b','c'),n,TRUE))
vs<-c('continuous','binary','nominal','ordinal'); w<-data.frame(participant_id=x$participant_id)
for(v in vs) {w[[v]]<-FALSE;w[[v]][sample(n,20)]<-TRUE;x[w[[v]],v]<-NA}
e<-w;for(v in vs)e[[v]]<-TRUE
sp<-data.frame(variable=vs,method=c('pmm','logreg','polyreg','polr'),block='clinical',levels_json=c('','["a","b"]','["a","b","c"]','["a","b","c"]'),reason='synthetic',verification_state='current_machine_dictionary_verified')
wr<-function(x,f){con<-gzfile(f,'wt');write.csv(x,con,row.names=FALSE,na='');close(con)}
for(st in c('SEH','POLI')) {wr(x,file.path(p,paste0('imputation_input_',st,'.csv.gz')));wr(w,file.path(p,paste0('where_',st,'.csv.gz')));wr(e,file.path(p,paste0('eligibility_',st,'.csv.gz')));write.csv(sp,file.path(p,paste0('imputation_spec_',st,'.csv')),row.names=FALSE,na='')}

status<-system2(file.path(R.home('bin'),'Rscript'),c('OPERATE/analysis/10_impute_first_visit_v2.R',shQuote(p)))
stopifnot(status==0L)
for(st in c('SEH','POLI')) {
 inv<-read.csv(file.path(p,paste0('imputation_inventory_',st,'.csv')))
 stopifnot(setequal(inv$variable,vs),all(inv$status=='imputed'))
 for(i in 1:50) {
  cc<-read.csv(gzfile(file.path(p,sprintf('imputed_%s_%02d.csv.gz',st,i))),na.strings='')
  stopifnot(identical(cc$participant_id,x$participant_id),identical(cc$gallstone_related_disease,x$gallstone_related_disease),all(complete.cases(cc)))
 }
}
cat('All four methods passed 50-dataset synthetic checks; fixture:',p,'\n')
