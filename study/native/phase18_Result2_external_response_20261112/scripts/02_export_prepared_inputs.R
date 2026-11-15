suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)]);out<-file.path(Sys.getenv('CRRM_WORKSPACE'),'phase18_Result2_external_response_20261112')
objects<-readRDS(file.path(out,'data/processed/external_log_expression.rds'))
for(name in names(objects)){
  x<-objects[[name]]
  write.csv(data.frame(human_ensembl=rownames(x$x),x$x,check.names=FALSE),gzfile(file.path(out,'data/processed',paste0(name,'_log_expression.csv.gz'))),row.names=FALSE,na='')
  write.csv(x$meta,file.path(out,'data/processed',paste0(name,'_metadata.csv')),row.names=FALSE,na='')
}
write_json(names(objects),file.path(out,'data/processed/datasets.json'),pretty=TRUE)
print(names(objects))
