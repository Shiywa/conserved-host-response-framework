suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db);library(jsonlite);library(digest)})
args<-commandArgs(trailingOnly=FALSE);f<-sub('^--file=','',args[grepl('^--file=',args)]);OUT<-normalizePath(file.path(dirname(f),'..'))
ROOT<-file.path(Sys.getenv('CRRM_DATA_ROOT'),"")
cfg<-fromJSON(file.path(OUT,'parameters/analysis_plan_v1.json'));receipt<-fromJSON(file.path(OUT,'parameters/declaration_receipt.json'))
stopifnot(digest(file=file.path(OUT,'parameters/analysis_plan_v1.json'),algo='sha256')==receipt$plan_sha256)
manifest<-read.csv(file.path(OUT,'data/clinical_sample_manifest.csv'),check.names=FALSE);manifest$primary<-tolower(as.character(manifest$primary))=='true';core<-read.csv(file.path(OUT,'parameters/frozen_CRRM388.csv'))
read_table<-function(path,begin,end){
  con<-gzfile(path,'rt');s<-readLines(con,warn=FALSE);close(con);b<-grep(begin,s)[1];e<-grep(end,s)[1]
  stopifnot(!is.na(b),!is.na(e));read.delim(textConnection(paste(s[(b+1):(e-1)],collapse='\n')),check.names=FALSE,stringsAsFactors=FALSE,quote='"')
}
official<-unique(AnnotationDbi::select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,keytype='SYMBOL'),keytype='SYMBOL',columns='ENSEMBL'))
aliases<-unique(AnnotationDbi::select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,keytype='ALIAS'),keytype='ALIAS',columns='ENSEMBL'))
official<-official[!is.na(official$ENSEMBL),];aliases<-aliases[!is.na(aliases$ENSEMBL),]
oe<-list2env(split(official$ENSEMBL,official$SYMBOL),hash=TRUE,parent=emptyenv());ae<-list2env(split(aliases$ENSEMBL,aliases$ALIAS),hash=TRUE,parent=emptyenv())
resolve<-function(s){
  if(is.na(s)||!nzchar(s))return('');tokens<-unique(trimws(strsplit(s,'\\s*(?:///|;|,)\\s*',perl=TRUE)[[1]]));tokens<-tokens[nzchar(tokens)]
  if(!length(tokens))return('');links<-lapply(tokens,function(t){v<-get0(t,envir=oe,inherits=FALSE);if(is.null(v))v<-get0(t,envir=ae,inherits=FALSE);unique(v)})
  if(any(lengths(links)!=1))return('');v<-unique(unlist(links));if(length(v)==1)v else ''
}
objects<-list();summaries<-list()
for(dataset in c('GSE38900_GPL6884','GSE42026')){
  path<-if(dataset=='GSE42026')file.path(OUT,'data/GSE42026_series_matrix.txt.gz')else file.path(ROOT,'data/raw/GSE38900/GSE38900-GPL6884_series_matrix.txt.gz')
  plat<-if(dataset=='GSE42026')'GPL6947'else'GPL6884'
  ap<-if(dataset=='GSE42026')file.path(OUT,'data/GPL6947.annot.gz')else file.path(ROOT,'data/metadata/cross_virus_context/platform_annotation/GPL6884.annot.gz')
  d<-read_table(path,'^!series_matrix_table_begin','^!series_matrix_table_end');x<-as.matrix(d[,-1]);storage.mode(x)<-'double';rownames(x)<-d[[1]]
  m<-manifest[manifest$dataset==dataset&manifest$primary,];x<-x[,m$geo_accession,drop=FALSE];stopifnot(!anyDuplicated(m$participant),all(is.finite(x)))
  q99<-unname(quantile(x,.99));neg<-sum(x<0);transformed<-q99>100
  if(transformed)x<-log2(pmax(x,0)+1)
  con<-gzfile(ap,'rt');s<-readLines(con,warn=FALSE);close(con);b<-grep('^ID\t',s)[1];e<-grep('^!annotation_table_end|^!platform_table_end',s)[1]
  if(is.na(e))e<-length(s)+1
  a<-read.delim(textConnection(paste(s[b:(e-1)],collapse='\n')),check.names=FALSE,quote='',fill=TRUE,stringsAsFactors=FALSE)
  symbols<-a[['Gene symbol']];unique_symbols<-unique(symbols);mapping<-vapply(unique_symbols,resolve,character(1));ens<-unname(mapping[match(symbols,unique_symbols)])
  audit<-data.frame(platform=plat,feature_id=a$ID,source_symbol=symbols,human_ensembl=ens,eligible=nzchar(ens),stringsAsFactors=FALSE);audit<-audit[audit$feature_id%in%rownames(x),]
  write.csv(audit,gzfile(file.path(OUT,'data',paste0(dataset,'_full_reference_mapping.csv.gz'))),row.names=FALSE,na='')
  groups<-split(audit$feature_id[audit$eligible],audit$human_ensembl[audit$eligible]);y<-do.call(rbind,lapply(groups,function(ids){v<-x[unique(ids),,drop=FALSE];if(nrow(v)==1)as.numeric(v[1,])else apply(v,2,median)}))
  rownames(y)<-names(groups);colnames(y)<-m$geo_accession
  write.csv(y,gzfile(file.path(OUT,'data',paste0(dataset,'_log_expression.csv.gz'))),row.names=TRUE)
  write.csv(m,file.path(OUT,'data',paste0(dataset,'_sample_metadata.csv')),row.names=FALSE,na='')
  objects[[dataset]]<-list(x=y,meta=m)
  summaries[[dataset]]<-data.frame(dataset=dataset,platform=plat,samples=ncol(y),genes=nrow(y),frozen_mapped=sum(core$human_ensembl%in%rownames(y)),q99_native=q99,negative_native=neg,transformation=if(transformed)'log2(max(x,0)+1)'else'Deposited normalized scale',probe_collapse='Median of unambiguous probes; official symbol priority then unique aliases; one Ensembl gene required')
  message(dataset,':',ncol(y),'samples;',sum(core$human_ensembl%in%rownames(y)),'frozen members')
}
saveRDS(objects,file.path(OUT,'data/array_expression.rds'));write.csv(do.call(rbind,summaries),file.path(OUT,'audit/array_preparation.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(OUT,'audit/array_preparation_R_session.txt'))
