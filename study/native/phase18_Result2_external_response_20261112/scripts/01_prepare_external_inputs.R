suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db);library(org.Mm.eg.db);library(edgeR);library(jsonlite);library(digest)})
args<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',args[grepl('^--file=',args)]);OUT<-file.path(Sys.getenv('CRRM_WORKSPACE'),'phase18_Result2_external_response_20261112')
PROJECT<-file.path(Sys.getenv('CRRM_DATA_ROOT'),"");REV<-file.path(PROJECT,'revisions/crrm_revision_20261111')
cfg<-fromJSON(file.path(OUT,'parameters/analysis_plan_v1.json'));receipt<-fromJSON(file.path(OUT,'parameters/declaration_receipt.json'))
stopifnot(digest(file=file.path(OUT,'parameters/analysis_plan_v1.json'),algo='sha256')==receipt$plan_sha256)
for(i in seq_len(nrow(cfg$inputs)))stopifnot(digest(file=cfg$inputs$path[i],algo='sha256')==cfg$inputs$sha256[i])
core<-read.csv(file.path(OUT,'parameters/frozen_CRRM388.csv'),check.names=FALSE);stopifnot(nrow(core)==388)
orth<-read.csv(file.path(REV,'phase02_crrm_reconstruction/mapping/babelgene22.9_mouse_human_one_to_one.csv'),stringsAsFactors=FALSE)
rawdir<-file.path(PROJECT,'data/raw/postfreeze_cross_system_validation');annodir<-file.path(PROJECT,'data/metadata/postfreeze_cross_system_validation')
objects<-list();mapping_rows<-list();coverage_rows<-list();manifest_rows<-list();norm_rows<-list()
write_tab<-function(x,name)write.csv(x,file.path(OUT,'tables',name),row.names=FALSE,na='')
cleanid<-function(x)sub('[.][0-9]+$','',as.character(x))
truth<-function(x)tolower(as.character(x))=='true'
series<-function(path){
  con<-gzfile(path,'rt');lines<-readLines(con,warn=FALSE);close(con)
  b<-grep('^!series_matrix_table_begin',lines);e<-grep('^!series_matrix_table_end',lines);stopifnot(length(b)==1,length(e)==1)
  d<-read.delim(textConnection(paste(lines[(b+1):(e-1)],collapse='\n')),check.names=FALSE,stringsAsFactors=FALSE,quote='"')
  x<-as.matrix(d[,-1]);storage.mode(x)<-'double';rownames(x)<-as.character(d[[1]])
  meta<-function(key){s<-lines[startsWith(lines,key)][1];v<-strsplit(s,'\t',fixed=TRUE)[[1]][-1];gsub('^"|"$','',v)}
  m<-data.frame(sample_id=meta('!Sample_geo_accession'),title=meta('!Sample_title'),stringsAsFactors=FALSE)
  stopifnot(identical(colnames(x),m$sample_id),!anyDuplicated(rownames(x)))
  list(x=x,meta=m)
}
platform<-function(path){
  con<-if(grepl('[.]gz$',path))gzfile(path,'rt')else file(path,'rt');s<-readLines(con,warn=FALSE);close(con)
  h<-grep('^ID\t',s)[1];end<-grep('^!platform_table_end|^!annotation_table_end',s);last<-if(length(end)&&end[1]>h)end[1]-1 else length(s)
  stopifnot(!is.na(h));read.delim(textConnection(paste(s[h:last],collapse='\n')),check.names=FALSE,stringsAsFactors=FALSE,quote='',fill=TRUE)
}
native_mouse<-AnnotationDbi::select(org.Mm.eg.db,keys=keys(org.Mm.eg.db,keytype='SYMBOL'),keytype='SYMBOL',columns='ENSEMBL')
native_mouse<-unique(native_mouse[!is.na(native_mouse$ENSEMBL),c('SYMBOL','ENSEMBL')]);official<-split(native_mouse$ENSEMBL,native_mouse$SYMBOL)
official_env<-list2env(official,hash=TRUE,parent=emptyenv())
alias_cache<-list()
mouse_project<-function(x,feature_ids,symbols,dataset){
  symbols[is.na(symbols)]<-'';pieces<-strsplit(symbols,'\\s*(?:///|;|,)\\s*',perl=TRUE);tokens<-unique(trimws(unlist(pieces)));tokens<-tokens[nzchar(tokens)]
  unknown<-setdiff(tokens,names(official));valid_alias<-intersect(unknown,keys(org.Mm.eg.db,keytype='ALIAS'))
  alias<-if(length(valid_alias))AnnotationDbi::select(org.Mm.eg.db,keys=valid_alias,keytype='ALIAS',columns='ENSEMBL')else data.frame(ALIAS=character(),ENSEMBL=character())
  alias<-alias[!is.na(alias$ENSEMBL),,drop=FALSE];alook<-split(alias$ENSEMBL,alias$ALIAS);alias_env<-list2env(alook,hash=TRUE,parent=emptyenv())
  getgenes<-function(token){v<-get0(token,envir=official_env,inherits=FALSE);if(is.null(v))v<-get0(token,envir=alias_env,inherits=FALSE);unique(v[!is.na(v)])}
  ids<-vapply(pieces,function(pp){pp<-unique(trimws(pp[nzchar(trimws(pp))]));if(!length(pp))return('');links<-lapply(pp,getgenes);if(any(lengths(links)!=1))return('');g<-unique(unlist(links));if(length(g)==1)g else ''},character(1))
  h<-orth$human_ensembl[match(ids,orth$ensembl)]
  audit<-data.frame(dataset=dataset,feature_id=as.character(feature_ids),source_symbol=symbols,native_ensembl=ids,human_ensembl=h,source_resolves_one_gene=nzchar(ids),projection_eligible=!is.na(h),stringsAsFactors=FALSE)
  audit<-audit[audit$feature_id%in%rownames(x),,drop=FALSE];mapping_rows[[dataset]]<<-audit
  complete_features<-rownames(x)[rowSums(is.finite(x))==ncol(x)];audit$complete_expression<-audit$feature_id%in%complete_features;mapping_rows[[dataset]]<<-audit
  audit<-unique(audit[audit$projection_eligible&audit$complete_expression&audit$human_ensembl%in%core$human_ensembl,c('feature_id','human_ensembl')])
  # Multiple unique probes/transcripts for the same native gene are collapsed without outcome information.
  groups<-split(audit$feature_id,audit$human_ensembl)
  out<-do.call(rbind,lapply(groups,function(ids){v<-x[unique(ids),,drop=FALSE];if(nrow(v)==1)as.numeric(v[1,])else apply(v,2,median,na.rm=FALSE)}));rownames(out)<-names(groups);colnames(out)<-colnames(x)
  out
}
register<-function(dataset,x,m,reference,method,projection){
  stopifnot(!anyDuplicated(rownames(x)),identical(colnames(x),m$sample_id),all(is.finite(x)),all(rownames(x)%in%core$human_ensembl),sum(reference)>=2)
  m$dataset<-dataset;m$reference_for_scaling<-reference
  mu<-rowMeans(x[,reference,drop=FALSE]);sdev<-apply(x[,reference,drop=FALSE],1,sd);present<-match(core$human_ensembl,rownames(x));usable<-!is.na(present);usable[usable]<-is.finite(sdev[present[usable]])&sdev[present[usable]]>0
  coverage_rows[[dataset]]<<-data.frame(dataset=dataset,human_ensembl=core$human_ensembl,display_symbol=core$display_symbol,frozen_direction=core$frozen_direction,present=!is.na(present),scored=usable,reference_sd=ifelse(is.na(present),NA,sdev[present]),projection=projection)
  objects[[dataset]]<<-list(x=x,meta=m,reference=reference,method=method,projection=projection)
  manifest_rows[[dataset]]<<-m
  norm_rows[[dataset]]<<-data.frame(dataset=dataset,method=method,n_native_features=attr(x,'n_native')%||%NA,n_mapped_core=nrow(x),n_scored=sum(usable),n_reference=sum(reference),projection=projection)
  message(dataset,':',ncol(x),'samples,',sum(usable),'/388variable mapped members')
}
`%||%`<-function(x,y)if(is.null(x))y else x
maybe_log<-function(x){q<-unname(quantile(x[is.finite(x)],.99));if(q>100)log2(pmax(x,0)+1)else x}

g<-series(file.path(rawdir,'GSE112956_series_matrix.txt.gz'));a<-platform(file.path(annodir,'GPL1261.annot.gz'))
m<-g$meta;m$strain<-sub('_.*','',m$title);m$condition<-ifelse(grepl('vehicle',m$title),'Control','RSV');m$time_d<-ifelse(grepl('5 day RSV',m$title),5,1);m$unit_id<-m$sample_id
counts112<-as.data.frame(table(m$strain,m$condition,m$time_d));write_tab(counts112,'00_GSE112956_actual_cell_counts.csv')
stopifnot(sum(m$condition=='Control')==14,sum(m$condition=='RSV'&m$time_d==1)==15,sum(m$time_d==5)==15)
write_json(list(declared_before_scoring=TRUE,reason='One DBA2J vehicle profile absent in deposited matrix relative to the paper nominal n=3 per cell',actual_counts=list(vehicle=14,RSV_D1=15,RSV_D5=15),primary_estimand='Equal-strain mean of5RSVD1-minus-vehicle differences',model='Saturated10-cell strain-by-D1-condition means with equal1/5contrast weights;pooled residual tCI',plan_core_and_thresholds_unchanged=TRUE),file.path(OUT,'parameters/metadata_amendment_GSE112956.json'),pretty=TRUE,auto_unbox=TRUE)
x<-mouse_project(maybe_log(g$x),a$ID,a[[intersect(c('Gene symbol','Gene Symbol'),names(a))[1]]],'GSE112956');register('GSE112956',x,m,m$time_d==1,'Deposited array expression;log2(x+1)onlyif99th percentile>100;median unambiguous probes','Frozen one-to-one mouse-human graph')

g<-series(file.path(rawdir,'GSE119849_series_matrix.txt.gz'));a<-platform(file.path(annodir,'GPL17543_platform.txt'));keep<-grepl('^RSV_',g$meta$title);g$x<-g$x[,keep,drop=FALSE];m<-g$meta[keep,,drop=FALSE];m$condition<-ifelse(grepl('RSV_Disease',m$title),'RSV','Control');m$strain<-'C57BL/6J';m$time_d<-NA_real_;m$unit_id<-m$sample_id
stopifnot(sum(m$condition=='RSV')==10,sum(m$condition=='Control')==10)
x<-mouse_project(maybe_log(g$x),a$ID,a$Symbol,'GSE119849');register('GSE119849',x,m,rep(TRUE,nrow(m)),'Deposited normalized array expression;median unambiguous probes','Frozen one-to-one mouse-human graph')

base<-file.path(PROJECT,'data/processed/next_batch/GSE140226_Mouse_lung_processed');raw<-read.csv(file.path(base,'GSE140226_Mouse_lung_processed_expression_matrix.csv'),check.names=FALSE);a<-read.csv(file.path(base,'GSE140226_Mouse_lung_processed_gene_annotation.csv'));m<-read.csv(file.path(base,'GSE140226_Mouse_lung_processed_sample_metadata.csv'))
m<-m[m$perturbation=='None_or_vehicle',];m$condition<-ifelse(m$condition=='RSV','RSV','Control');m$time_d<-ifelse(grepl('_5d_',m$sample_id),5,7);m$strain<-'BALB/c';m$unit_id<-m$sample_id
x<-as.matrix(raw[,-1]);storage.mode(x)<-'double';rownames(x)<-raw[[1]];x<-log2(x[,m$sample_id,drop=FALSE]+1)
x<-mouse_project(x,a$feature_id,a$gene_symbol,'GSE140226');register('GSE140226',x,m,m$time_d==5,'log2(processed transcript abundance+1);median unambiguous transcripts','Frozen one-to-one mouse-human graph')

g<-series(file.path(rawdir,'GSE18170_series_matrix.txt.gz'));a<-platform(file.path(annodir,'GPL7202.annot.gz'));keep<-grepl('^WT',g$meta$title);m<-g$meta[keep,,drop=FALSE];m$condition<-ifelse(grepl('RSV',m$title),'RSV','Control');m$strain<-'WT';m$time_d<-NA_real_;m$unit_id<-m$sample_id
x<-mouse_project(maybe_log(g$x[,keep,drop=FALSE]),a$ID,a[[intersect(c('Gene symbol','Gene Symbol'),names(a))[1]]],'GSE18170');register('GSE18170',x,m,rep(TRUE,nrow(m)),'Deposited array expression;log2(x+1)onlyif99th percentile>100;median unambiguous probes','Frozen one-to-one mouse-human graph')

base<-file.path(PROJECT,'data/processed/postfreeze_cross_system_validation/PRJNA1271255_Zenodo15953852');raw<-read.csv(gzfile(file.path(base,'raw_counts_gene_id.csv.gz')),check.names=FALSE);a<-read.csv(file.path(base,'gene_annotation.csv'),check.names=FALSE);m0<-read.csv(file.path(base,'sample_metadata.csv'))
cc<-as.matrix(raw[,-1]);storage.mode(cc)<-'double';rownames(cc)<-as.character(raw[[1]])
symbols<-toupper(trimws(a$gene_name));human<-toupper(core$display_symbol);source_degree<-table(symbols);eligible<-symbols%in%human & as.integer(source_degree[symbols])==1 & !duplicated(toupper(core$display_symbol))[match(symbols,human)]
ha<-data.frame(feature_id=as.character(a$gene_id),source_symbol=a$gene_name,human_ensembl=core$human_ensembl[match(symbols,human)],projection_eligible=eligible);ha<-ha[!is.na(ha$projection_eligible)&ha$projection_eligible,];mapping_rows[['Hamster_symbols']]<-ha
for(tissue in c('Lung','Whole blood')){
  m<-m0[m0$measurement=='count'&m0$tissue==tissue&m0$condition%in%c('NC','MOCK'),];m$condition<-ifelse(m$condition=='NC','Control','RSV');m$time_d<-5;m$strain<-'Syrian hamster';m$unit_id<-paste(m$condition,m$animal_label,sep='_')
  d<-calcNormFactors(DGEList(cc[,m$sample_id,drop=FALSE]),method='TMM');lcpm<-cpm(d,log=TRUE,prior.count=1);x<-lcpm[ha$feature_id,,drop=FALSE];rownames(x)<-ha$human_ensembl
  dataset<-if(tissue=='Lung')'PRJNA1271255_Lung'else'PRJNA1271255_Blood';register(dataset,x,m,rep(TRUE,nrow(m)),'TMMlog2CPM(prior.count1);NC/MOCKonly','Unique author-symbol projection;not direct evolutionary orthology')
}

tarfile<-file.path(rawdir,'GSE247298_RAW.tar');files<-untar(tarfile,list=TRUE);tmp<-file.path(OUT,'data/processed/GSE247298_transcripts');dir.create(tmp,showWarnings=FALSE);untar(tarfile,exdir=tmp)
first<-read.delim(gzfile(file.path(tmp,files[1])),check.names=FALSE);ids<-unique(cleanid(first$target_id));valid<-intersect(ids,keys(org.Hs.eg.db,keytype='REFSEQ'))
links<-AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='REFSEQ',columns='ENSEMBL');links<-unique(links[!is.na(links$ENSEMBL),c('REFSEQ','ENSEMBL')]);degree<-table(links$REFSEQ);links<-links[degree[links$REFSEQ]==1,];mapping_rows[['GSE247298_RefSeq']]<-links
mat<-vapply(files,function(filename){d<-read.delim(gzfile(file.path(tmp,filename)));ref<-cleanid(d$target_id);ens<-links$ENSEMBL[match(ref,links$REFSEQ)];ok<-!is.na(ens);v<-rowsum(matrix(d$tpm[ok],ncol=1),group=ens[ok],reorder=FALSE);out<-rep(NA_real_,nrow(core));names(out)<-core$human_ensembl;hit<-match(names(out),rownames(v));out[!is.na(hit)]<-log2(v[hit[!is.na(hit)],1]+1);out},numeric(nrow(core)))
rownames(mat)<-core$human_ensembl;colnames(mat)<-sub('[.]tsv[.]gz$','',sub('^GSM[0-9]+_','',files));mat<-mat[rowSums(is.finite(mat))==ncol(mat),,drop=FALSE]
m<-data.frame(sample_id=colnames(mat),title=colnames(mat),stringsAsFactors=FALSE);m$donor_id<-sub('_([0-9]+)hrs$','',sub('^RSV_','',m$sample_id));m$background<-ifelse(grepl('^CF',m$donor_id),'CF','Non-CF');m$time_h<-as.integer(sub('.*_([0-9]+)hrs$','\\1',m$sample_id));m$condition<-ifelse(m$time_h==0,'Baseline','RSV');m$unit_id<-m$donor_id
stopifnot(nrow(m)==48,length(unique(m$donor_id))==8,all(table(m$donor_id)==6),setequal(unique(m$time_h),c(0,12,24,48,72,120)))
register('GSE247298',mat,m,m$background=='Non-CF','Sum TPM over unique RefSeq-to-human Ensembl links;log2(TPM+1)','Unique full-reference human RefSeq links')

paths<-setNames(file.path(PROJECT,'data/raw',c('GSE155237','GSE166161'),paste0(c('GSE155237','GSE166161'),'_gene_counts.csv.gz')),c('GSE155237','GSE166161'))
native<-lapply(paths,function(p){d<-read.csv(gzfile(p),check.names=FALSE);x<-as.matrix(d[,-1]);storage.mode(x)<-'double';rownames(x)<-cleanid(d[[1]]);stopifnot(!anyDuplicated(rownames(x)));x})
stopifnot(setequal(rownames(native[[1]]),rownames(native[[2]])));native[[2]]<-native[[2]][rownames(native[[1]]),,drop=FALSE]
m<-read.csv(file.path(REV,'phase01_sample_manifest/tables/07_human_challenge_91_unique_specimens.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM');m$sample_id<-m$record_id;m$donor_id<-m$donor_id_canonical;m$unit_id<-m$donor_id;m$time_d<-m$timepoint_day;m$condition<-m$eventual_infection_outcome;m$primary<-truth(m$include_primary_D3)
cc<-vapply(seq_len(nrow(m)),function(i)native[[m$accession[i]]][,m$expression_column[i]],numeric(nrow(native[[1]])));rownames(cc)<-rownames(native[[1]]);colnames(cc)<-m$sample_id
refcol<-which(m$primary&m$time_d==0)[1];d<-calcNormFactors(DGEList(cc),method='TMM',refColumn=refcol);lcpm<-cpm(d,log=TRUE,prior.count=1);hit<-intersect(core$human_ensembl,rownames(lcpm));x<-lcpm[hit,,drop=FALSE]
stopifnot(nrow(m)==91,sum(m$primary)==76,length(unique(m$donor_id[m$primary]))==38)
register('Human_challenge',x,m,m$primary,'TMMlog2CPM(prior.count1),fixed primary-baseline TMMreference;91unique specimens,76primaryscalingprofiles','Direct human Ensembl IDs;13shared baselines deduplicated')
manifest<-function(rows){cols<-unique(unlist(lapply(rows,names)));do.call(rbind,lapply(rows,function(x){for(k in setdiff(cols,names(x)))x[[k]]<-NA;x[,cols,drop=FALSE]}))}
write_tab(manifest(manifest_rows),'01_external_sample_manifest.csv');write_tab(do.call(rbind,coverage_rows),'02_member_coverage.csv');write_tab(do.call(rbind,norm_rows),'03_normalization_and_projection.csv')
for(name in names(mapping_rows))write.csv(mapping_rows[[name]],gzfile(file.path(OUT,'mapping',paste0(name,'_mapping.csv.gz'))),row.names=FALSE,na='')
coverage<-do.call(rbind,coverage_rows);summary<-do.call(rbind,lapply(split(coverage,coverage$dataset),function(c){data.frame(dataset=c$dataset[1],n_present=sum(c$present),n_scored=sum(c$scored),n_up=sum(c$scored&c$frozen_direction=='up'),n_down=sum(c$scored&c$frozen_direction=='down'),eligible_primary=sum(c$scored)/388>=.7&&sum(c$scored&c$frozen_direction=='up')/280>=.5&&sum(c$scored&c$frozen_direction=='down')/108>=.5)}));write_tab(summary,'04_coverage_summary.csv')
saveRDS(objects,file.path(OUT,'data/processed/external_log_expression.rds'));capture.output(sessionInfo(),file=file.path(OUT,'logs/input_preparation_R_sessionInfo.txt'))
write_json(list(status='PREPARED',datasets=names(objects),coverage=summary,no_CRRM_gene_selection=TRUE,score_outcomes_not_used_for_mapping=TRUE,plan_sha256=receipt$plan_sha256),file.path(OUT,'logs/input_preparation_summary.json'),pretty=TRUE,auto_unbox=TRUE)
print(summary)
