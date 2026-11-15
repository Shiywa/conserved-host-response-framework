#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(edgeR);library(jsonlite);library(digest)})
R<-file.path(Sys.getenv('CRRM_WORKSPACE'),"phase06_external_validation_figure1")
cfg<-fromJSON(file.path(R,'parameters/analysis_plan_v1.json'));receipt<-fromJSON(file.path(R,'parameters/declaration_receipt.json'));impl<-fromJSON(file.path(R,'parameters/count_implementation_receipt.json'))
stopifnot(digest(file=file.path(R,'parameters/analysis_plan_v1.json'),algo='sha256')==receipt$config_sha256,receipt$config_sha256==impl$analysis_plan_sha256)
for(k in names(cfg$inputs))stopifnot(digest(file=cfg$inputs[[k]]$path,algo='sha256')==cfg$inputs[[k]]$sha256)
for(i in seq_len(nrow(impl$processed_matrices))){x<-impl$processed_matrices[i,];stopifnot(digest(file=x$path,algo='sha256')==x$sha256,digest(file=x$manifest,algo='sha256')==x$manifest_sha256)}
members<-list(priority388=read.csv(cfg$inputs$priority388$path,check.names=FALSE),broad1980=read.csv(cfg$inputs$broad1980$path,check.names=FALSE))
stopifnot(nrow(members$priority388)==388,nrow(members$broad1980)==1980)
scores_all<-list();pairs_all<-list();coverage_all<-list();norm_all<-list();gene_all<-list();objects<-list()
write_table<-function(x,f)write.csv(x,file.path(R,'tables',f),row.names=FALSE,na='')
cohorts<-strsplit(Sys.getenv('CRRM_EXTERNAL_COHORTS',unset='GSE270463,GSE309353'),',',fixed=TRUE)[[1]]
for(cohort in cohorts){
 raw<-read.csv(gzfile(file.path(R,'data/processed',paste0(cohort,'_native_counts.csv.gz'))),check.names=FALSE)
 counts<-as.matrix(raw[,-1]);storage.mode(counts)<-'double';rownames(counts)<-raw$feature_id
 meta<-read.csv(file.path(R,'data/processed',paste0(cohort,'_selected_manifest.csv')),check.names=FALSE);meta$donor_id<-as.character(meta$donor_id)
 stopifnot(identical(colnames(counts),meta$gsm),!anyDuplicated(rownames(counts)),all(is.finite(counts)),all(counts>=0))
 group<-paste(meta$donor_id,meta$condition,sep='__');ug<-unique(group)
 pooled<-vapply(ug,function(g)rowSums(counts[,group==g,drop=FALSE]),numeric(nrow(counts)));rownames(pooled)<-rownames(counts)
 pm<-meta[match(ug,group),c('dataset','donor_id','condition','time_h')];pm$profile_id<-ug;pm$n_culture_profiles<-as.integer(table(group)[ug]);rownames(pm)<-NULL
 methods<-list(primary=list(counts=pooled,meta=pm))
 if(cohort=='GSE309353'){
  sm<-meta[,c('dataset','donor_id','condition','time_h')];sm$profile_id<-meta$gsm;sm$n_culture_profiles<-1
  methods$culture_profile_mean<-list(counts=counts,meta=sm)
 }
 for(method in names(methods)){
  cc<-methods[[method]]$counts;mm<-methods[[method]]$meta
  dge<-calcNormFactors(DGEList(cc),method='TMM');lcpm<-cpm(dge,log=TRUE,prior.count=1)
  norm_all[[paste(cohort,method)]]<-cbind(data.frame(cohort=cohort,method=method),mm,dge$samples)
  for(set in names(members)){
   genes<-members[[set]];sign<-ifelse(genes$frozen_direction=='up',1,-1);hit<-match(genes$human_ensembl,rownames(lcpm));present<-!is.na(hit)
   sdev<-rep(NA_real_,nrow(genes));nz<-rep(NA_integer_,nrow(genes));sdev[present]<-apply(lcpm[hit[present],,drop=FALSE],1,sd);nz[present]<-rowSums(cc[hit[present],,drop=FALSE]>0)
   usable<-present&is.finite(sdev)&sdev>0&!is.na(nz)&nz>0;stopifnot(any(usable))
   coverage_all[[paste(cohort,method,set)]]<-data.frame(cohort=cohort,method=method,set=set,human_ensembl=genes$human_ensembl,display_symbol=genes$display_symbol,frozen_direction=genes$frozen_direction,direction_multiplier=sign,present=present,n_profiles_nonzero=nz,sd=sdev,scored=usable)
   xx<-lcpm[hit[usable],,drop=FALSE];z<-sweep(sweep(xx,1,rowMeans(xx),'-'),1,sdev[usable],'/');score<-colMeans(sweep(z,1,sign[usable],'*'))
   ss<-cbind(data.frame(cohort=cohort,method=method,set=set),mm,score=unname(score));scores_all[[paste(cohort,method,set)]]<-ss
   pair<-lapply(sort(unique(mm$donor_id)),function(d){
    a<-ss$score[ss$donor_id==d&ss$condition=='mock'];b<-ss$score[ss$donor_id==d&ss$condition=='RSV'];stopifnot(length(a)==length(b),length(a)%in%c(1,2))
    data.frame(cohort=cohort,method=method,set=set,donor_id=d,time_h=unique(mm$time_h),n_mock_profiles=length(a),n_RSV_profiles=length(b),score_mock=mean(a),score_RSV=mean(b),delta=mean(b)-mean(a))
   });pairs_all[[paste(cohort,method,set)]]<-do.call(rbind,pair)
   if(method=='primary'){
    write_table(data.frame(human_ensembl=genes$human_ensembl[usable],direction_multiplier=sign[usable],z,check.names=FALSE),paste0('gene_zscores_',cohort,'_',set,'.csv'))
    write_table(data.frame(human_ensembl=genes$human_ensembl[usable],xx,check.names=FALSE),paste0('gene_logCPM_',cohort,'_',set,'.csv'))
    dif<-vapply(sort(unique(mm$donor_id)),function(d)xx[,mm$donor_id==d&mm$condition=='RSV']-xx[,mm$donor_id==d&mm$condition=='mock'],numeric(nrow(xx)))
    avg<-rowMeans(dif);gene_all[[paste(cohort,set)]]<-data.frame(cohort=cohort,set=set,human_ensembl=genes$human_ensembl[usable],display_symbol=genes$display_symbol[usable],direction_multiplier=sign[usable],mean_paired_logCPM_delta=avg,direction_concordant=avg*sign[usable]>0,n_donors_in_frozen_direction=rowSums(sweep(dif,1,sign[usable],'*')>0))
   }
   objects[[paste(cohort,method,set)]]<-list(logCPM=xx,z=z,score=ss,measured_genes=genes[usable,])
  }
 }
}
scores<-do.call(rbind,scores_all);pairs<-do.call(rbind,pairs_all);coverage<-do.call(rbind,coverage_all);norms<-do.call(rbind,norm_all);gene<-do.call(rbind,gene_all)
rownames(scores)<-rownames(pairs)<-rownames(coverage)<-rownames(norms)<-rownames(gene)<-NULL
write_table(scores,'04_external_sample_scores.csv');write_table(pairs,'05_external_donor_pairs.csv');write_table(coverage,'06_external_member_coverage.csv');write_table(norms,'07_external_normalization.csv');write_table(gene,'08_external_gene_direction.csv')
stat_rows<-list()
for(key in unique(paste(pairs$cohort,pairs$method,pairs$set,sep='|'))){
 pp<-pairs[paste(pairs$cohort,pairs$method,pairs$set,sep='|')==key,];d<-pp$delta;n<-length(d);stopifnot(n==3,!anyDuplicated(pp$donor_id))
 signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));null<-drop(signs%*%d)/n;p<-mean(abs(null)>=abs(mean(d))-1e-12)
 tt<-t.test(d);stat_rows[[key]]<-data.frame(cohort=pp$cohort[1],method=pp$method[1],set=pp$set[1],time_h=pp$time_h[1],n_donors=n,n_positive=sum(d>0),mean_delta=mean(d),CI95_low=tt$conf.int[1],CI95_high=tt$conf.int[2],paired_dz=mean(d)/sd(d),exact_signflip_p=p,paired_t_p=tt$p.value,min_attainable_exact_two_sided_P=2/2^n)
}
statistics<-do.call(rbind,stat_rows);rownames(statistics)<-NULL;statistics$adjusted_p<-NA_real_;statistics$multiplicity_role<-'Normalization sensitivity; descriptive'
i<-statistics$method=='primary'&statistics$set=='priority388';statistics$adjusted_p[i]<-p.adjust(statistics$exact_signflip_p[i],'holm',n=2);statistics$multiplicity_role[i]<-'PRIMARY: Holm across2 declared new independent cohort tests'
i<-statistics$method=='primary'&statistics$set=='broad1980';statistics$adjusted_p[i]<-p.adjust(statistics$exact_signflip_p[i],'BH',n=2);statistics$multiplicity_role[i]<-'SECONDARY: BH across2 declared new independent cohort tests'
write_table(statistics,'09_external_paired_statistics.csv')
summary<-list(completed_at=as.character(Sys.time()),status='Frozen388/1980 external checking; completed cohorts explicitly listed',completed_cohorts=cohorts,pending_cohorts=setdiff(c('GSE270463','GSE309353'),cohorts),results=statistics,coverage=aggregate(cbind(n_frozen=rep(1,nrow(coverage)),present=as.integer(coverage$present),scored=as.integer(coverage$scored))~cohort+method+set,coverage,sum),config_sha256=digest(file=file.path(R,'parameters/analysis_plan_v1.json'),algo='sha256'),implementation_receipt_sha256=digest(file=file.path(R,'parameters/count_implementation_receipt.json'),algo='sha256'),members_changed=FALSE)
write_json(summary,file.path(R,'external_validation_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
saveRDS(list(config=cfg,objects=objects,scores=scores,pairs=pairs,statistics=statistics,coverage=coverage,normalization=norms,gene_direction=gene),file.path(R,'models/external_validation.rds'))
capture.output(sessionInfo(),file=file.path(R,'logs/external_R_sessionInfo.txt'))
file.copy('/tmp/crrm_phase06_score.R',file.path(R,'scripts/03_frozen_external_scoring.R'),overwrite=TRUE)
print(statistics);print(summary$coverage)
