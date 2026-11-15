#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(edgeR);library(jsonlite);library(digest);library(ggplot2)})
P<-file.path(Sys.getenv('CRRM_DATA_ROOT'),"")
REV<-Sys.getenv('CRRM_WORKSPACE')
ROOT<-Sys.getenv('CRRM_PHASE05_VALIDATION_ROOT',unset=file.path(REV,'phase05_mapping_corrected_reconstruction/validation/GSE231666'))
cfg<-fromJSON(file.path(ROOT,'parameters/human_validation_v2.json'))
stopifnot(file.exists(file.path(ROOT,'evidence/source_eligibility_audit.json')))
stopifnot(identical(digest(file=cfg$construction_freeze,algo='sha256'),cfg$construction_freeze_sha256))
for(k in names(cfg$inputs))stopifnot(identical(digest(file=cfg$inputs[[k]]$path,algo='sha256'),cfg$inputs[[k]]$sha256))
write_table<-function(x,name)write.csv(x,file.path(ROOT,'tables',name),row.names=FALSE,na='')
meta<-read.csv(cfg$inputs$manifest$path,check.names=FALSE,stringsAsFactors=FALSE)
meta$complete_mock_RSV_pair<-tolower(as.character(meta$complete_mock_RSV_pair))%in%c('true','t','1')
raw<-read.csv(gzfile(cfg$inputs$counts$path),check.names=FALSE,stringsAsFactors=FALSE)
ensg<-sub('\\.[0-9]+$','',raw$GeneID)
stopifnot(!anyDuplicated(ensg),!anyDuplicated(meta$expression_column),all(meta$expression_column%in%names(raw)))
cm<-as.matrix(raw[,meta$expression_column]);storage.mode(cm)<-'double';rownames(cm)<-ensg
stopifnot(ncol(cm)==47,nrow(cm)==60649,all(is.finite(cm)),all(cm>=0),all(cm==round(cm)))
symbol<-toupper(trimws(raw$gene_name));keep<-!is.na(symbol)&nzchar(symbol)
symbol_cm<-rowsum(cm[keep,,drop=FALSE],group=symbol[keep],reorder=FALSE)
primary_core<-read.csv(cfg$inputs$primary_core$path,check.names=FALSE)
broad_consensus<-read.csv(cfg$inputs$broad_consensus$path,check.names=FALSE)
legacy<-read.csv(cfg$inputs$legacy70$path,check.names=FALSE)
sets<-list(priority_core=list(ids=primary_core$human_ensembl,direction=ifelse(primary_core$frozen_direction=='up',1,-1),type='native'),broad_consensus=list(ids=broad_consensus$human_ensembl,direction=ifelse(broad_consensus$frozen_direction=='up',1,-1),type='native'),legacy70=list(ids=toupper(legacy$gene),direction=rep(1,70),type='symbol'))
stopifnot(nrow(primary_core)==cfg$expected_members$primary_core,nrow(broad_consensus)==cfg$expected_members$broad_consensus,nrow(legacy)==70)
complete<-which(meta$complete_mock_RSV_pair)
stopifnot(length(complete)==44,all(table(meta$chronological_age[complete])==22))
pair_ids<-unique(meta$pair_id[complete])
stopifnot(length(pair_ids)==22)
covlist<-list();normlist<-list();scorelist<-list();pairlist<-list();genesummary<-list();zmatrices<-list();lmats<-list()
make_scores<-function(idx,method){
  dge<-calcNormFactors(DGEList(cm[,idx,drop=FALSE]),method='TMM')
  logcpm<-cpm(dge,log=TRUE,prior.count=1)
  sl<-DGEList(symbol_cm[,idx,drop=FALSE],lib.size=dge$samples$lib.size,norm.factors=dge$samples$norm.factors)
  sym_logcpm<-cpm(sl,log=TRUE,prior.count=1)
  normlist[[method]]<<-data.frame(method=method,gsm=meta$gsm[idx],dge$samples,row.names=NULL)
  for(setname in names(sets)){
    obj<-sets[[setname]];mat<-if(obj$type=='native')logcpm else sym_logcpm
    hit<-match(obj$ids,rownames(mat));available<-!is.na(hit)
    sdev<-rep(NA_real_,length(hit));sdev[available]<-apply(mat[hit[available],,drop=FALSE],1,sd)
    source_counts<-if(obj$type=='native')cm[,idx,drop=FALSE] else symbol_cm[,idx,drop=FALSE]
    nonzero<-rep(NA_integer_,length(hit));nonzero[available]<-rowSums(source_counts[hit[available],,drop=FALSE]>0)
    # All-zero count rows can acquire SD around1e-16 from logCPM arithmetic.
    # They contain no measurable gene variation and must not become amplified z noise.
    usable<-available&is.finite(sdev)&sdev>0&!is.na(nonzero)&nonzero>0
    covlist[[paste(method,setname,sep='__')]]<<-data.frame(method=method,set=setname,frozen_member=obj$ids,frozen_direction=obj$direction,present=available,n_profiles_nonzero=nonzero,sd=sdev,scored=usable)
    stopifnot(any(usable))
    xx<-mat[hit[usable],,drop=FALSE]
    z<-sweep(sweep(xx,1,rowMeans(xx),'-'),1,sdev[usable],'/')
    signed<-sweep(z,1,obj$direction[usable],'*')
    score<-colMeans(signed)
    ss<-meta[idx,];ss$method<-method;ss$set<-setname;ss$score<-unname(score)
    scorelist[[paste(method,setname,sep='__')]]<<-ss
    if(method=='complete44'){
      zm<-data.frame(frozen_member=obj$ids[usable],frozen_direction=obj$direction[usable],z,check.names=FALSE)
      write_table(zm,paste0('gene_zscores_',setname,'_complete44.csv'))
      write_table(data.frame(frozen_member=obj$ids[usable],xx,check.names=FALSE),paste0('gene_logCPM_',setname,'_complete44.csv'))
      zmatrices[[setname]]<<-z;lmats[[setname]]<<-xx
    }
    for(pair in intersect(unique(ss$pair_id),pair_ids)){
      m<-which(ss$pair_id==pair&ss$condition=='mock');r<-which(ss$pair_id==pair&ss$condition=='RSV')
      if(length(m)!=1||length(r)!=1)next
      stopifnot(ss$batch[m]==ss$batch[r],ss$chronological_age[m]==ss$chronological_age[r],ss$child_id_provisional[m]==ss$child_id_provisional[r])
      pr<-data.frame(method=method,set=setname,pair_id=pair,age=ss$chronological_age[m],child_provisional=ss$child_id_provisional[m],raw_origin=ss$donor_sample_origin_raw[m],batch=ss$batch[m],gestational_age=ss$gestational_age[m],mock_gsm=ss$gsm[m],RSV_gsm=ss$gsm[r],score_mock=score[m],score_RSV=score[r],delta=score[r]-score[m])
      pairlist[[paste(method,setname,pair,sep='__')]]<<-pr
    }
  }
}
make_scores(complete,'complete44')
make_scores(seq_len(47),'all47')
for(age in c('Y0','Y1'))make_scores(complete[meta$chronological_age[complete]==age],paste0('within_age_',age))
scores<-do.call(rbind,scorelist);rownames(scores)<-NULL
pairs<-do.call(rbind,pairlist);rownames(pairs)<-NULL
coverage<-do.call(rbind,covlist);rownames(coverage)<-NULL
normalization<-do.call(rbind,normlist);rownames(normalization)<-NULL
write_table(scores,'01_sample_scores.csv');write_table(pairs,'02_pair_scores.csv');write_table(coverage,'03_member_coverage.csv');write_table(normalization,'04_normalization.csv')
signs<-new.env()
exact_flip<-function(delta){
 n<-length(delta);key<-as.character(n)
 if(!exists(key,envir=signs,inherits=FALSE))assign(key,sapply(seq_len(n)-1L,function(j)2L*as.integer(bitwAnd(seq_len(2^n)-1L,bitwShiftL(1L,j))>0L)-1L),envir=signs)
 null<-drop(get(key,envir=signs)%*%delta)/n
 sum(abs(null)>=abs(mean(delta))-1e-12)/length(null)
}
stats<-function(d,flip=TRUE){
 tt<-t.test(d);n<-length(d)
 data.frame(n_pairs=n,n_positive=sum(d>0),mean_delta=mean(d),CI95_low=unname(tt$conf.int[1]),CI95_high=unname(tt$conf.int[2]),median_delta=median(d),paired_dz=mean(d)/sd(d),exact_signflip_p=if(flip)exact_flip(d) else NA_real_,paired_t_p=tt$p.value,wilcoxon_p=suppressWarnings(wilcox.test(d,exact=FALSE)$p.value),sign_p=binom.test(sum(d>0),sum(d!=0))$p.value)
}
out<-list()
for(method in unique(pairs$method))for(set in names(sets))for(age in c('Y0','Y1')){
 x<-pairs[pairs$method==method&pairs$set==set&pairs$age==age,]
 if(!nrow(x))next
 stopifnot(nrow(x)==11,!anyDuplicated(x$child_provisional))
 out[[paste(method,set,age,sep='__')]]<-cbind(data.frame(method=method,set=set,age=age),stats(x$delta))
}
age_stats<-do.call(rbind,out);rownames(age_stats)<-NULL
age_stats$adjusted_p<-NA_real_;age_stats$multiplicity_role<-'Normalization sensitivity, descriptive'
i<-which(age_stats$method=='complete44'&age_stats$set=='priority_core')
age_stats$adjusted_p[i]<-p.adjust(age_stats$exact_signflip_p[i],'holm');age_stats$multiplicity_role[i]<-'PRIMARY: Holm over 2 age tests'
i<-which(age_stats$method=='complete44'&age_stats$set!='priority_core')
age_stats$adjusted_p[i]<-p.adjust(age_stats$exact_signflip_p[i],'BH');age_stats$multiplicity_role[i]<-'SECONDARY: BH over 4 set-age tests'
write_table(age_stats,'05_age_specific_statistics.csv')
child_out<-list();child_values<-list()
for(set in names(sets)){
 pp<-pairs[pairs$method=='complete44'&pairs$set==set,]
 means<-aggregate(delta~child_provisional,pp,mean);means$n_visits<-as.integer(table(pp$child_provisional)[means$child_provisional]);means$set<-set
 child_values[[set]]<-means
 early<-pp[order(pp$age),];early<-early[!duplicated(early$child_provisional),]
 stopifnot(nrow(means)==17,nrow(early)==17)
 child_out[[paste0(set,'_mean')]]<-cbind(data.frame(set=set,repeat_handling='Provisional child mean of observed age deltas'),stats(means$delta))
 child_out[[paste0(set,'_earliest')]]<-cbind(data.frame(set=set,repeat_handling='Earliest observed age per provisional child'),stats(early$delta))
}
child_stats<-do.call(rbind,child_out);rownames(child_stats)<-NULL
write_table(child_stats,'06_provisional_child_sensitivity.csv');write_table(do.call(rbind,child_values),'07_child_mean_deltas.csv')
gene_rows<-list()
for(set in names(sets)){
 obj<-sets[[set]];xx<-lmats[[set]]
 available<-match(rownames(xx),obj$ids);dirs<-obj$direction[available]
 for(age in c('Y0','Y1')){
  pp<-pairs[pairs$method=='complete44'&pairs$set==set&pairs$age==age,]
  dif<-vapply(seq_len(nrow(pp)),function(i)xx[,match(meta$expression_column[match(pp$RSV_gsm[i],meta$gsm)],colnames(xx))]-xx[,match(meta$expression_column[match(pp$mock_gsm[i],meta$gsm)],colnames(xx))],numeric(nrow(xx)))
  dd<-rowMeans(dif)
  gene_rows[[paste(set,age)]]<-data.frame(set=set,age=age,frozen_member=rownames(xx),frozen_direction=dirs,mean_paired_logCPM_delta=dd,direction_concordant=dd*dirs>0,n_pairs_in_frozen_direction=rowSums(sweep(dif,1,dirs,'*')>0))
 }
}
gene_stats<-do.call(rbind,gene_rows);rownames(gene_stats)<-NULL
write_table(gene_stats,'08_gene_direction_concordance.csv')
# Exact historical replay uses symbol-summed libraries, not native Ensembl normalization.
hdge<-calcNormFactors(DGEList(symbol_cm),method='TMM')
hlcpm<-cpm(hdge,log=TRUE,prior.count=1)
hz<-t(scale(t(hlcpm[toupper(legacy$gene),,drop=FALSE])));hscore<-colMeans(hz)
old<-read.csv(file.path(P,'results/Final_Figure2/Figure2B_GSE231666_CRRM_validation/source_tables/Figure2B_sample_metadata_and_CRRM_scores.csv'),check.names=FALSE)
stopifnot(all(meta$expression_column%in%old$sample_id))
cmp<-data.frame(gsm=meta$gsm,sample_id=meta$expression_column,replayed_score=unname(hscore),archived_score=old$CRRM_score[match(meta$expression_column,old$sample_id)])
cmp$absolute_difference<-abs(cmp$replayed_score-cmp$archived_score);stopifnot(max(cmp$absolute_difference)<1e-12)
write_table(cmp,'09_legacy_exact_replay.csv')
hdelta<-vapply(pair_ids,function(pair){ii<-which(meta$pair_id==pair);hscore[ii[meta$condition[ii]=='RSV']]-hscore[ii[meta$condition[ii]=='mock']]},numeric(1))
historical_summary<-stats(hdelta,flip=FALSE)
write_table(historical_summary,'10_legacy_replay_statistics.csv')
summary<-list(completed_at=as.character(Sys.time()),status='Frozen candidate human external checking, discovery qualifications retained',n_raw_profiles=47,n_complete_profiles=44,n_pairs=22,n_pairs_per_age=11,n_provisional_children=17,primary_results=age_stats[age_stats$method=='complete44'&age_stats$set=='priority_core',],secondary_results=age_stats[age_stats$method=='complete44'&age_stats$set!='priority_core',],historical_max_abs_score_difference=max(cmp$absolute_difference),historical_replay=historical_summary,config_sha256=digest(file=file.path(ROOT,'parameters/human_validation_v2.json'),algo='sha256'),frozen_members_unchanged=TRUE,clinical_or_single_cell_run=FALSE)
write_json(summary,file.path(ROOT,'validation_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
saveRDS(list(config=cfg,scores=scores,pairs=pairs,age_statistics=age_stats,coverage=coverage,normalization=normalization,priority_core_logCPM=lmats$priority_core,priority_core_z=zmatrices$priority_core,primary_genes=primary_core,broad_genes=broad_consensus),file.path(ROOT,'models/GSE231666_frozen_validation.rds'))
capture.output(sessionInfo(),file=file.path(ROOT,'logs/R_sessionInfo.txt'))
print(age_stats[age_stats$method=='complete44',c('set','age','mean_delta','CI95_low','CI95_high','exact_signflip_p','adjusted_p','n_positive')])
print(paste('Historical score max deviation',max(cmp$absolute_difference)))
