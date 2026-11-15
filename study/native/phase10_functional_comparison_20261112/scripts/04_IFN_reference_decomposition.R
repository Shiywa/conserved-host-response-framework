#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(jsonlite))
OUT<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]))))
REV<-file.path(Sys.getenv('CRRM_WORKSPACE'),"")
a<-readRDS(file.path(OUT,'data/processed/functional_atlas.rds'));q<-a$members;atlas<-a$atlas
membership<-list(CRRM388=q$human_ensembl,IFN_overlap=atlas$human_ensembl[atlas$Hallmark_IFN_union],Outside_IFN_reference=atlas$human_ensembl[!atlas$Hallmark_IFN_union])
y<-readRDS(file.path(REV,'phase05_mapping_corrected_reconstruction/validation/GSE231666/models/GSE231666_frozen_validation.rds'))
ss<-y$scores[y$scores$method=='complete44'&y$scores$set=='priority_core',];ss<-ss[match(colnames(y$priority_core_z),ss$expression_column),]
objects<-list(GSE231666=list(z=y$priority_core_z,metadata=data.frame(cohort='GSE231666',stratum=ss$chronological_age,unit_id=ss$pair_id,condition=ss$condition,profile_id=ss$expression_column,time_h=96),old_score=ss$score))
x<-readRDS(file.path(REV,'phase06_external_validation_figure1/models/external_validation.rds'))
for(cohort in c('GSE270463','GSE309353')){
 ob<-x$objects[[paste(cohort,'primary','priority388')]];m<-ob$score
 stopifnot(identical(colnames(ob$z),as.character(m$profile_id)))
 objects[[cohort]]<-list(z=ob$z,metadata=data.frame(cohort=cohort,stratum='All',unit_id=m$donor_id,condition=m$condition,profile_id=m$profile_id,time_h=m$time_h),old_score=m$score)
}
scores<-list();pairs<-list();stats<-list();checks<-list()
for(cohort in names(objects)){
 ob<-objects[[cohort]];z<-ob$z;mult<-ifelse(q$frozen_direction[match(rownames(z),q$human_ensembl)]=='up',1,-1);signed<-sweep(z,1,mult,'*')
 checks[[cohort]]<-data.frame(cohort=cohort,max_abs_replayed_total_score_error=max(abs(colMeans(signed)-ob$old_score)))
 stopifnot(max(abs(colMeans(signed)-ob$old_score))<1e-12)
 comp<-list()
 for(set in names(membership)){
  idx<-rownames(signed)%in%membership[[set]];stopifnot(sum(idx)>0);score<-colMeans(signed[idx,,drop=FALSE]);comp[[set]]<-list(score=score,n=sum(idx))
  mm<-ob$metadata;mm$component<-set;mm$n_measured<-sum(idx);mm$n_fixed<-length(membership[[set]]);mm$score<-unname(score);scores[[paste(cohort,set)]]<-mm
  pp<-do.call(rbind,lapply(unique(mm$unit_id),function(id){s<-mm[mm$unit_id==id,];stopifnot(nrow(s)==2);data.frame(cohort=cohort,stratum=unique(s$stratum),unit_id=id,component=set,n_measured=sum(idx),n_fixed=length(membership[[set]]),delta=s$score[s$condition=='RSV']-s$score[s$condition=='mock'])}));pairs[[paste(cohort,set)]]<-pp
  for(stratum in unique(pp$stratum)){
   v<-pp$delta[pp$stratum==stratum];n<-length(v);se<-sd(v)/sqrt(n);ci<-mean(v)+c(-1,1)*qt(.975,n-1)*se
   # Descriptive follow-up; exact P retained with BH over the 12 component-by-stratum tests.
   signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));pv<-mean(abs(drop(signs%*%v)/n)>=abs(mean(v))-1e-12)
   stats[[paste(cohort,set,stratum)]]<-data.frame(cohort=cohort,stratum=stratum,component=set,n_units=n,n_measured=sum(idx),n_fixed=length(membership[[set]]),n_positive=sum(v>0),mean_delta=mean(v),CI95_low=ci[1],CI95_high=ci[2],exact_signflip_p=pv)
  }
 }
 reconstructed<-(comp$IFN_overlap$score*comp$IFN_overlap$n+comp$Outside_IFN_reference$score*comp$Outside_IFN_reference$n)/comp$CRRM388$n
 stopifnot(max(abs(reconstructed-comp$CRRM388$score))<1e-12)
 checks[[cohort]]$max_abs_component_reconstruction_error<-max(abs(reconstructed-comp$CRRM388$score))
}
sco<-do.call(rbind,scores);pp<-do.call(rbind,pairs);st<-do.call(rbind,stats);st$FDR_secondary_BH<-p.adjust(st$exact_signflip_p,'BH')
write.csv(sco,file.path(OUT,'tables/13_IFN_component_sample_scores.csv'),row.names=FALSE);write.csv(pp,file.path(OUT,'tables/14_IFN_component_paired_deltas.csv'),row.names=FALSE);write.csv(st,file.path(OUT,'tables/15_IFN_component_summary.csv'),row.names=FALSE);write.csv(do.call(rbind,checks),file.path(OUT,'tables/16_score_replay_checks.csv'),row.names=FALSE)
write_json(list(IFN_reference='Exact Hallmark alpha/gamma union; not an exhaustive IFN-dependent gene inventory',n_IFN_overlap=79,n_outside_IFN_reference=309,followup='Descriptive post-freeze decomposition; same precomputed gene z scores and frozen signs; no refitting or gene tuning; cohort-relative score magnitudes are not cross-cohort calibration',results=st),file.path(OUT,'IFN_decomposition_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
print(st);print(do.call(rbind,checks))
