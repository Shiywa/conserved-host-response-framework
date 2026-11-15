a<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',a[grepl('^--file=',a)]);OUT<-file.path(Sys.getenv('CRRM_WORKSPACE'),'Result4_analysis')
effects<-read.csv(file.path(OUT,'results/primary_cellular_effects.csv'),check.names=FALSE,stringsAsFactors=FALSE)
units<-read.csv(file.path(OUT,'results/primary_unit_effects.csv'),check.names=FALSE,stringsAsFactors=FALSE)
records<-list()
for(i in seq_len(nrow(effects))){
 r<-effects[i,];d<-units[units$comparison_id==r$comparison_id&units$score==r$score,]
 if(r$status!='OK')next
 if(r$design=='unpaired_participants')test<-t.test(d$case_score[!is.na(d$case_score)],d$control_score[!is.na(d$control_score)],var.equal=FALSE)
 else test<-t.test(d$case_score,d$control_score,paired=TRUE)
 estimate<-if(length(test$estimate)==2)unname(test$estimate[1]-test$estimate[2])else unname(test$estimate)
 observed<-c(estimate,test$conf.int[1],test$conf.int[2],test$p.value);expected<-c(r$estimate,r$lower,r$upper,r$P)
 stopifnot(max(abs(observed-expected))<1e-10)
 records[[length(records)+1]]<-data.frame(comparison=r$comparison_id,score=r$score,estimate_abs_error=abs(observed[1]-expected[1]),CI_max_abs_error=max(abs(observed[2:3]-expected[2:3])),P_abs_error=abs(observed[4]-expected[4]),status='PASS')
}
for(family in unique(effects$correction_family)){
 ix<-which(effects$correction_family==family);q<-p.adjust(effects$P[ix],method='BH',n=length(ix));stopifnot(isTRUE(all.equal(q,effects$adjusted_P[ix],tolerance=1e-12,check.attributes=FALSE)))
}
result<-do.call(rbind,records);write.csv(result,file.path(OUT,'audit/independent_R_primary_statistics.csv'),row.names=FALSE)
cat('PASS:',nrow(result),'independently recomputed effects/intervals/P values; all BH families agree\n')
