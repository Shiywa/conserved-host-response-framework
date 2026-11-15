suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)]);out<-file.path(Sys.getenv('CRRM_WORKSPACE'),'phase18_Result2_external_response_20261112')
objects<-readRDS(file.path(out,'data/processed/external_log_expression.rds'))
coverage<-read.csv(file.path(out,'tables/02_member_coverage.csv'))
norm<-read.csv(file.path(out,'tables/03_normalization_and_projection.csv'))
dir.create(file.path(out,'evidence/before_numerical_variance_guard'),showWarnings=FALSE)
for(f in c('02_member_coverage.csv','03_normalization_and_projection.csv','04_coverage_summary.csv','21_additional_context_effects.csv')){
 p<-file.path(out,'tables',f);b<-file.path(out,'evidence/before_numerical_variance_guard',f);if(file.exists(p)&&!file.exists(b))file.copy(p,b)
}
audit<-list()
for(name in names(objects)){
 o<-objects[[name]];sd0<-apply(o$x[,o$reference,drop=FALSE],1,sd);bad<-names(sd0)[is.finite(sd0)&sd0>0&sd0<=1e-8]
 if(length(bad))audit[[name]]<-data.frame(dataset=name,human_ensembl=bad,reference_sd=sd0[bad],reason='Machine-precision variance; not suitable for z-score standardization')
 ix<-which(coverage$dataset==name);v<-sd0[match(coverage$human_ensembl[ix],names(sd0))];coverage$scored[ix]<-is.finite(v)&v>1e-8;norm$n_scored[norm$dataset==name]<-sum(coverage$scored[ix])
}
audit<-do.call(rbind,audit);stopifnot(all(audit$dataset=='PRJNA1271255_Blood'))
write.csv(audit,file.path(out,'tables/22_numerical_variance_exclusions.csv'),row.names=FALSE)
write.csv(coverage,file.path(out,'tables/02_member_coverage.csv'),row.names=FALSE,na='')
write.csv(norm,file.path(out,'tables/03_normalization_and_projection.csv'),row.names=FALSE,na='')
summary<-do.call(rbind,lapply(split(coverage,coverage$dataset),function(c){data.frame(dataset=c$dataset[1],n_present=sum(c$present),n_scored=sum(c$scored),n_up=sum(c$scored&c$frozen_direction=='up'),n_down=sum(c$scored&c$frozen_direction=='down'),eligible_primary=sum(c$scored)/388>=.7&&sum(c$scored&c$frozen_direction=='up')/280>=.5&&sum(c$scored&c$frozen_direction=='down')/108>=.5)}))
write.csv(summary,file.path(out,'tables/04_coverage_summary.csv'),row.names=FALSE)
write_json(list(amendment='Numerical variance guard',threshold_reference_SD=1e-8,scale='Native normalized log expression',declared_after_initial_scores=TRUE,reason='Independent R/Python reconstruction detected numerical instability in machine-precision SDs of zero-count hamster blood genes',excluded_count=nrow(audit),affected_datasets=unique(audit$dataset),primary_contrasts_affected=0,member_list_changed=FALSE,outcome_based_selection=FALSE,eligibility_gates_unchanged=TRUE,replay_order='01_prepare;01b_variance_guard;02_export;03_analysis;04_checks'),file.path(out,'parameters/numerical_variance_amendment.json'),pretty=TRUE,auto_unbox=TRUE)
print(summary);print(audit)
