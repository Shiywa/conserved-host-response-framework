suppressPackageStartupMessages({library(jsonlite)})
args<-commandArgs(trailingOnly=FALSE);f<-sub('^--file=','',args[grepl('^--file=',args)]);OUT<-normalizePath(file.path(dirname(f),'..'))
effects<-read.csv(file.path(OUT,'results/virus_comparison_effects.csv'),check.names=FALSE);values<-read.csv(file.path(OUT,'results/participant_donor_values.csv'),check.names=FALSE)
delta<-list();checks<-list()
for(i in seq_len(nrow(effects))){
  e<-effects[i,];v<-values[values$key==e$key&values$score==e$score,];a<-v[v$role=='case',];b<-v[v$role=='control',]
  paired<-tolower(as.character(e$paired))=='true'
  if(paired){b<-b[match(a$analysis_unit,b$analysis_unit),];stopifnot(identical(a$analysis_unit,b$analysis_unit))}
  aa<-a$value[is.finite(a$value)];bb<-b$value[is.finite(b$value)]
  if(min(length(aa),length(bb))<3||!is.finite(e$P))next
  t<-t.test(aa,bb,paired=paired,var.equal=FALSE);est<-mean(aa)-mean(bb)
  err<-max(abs(c(est-e$estimate,t$conf.int[1]-e$lower,t$conf.int[2]-e$upper,t$p.value-e$P)))
  checks[[length(checks)+1]]<-data.frame(key=e$key,score=e$score,max_error=err);stopifnot(err<1e-10)
}
qerrs<-c()
for(family in unique(effects$correction_family)){
  ix<-which(effects$correction_family==family);p<-effects$P[ix];p[!is.finite(p)]<-1;q<-p.adjust(p,'BH');ok<-is.finite(effects$adjusted_P[ix]);err<-max(abs(q[ok]-effects$adjusted_P[ix][ok]));qerrs<-c(qerrs,err);stopifnot(err<1e-10)
}
core<-read.csv(file.path(OUT,'parameters/frozen_CRRM388.csv'));refs<-read.csv(file.path(OUT,'data/gene_reference_parameters.csv'));arr<-readRDS(file.path(OUT,'data/array_expression.rds'));gp<-read.csv(file.path(OUT,'results/profile_scores.csv'),check.names=FALSE)
scoreerrs<-c()
for(source in names(arr)){
  a<-arr[[source]];r<-refs[refs$source==source&tolower(as.character(refs$variable))=='true',];x<-a$x[r$human_ensembl,,drop=FALSE];sign<-ifelse(core$frozen_direction[match(r$human_ensembl,core$human_ensembl)]=='up',1,-1);z<-(x-r$mean)/r$SD;total<-colMeans(z*sign);g<-gp[gp$source_label==source,];g<-g[match(colnames(x),g$geo_accession),];err<-max(abs(total-g$CRRM388));scoreerrs<-c(scoreerrs,err);stopifnot(err<1e-10)
}
write.csv(do.call(rbind,checks),file.path(OUT,'audit/independent_R_statistics.csv'),row.names=FALSE)
genechecks<-list()
for(source in names(arr)){
  a<-arr[[source]];r<-refs[refs$source==source,];r<-r[match(core$human_ensembl,r$human_ensembl),];
  file<-if(source=='GSE42026')'replication_all_gene_differences.csv'else'discovery_gene_differences.csv'
  recorded<-read.csv(file.path(OUT,'results',file),check.names=FALSE)
  for(i in seq_len(nrow(recorded))){
    z<-recorded[i,];id<-z$human_ensembl;j<-match(id,r$human_ensembl)
    if(!is.finite(z$P)||!id%in%rownames(a$x))next
    sg<-ifelse(core$frozen_direction[match(id,core$human_ensembl)]=='up',1,-1)
    y<-(a$x[id,]-r$mean[j])/r$SD[j]*sg;aa<-y[a$meta$condition==z$case];bb<-y[a$meta$condition==z$control]
    tt<-tryCatch(t.test(aa,bb,var.equal=FALSE),error=function(e)NULL);if(is.null(tt))stop('FinitePython gene P but independentR comparison non-estimable: ',id,' ',z$key)
    err<-max(abs(c(mean(aa)-mean(bb)-z$aligned_gene_effect,tt$p.value-z$P,tt$conf.int[1]-z$lower,tt$conf.int[2]-z$upper)));stopifnot(err<1e-9)
    genechecks[[length(genechecks)+1]]<-data.frame(source=source,gene=id,key=z$key,error=err)
  }
  p<-recorded$P;p[!is.finite(p)]<-1;qq<-p.adjust(p,'BH');ok<-is.finite(recorded$adjusted_P);stopifnot(max(abs(qq[ok]-recorded$adjusted_P[ok]))<1e-10)
}
write.csv(do.call(rbind,genechecks),file.path(OUT,'audit/independent_R_gene_tests.csv'),row.names=FALSE)
write_json(list(status='PASS',tests_checked=length(checks),max_statistic_error=max(do.call(rbind,checks)$max_error),BH_families_checked=length(qerrs),max_BH_error=max(qerrs),array_sources_reconstructed=length(scoreerrs),max_array_score_error=max(scoreerrs),gene_tests_checked=length(genechecks),max_gene_test_error=max(do.call(rbind,genechecks)$error),gene_BH_families_checked=c(776,388)),file.path(OUT,'audit/independent_R_receipt.json'),pretty=TRUE,auto_unbox=TRUE)
message('Independent R replay passed: ',length(checks),'tests')
