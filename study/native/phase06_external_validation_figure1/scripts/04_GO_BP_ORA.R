#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db);library(GO.db);library(jsonlite);library(digest)})
R<-file.path(Sys.getenv('CRRM_WORKSPACE'),"phase06_external_validation_figure1")
cfg<-fromJSON(file.path(R,'parameters/analysis_plan_v1.json'))
for(k in names(cfg$inputs))stopifnot(digest(file=cfg$inputs[[k]]$path,algo='sha256')==cfg$inputs[[k]]$sha256)
u<-read.csv(cfg$inputs$universe$path,check.names=FALSE);q<-read.csv(cfg$inputs$priority388$path,check.names=FALSE)
flag<-function(x)tolower(as.character(x))%in%c('true','t','1')
u<-u[flag(u$A549_mapped_and_tested)&flag(u$mouse_mapped_and_tested)&flag(u$cotton_rat_mapped_and_tested),]
stopifnot(all(q$human_ensembl%in%u$human_ensembl),!anyDuplicated(u$human_ensembl))
valid<-intersect(u$human_ensembl,keys(org.Hs.eg.db,keytype='ENSEMBL'))
ann<-AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='ENSEMBL',columns=c('GOALL','ONTOLOGYALL'))
ann<-unique(ann[!is.na(ann$GOALL)&ann$ONTOLOGYALL=='BP',c('ENSEMBL','GOALL')])
go<-AnnotationDbi::select(GO.db,keys=unique(ann$GOALL),keytype='GOID',columns=c('TERM','ONTOLOGY'))
go<-go[!is.na(go$TERM)&go$ONTOLOGY=='BP',];ann<-ann[ann$GOALL%in%go$GOID,]
background<-sort(unique(ann$ENSEMBL));N<-length(background)
bg<-data.frame(human_ensembl=u$human_ensembl,GO_BP_annotated=u$human_ensembl%in%background,in_joint_tested_background=TRUE)
write.csv(bg,file.path(R,'tables/15_GO_background_audit.csv'),row.names=FALSE)
edgecon<-gzfile(file.path(R,'data/processed/GO_BP_background_edges.csv.gz'),'wt');write.csv(ann,edgecon,row.names=FALSE);close(edgecon)
sets<-split(ann$ENSEMBL,ann$GOALL);sets<-lapply(sets,unique);sets<-sets[lengths(sets)>=10&lengths(sets)<=500]
out<-list()
for(direction in c('up','down')){
 query<-intersect(q$human_ensembl[q$frozen_direction==direction],background);n<-length(query)
 out[[direction]]<-do.call(rbind,lapply(names(sets),function(term){
  genes<-sets[[term]];hits<-intersect(query,genes);K<-length(genes);k<-length(hits)
  data.frame(direction=direction,GOID=term,term=go$TERM[match(term,go$GOID)],N_background=N,n_query=n,K_term=K,k_overlap=k,p_value=phyper(k-1,K,N-K,n,lower.tail=FALSE),fold_enrichment=if(n)k/n/(K/N)else NA_real_,overlap_genes=paste(sort(hits),collapse=';'))
 }))
}
all<-do.call(rbind,out);rownames(all)<-NULL;all$FDR_joint_BH<-p.adjust(all$p_value,'BH');all<-all[order(all$FDR_joint_BH,all$p_value,-all$k_overlap,all$GOID),]
write.csv(all,file.path(R,'tables/16_GO_BP_ORA_full.csv'),row.names=FALSE)
selected<-list()
for(direction in c('up','down')){
 x<-all[all$direction==direction&all$FDR_joint_BH<.05,];picked<-character()
 if(nrow(x))for(i in seq_len(nrow(x))){
  id<-x$GOID[i];g<-sets[[id]]
  too_similar<-length(picked)>0&&any(vapply(picked,function(p)length(intersect(g,sets[[p]]))/length(union(g,sets[[p]]))>0.5,logical(1)))
  if(!too_similar)picked<-c(picked,id)
  if(length(picked)==5)break
 }
 selected[[direction]]<-x[x$GOID%in%picked,]
}
display<-do.call(rbind,selected);write.csv(display,file.path(R,'tables/17_GO_display_terms.csv'),row.names=FALSE)
write.csv(data.frame(human_ensembl=q$human_ensembl,display_symbol=q$display_symbol,direction=q$frozen_direction,in_GO_background=q$human_ensembl%in%background),file.path(R,'tables/18_priority_GO_coverage.csv'),row.names=FALSE)
metadata_out<-list(org_Hs_eg_db=as.character(packageVersion('org.Hs.eg.db')),GO_db=as.character(packageVersion('GO.db')),org_metadata=metadata(org.Hs.eg.db),GO_metadata=metadata(GO.db),source_org_sqlite_sha256=digest(file=system.file('extdata','org.Hs.eg.sqlite',package='org.Hs.eg.db'),algo='sha256'),source_GO_sqlite_sha256=digest(file=system.file('extdata','GO.sqlite',package='GO.db'),algo='sha256'))
write_json(metadata_out,file.path(R,'evidence/GO_annotation_provenance.json'),pretty=TRUE,auto_unbox=TRUE)
summary<-list(n_joint_tested_genes=nrow(u),n_GO_background=N,n_eligible_terms=length(sets),n_tests_joint_BH=nrow(all),n_up_query=length(intersect(q$human_ensembl[q$frozen_direction=='up'],background)),n_down_query=length(intersect(q$human_ensembl[q$frozen_direction=='down'],background)),n_significant_up=sum(all$direction=='up'&all$FDR_joint_BH<.05),n_significant_down=sum(all$direction=='down'&all$FDR_joint_BH<.05),display_terms=display,selection_uses_external_response=FALSE)
write_json(summary,file.path(R,'GO_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
saveRDS(list(config=cfg,background=background,annotation=ann,term_sets=sets,results=all,display=display),file.path(R,'models/GO_BP_ORA.rds'))
capture.output(sessionInfo(),file=file.path(R,'logs/GO_R_sessionInfo.txt'))
file.copy('/tmp/crrm_phase06_GO.R',file.path(R,'scripts/04_GO_BP_ORA.R'),overwrite=TRUE)
print(summary[c('n_joint_tested_genes','n_GO_background','n_eligible_terms','n_tests_joint_BH','n_up_query','n_down_query','n_significant_up','n_significant_down')]);print(display[,c('direction','term','k_overlap','FDR_joint_BH')])
