#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db);library(GO.db);library(jsonlite);library(digest)})
args<-commandArgs(trailingOnly=FALSE);OUT<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
cfg<-fromJSON(file.path(OUT,'parameters/functional_comparison_plan_v1.json'))
for(k in names(cfg$inputs))stopifnot(digest(file=cfg$inputs[[k]]$path,algo='sha256')==cfg$inputs[[k]]$sha256)
snapshots<-fromJSON(file.path(OUT,'evidence/public_source_downloads.json'))
for(name in c('hallmark_2026_1.gmt','Reactome_pathways.gmt.zip','Reactome_release.json')){
 row<-snapshots[snapshots$name==name,];stopifnot(nrow(row)==1,digest(file=file.path(OUT,'data/public',name),algo='sha256')==row$sha256)
}
stopifnot(trimws(paste(readLines(file.path(OUT,'data/public/Reactome_release.json'),warn=FALSE),collapse=''))=='98')
un<-read.csv(cfg$inputs$universe$path,check.names=FALSE);q<-read.csv(cfg$inputs$priority388$path,check.names=FALSE)
flag<-function(x)tolower(as.character(x))%in%c('true','t','1')
u<-un[flag(un$A549_mapped_and_tested)&flag(un$mouse_mapped_and_tested)&flag(un$cotton_rat_mapped_and_tested),]
stopifnot(nrow(u)==8607,nrow(q)==388,all(q$human_ensembl%in%u$human_ensembl),!anyDuplicated(u$human_ensembl))
valid<-intersect(u$human_ensembl,keys(org.Hs.eg.db,keytype='ENSEMBL'))
ann<-AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='ENSEMBL',columns=c('SYMBOL','ENTREZID','GENENAME','GOALL','ONTOLOGYALL'))
go<-AnnotationDbi::select(GO.db,keys=unique(na.omit(ann$GOALL)),keytype='GOID',columns=c('TERM','ONTOLOGY'))
edges<-list();terms<-list();sets<-list();backgrounds<-list();mapping_audit<-list()
for(ont in c('BP','MF','CC')){
 lib<-paste0('GO_',ont);a<-unique(ann[!is.na(ann$GOALL)&ann$ONTOLOGYALL==ont,c('ENSEMBL','GOALL')]);a<-a[a$GOALL%in%go$GOID,];names(a)<-c('human_ensembl','term_id')
 edges[[lib]]<-data.frame(library=lib,a,stringsAsFactors=FALSE)
 terms[[lib]]<-unique(data.frame(library=lib,term_id=a$term_id,term=go$TERM[match(a$term_id,go$GOID)],stringsAsFactors=FALSE))
}
read_gmt<-function(path,lib){
 x<-strsplit(readLines(path,warn=FALSE),'\t',fixed=TRUE);out<-list();meta<-list()
 # Exact frozen display-symbol matching; ambiguous symbols do not silently duplicate identifiers.
 su<-split(u$human_ensembl,u$display_symbol);su<-su[lengths(su)==1]
 for(i in seq_along(x)){
  z<-x[[i]];id<-if(lib=='Reactome98')z[2]else z[1];sy<-unique(z[-c(1,2)]);mapped<-unlist(su[intersect(sy,names(su))],use.names=FALSE)
  if(length(mapped))out[[i]]<-data.frame(library=lib,human_ensembl=unique(mapped),term_id=id)
  meta[[i]]<-data.frame(library=lib,term_id=id,term=z[1],n_source_symbols=length(sy),n_mapped_common=length(unique(mapped)),n_unmapped_symbols=sum(!sy%in%names(su)))
 }
 edges[[lib]]<<-unique(do.call(rbind,out));terms[[lib]]<<-unique(do.call(rbind,meta))[,c('library','term_id','term')]
 mapping_audit[[lib]]<<-do.call(rbind,meta)
}
read_gmt(file.path(OUT,'data/public/ReactomePathways.gmt'),'Reactome98')
read_gmt(file.path(OUT,'data/public/hallmark_2026_1.gmt'),'Hallmark2026.1.Hs')
mapping<-do.call(rbind,mapping_audit);write.csv(mapping,file.path(OUT,'tables/02_reference_mapping_audit.csv'),row.names=FALSE)
all_edges<-do.call(rbind,edges);all_terms<-unique(do.call(rbind,terms));rownames(all_edges)<-NULL
con<-gzfile(file.path(OUT,'data/processed/all_background_annotation_edges.csv.gz'),'wt');write.csv(all_edges,con,row.names=FALSE);close(con)
write.csv(all_terms,file.path(OUT,'tables/03_term_dictionary.csv'),row.names=FALSE)
res<-list();coverage<-list();eligible<-list()
for(lib in names(edges)){
 a<-edges[[lib]];bg<-sort(unique(a$human_ensembl));backgrounds[[lib]]<-bg;ss<-lapply(split(a$human_ensembl,a$term_id),unique);sets[[lib]]<-ss;ss<-ss[lengths(ss)>=10&lengths(ss)<=500];eligible[[lib]]<-names(ss)
 coverage[[lib]]<-data.frame(library=lib,background_common=8607,background_annotated=length(bg),members_annotated=sum(q$human_ensembl%in%bg),up_annotated=sum(q$human_ensembl[q$frozen_direction=='up']%in%bg),down_annotated=sum(q$human_ensembl[q$frozen_direction=='down']%in%bg),n_eligible_terms=length(ss))
 for(dr in c('up','down')){
  query_all<-q$human_ensembl[q$frozen_direction==dr];query<-intersect(query_all,bg);N<-length(bg);n<-length(query)
  rr<-do.call(rbind,lapply(names(ss),function(id){
   hit<-intersect(query,ss[[id]]);K<-length(ss[[id]]);k<-length(hit)
   data.frame(library=lib,direction=dr,term_id=id,term=all_terms$term[match(paste(lib,id),paste(all_terms$library,all_terms$term_id))],N_background=N,n_query=n,K_term=K,k_overlap=k,
              p_value=phyper(k-1,K,N-K,n,lower.tail=FALSE),fold_enrichment=k/n/(K/N),
              p_full8607=phyper(k-1,K,8607-K,length(query_all),lower.tail=FALSE),
              overlap_ensembl=paste(sort(hit),collapse=';'),overlap_symbols=paste(sort(u$display_symbol[match(hit,u$human_ensembl)]),collapse=';'))
  }));res[[paste(lib,dr)]]<-rr
 }
}
r<-do.call(rbind,res);rownames(r)<-NULL;r$FDR_global_BH<-p.adjust(r$p_value,'BH');r$FDR_full8607_global_BH<-p.adjust(r$p_full8607,'BH');r$FDR_library_BH<-ave(r$p_value,r$library,FUN=function(z)p.adjust(z,'BH'));r<-r[order(r$FDR_global_BH,r$p_value,-r$k_overlap,r$library,r$term_id),]
write.csv(r,file.path(OUT,'tables/04_functional_ORA_all_terms.csv'),row.names=FALSE)
cv<-do.call(rbind,coverage);cv$significant_up<-vapply(cv$library,function(l)sum(r$library==l&r$direction=='up'&r$FDR_global_BH<.05),integer(1));cv$significant_down<-vapply(cv$library,function(l)sum(r$library==l&r$direction=='down'&r$FDR_global_BH<.05),integer(1));write.csv(cv,file.path(OUT,'tables/05_annotation_coverage.csv'),row.names=FALSE)
# Complete member-level annotation table, and a compact all388 atlas.
qe<-all_edges[all_edges$human_ensembl%in%q$human_ensembl,];qe$term<-all_terms$term[match(paste(qe$library,qe$term_id),paste(all_terms$library,all_terms$term_id))];qe$display_symbol<-q$display_symbol[match(qe$human_ensembl,q$human_ensembl)];qe$direction<-q$frozen_direction[match(qe$human_ensembl,q$human_ensembl)];write.csv(qe,file.path(OUT,'tables/06_all388_annotation_edges.csv'),row.names=FALSE)
atlas<-u[match(q$human_ensembl,u$human_ensembl),];atlas$frozen_direction<-q$frozen_direction;atlas$direction_multiplier<-q$direction_multiplier
uniq_join<-function(z)paste(sort(unique(na.omit(z))),collapse=';')
atlas$gene_name<-vapply(atlas$human_ensembl,function(id)uniq_join(ann$GENENAME[ann$ENSEMBL==id]),character(1));atlas$EntrezID<-vapply(atlas$human_ensembl,function(id)uniq_join(ann$ENTREZID[ann$ENSEMBL==id]),character(1))
for(lib in names(edges)){
 a<-qe[qe$library==lib,];atlas[[paste0(lib,'_term_count')]]<-vapply(atlas$human_ensembl,function(id)sum(a$human_ensembl==id),integer(1))
 atlas[[paste0(lib,'_terms')]]<-vapply(atlas$human_ensembl,function(id)uniq_join(paste0(a$term_id[a$human_ensembl==id],': ',a$term[a$human_ensembl==id])),character(1))
}
ifna<-sets[['Hallmark2026.1.Hs']][['HALLMARK_INTERFERON_ALPHA_RESPONSE']];ifng<-sets[['Hallmark2026.1.Hs']][['HALLMARK_INTERFERON_GAMMA_RESPONSE']];ifn<-union(ifna,ifng)
atlas$Hallmark_IFNA<-atlas$human_ensembl%in%ifna;atlas$Hallmark_IFNG<-atlas$human_ensembl%in%ifng;atlas$Hallmark_IFN_union<-atlas$human_ensembl%in%ifn
write.csv(atlas,file.path(OUT,'tables/01_all388_gene_function_atlas.csv'),row.names=FALSE)
hm<-mapping[mapping$library=='Hallmark2026.1.Hs'&grepl('INTERFERON_(ALPHA|GAMMA)_RESPONSE',mapping$term_id),];hm$n_CRRM<-vapply(hm$term_id,function(id)sum(q$human_ensembl%in%sets[['Hallmark2026.1.Hs']][[id]]),integer(1));write.csv(hm,file.path(OUT,'tables/07_IFN_reference_overlap.csv'),row.names=FALSE)
# Redundancy reduction for the display only. These representative terms are not inferred biological modules.
picked<-list()
for(lib in names(edges))for(dr in c('up','down')){
 x<-r[r$library==lib&r$direction==dr&r$FDR_global_BH<.05,];ids<-character();lim<-if(lib=='GO_BP')8 else if(lib=='Reactome98')6 else 4
 if(nrow(x))for(i in seq_len(nrow(x))){id<-x$term_id[i];g<-sets[[lib]][[id]];similar<-length(ids)>0&&any(vapply(ids,function(j)length(intersect(g,sets[[lib]][[j]]))/length(union(g,sets[[lib]][[j]]))>.5,logical(1)));if(!similar)ids<-c(ids,id);if(length(ids)>=lim)break}
 picked[[paste(lib,dr)]]<-x[x$term_id%in%ids,]
}
display<-do.call(rbind,picked);write.csv(display,file.path(OUT,'tables/08_representative_functional_terms.csv'),row.names=FALSE)
# Publicly described Wang2026 biology: identify matching pathway definitions without creating an author signature.
focus<-r[grepl('interleukin.?17|IL.?17|Toll.?like|NF.?kappa|NF.?KB|neutrophil degranulation|ZBP1|Z.?DNA|type I interferon|interferon.alpha|interferon.gamma',r$term,ignore.case=TRUE),];write.csv(focus,file.path(OUT,'tables/09_Wang2026_related_pathway_results.csv'),row.names=FALSE)
metadata_out<-list(org_Hs_eg_db=as.character(packageVersion('org.Hs.eg.db')),GO_db=as.character(packageVersion('GO.db')),org_metadata=metadata(org.Hs.eg.db),GO_metadata=metadata(GO.db),org_sqlite_sha256=digest(file=system.file('extdata','org.Hs.eg.sqlite',package='org.Hs.eg.db'),algo='sha256'),GO_sqlite_sha256=digest(file=system.file('extdata','GO.sqlite',package='GO.db'),algo='sha256'),Reactome_release=readLines(file.path(OUT,'data/public/Reactome_release.json')),Hallmark_release='2026.1.Hs',symbol_mapping='Exact frozen display-symbol match to unique symbol-to-Ensembl links within joint8607; all absent mappings recorded')
write_json(metadata_out,file.path(OUT,'evidence/annotation_provenance.json'),pretty=TRUE,auto_unbox=TRUE)
summary<-list(n_members=388,n_up=280,n_down=108,n_common_background=8607,n_primary_tests=nrow(r),primary_FDR='Global BH over all five libraries and both directions',coverage=cv,IFN_union_members=sum(atlas$Hallmark_IFN_union),IFN_union_up=sum(atlas$Hallmark_IFN_union&atlas$frozen_direction=='up'),IFN_union_down=sum(atlas$Hallmark_IFN_union&atlas$frozen_direction=='down'),outside_IFN_union=sum(!atlas$Hallmark_IFN_union),membership_changed=FALSE)
write_json(summary,file.path(OUT,'functional_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
saveRDS(list(config=cfg,universe=u,members=q,atlas=atlas,edges=all_edges,dictionary=all_terms,sets=sets,backgrounds=backgrounds,results=r,display=display),file.path(OUT,'data/processed/functional_atlas.rds'))
capture.output(sessionInfo(),file=file.path(OUT,'logs/functional_R_sessionInfo.txt'))
print(cv);print(summary[c('n_primary_tests','IFN_union_members','outside_IFN_union')]);print(display[,c('library','direction','term','k_overlap','FDR_global_BH')],row.names=FALSE)
