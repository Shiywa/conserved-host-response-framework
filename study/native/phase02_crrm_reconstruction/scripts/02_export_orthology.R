#!/usr/bin/env Rscript
# R 4.4 / babelgene 22.9: export a frozen full mouse-human evidence graph.
suppressPackageStartupMessages(library(babelgene))
args<-commandArgs(trailingOnly=TRUE);root<-args[1]
dir.create(file.path(root,"mapping"),recursive=TRUE,showWarnings=FALSE)
raw<-getFromNamespace("orthologs_df","babelgene")
mouse_ids<-unique(as.character(raw$ensembl[raw$taxon_id==10090 & !is.na(raw$ensembl)]))
mouse_ids<-mouse_ids[nzchar(mouse_ids)]
all<-orthologs(genes=mouse_ids,species="mouse",human=FALSE,min_support=3,top=FALSE)
stopifnot(all$taxon_id==10090,all$support_n>=3)
write.csv(all,file.path(root,"mapping/babelgene22.9_mouse_human_all_min3.csv"),row.names=FALSE)
edges<-unique(all[,c("ensembl","human_ensembl")])
edges<-edges[!is.na(edges$ensembl)&nzchar(edges$ensembl)&!is.na(edges$human_ensembl)&nzchar(edges$human_ensembl),]
mouse_degree<-table(edges$ensembl);human_degree<-table(edges$human_ensembl)
edges$mouse_degree<-as.integer(mouse_degree[edges$ensembl])
edges$human_degree<-as.integer(human_degree[edges$human_ensembl])
edges$eligible_one_to_one<-edges$mouse_degree==1 & edges$human_degree==1
edges<-merge(edges,all,by=c("ensembl","human_ensembl"),all.x=TRUE,sort=FALSE)
edges<-unique(edges)
write.csv(edges,file.path(root,"mapping/babelgene22.9_mouse_human_graph_audit.csv"),row.names=FALSE)
keep<-edges[edges$eligible_one_to_one,]
# One gene pair may have multiple historical symbol/Entrez metadata rows.
# Collapse metadata for the same stable-ID edge; never select by DE evidence.
collapse<-function(x)paste(sort(unique(x[!is.na(x)&nzchar(as.character(x))])),collapse=";")
groups<-split(seq_len(nrow(keep)),paste(keep$ensembl,keep$human_ensembl,sep="|"))
keep<-do.call(rbind,lapply(groups,function(i){
  z<-keep[i,,drop=FALSE]
  data.frame(ensembl=z$ensembl[1],human_ensembl=z$human_ensembl[1],
    mouse_degree=1L,human_degree=1L,eligible_one_to_one=TRUE,
    human_symbol=collapse(z$human_symbol),human_entrez=collapse(z$human_entrez),
    symbol=collapse(z$symbol),entrez=collapse(z$entrez),support=collapse(z$support),
    support_n=min(z$support_n),support_n_max=max(z$support_n),n_metadata_rows=nrow(z),stringsAsFactors=FALSE)
}))
stopifnot(!anyDuplicated(keep$ensembl),!anyDuplicated(keep$human_ensembl))
write.csv(keep,file.path(root,"mapping/babelgene22.9_mouse_human_one_to_one.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,"evidence/R44_orthology_sessionInfo.txt"))
capture.output(packageDescription("babelgene"),file=file.path(root,"evidence/babelgene_package_description.txt"))
cat("All supported mouse-human pairs:",nrow(all),"; eligible one-to-one:",nrow(keep),"\n")
