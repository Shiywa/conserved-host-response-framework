#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db);library(org.Mm.eg.db);library(jsonlite);library(digest)})
args<-commandArgs(trailingOnly=TRUE);project<-args[1];root<-args[2]
clean<-function(x)sub("[.][0-9]+$","",as.character(x))
one_value<-function(x){x<-sort(unique(x[!is.na(x)&nzchar(x)]));if(length(x)==1)x else if(length(x)>1)paste(x,collapse=";") else ""}
source_hashes<-list()
record_input<-function(p)source_hashes[[length(source_hashes)+1]]<<-data.frame(path=p,bytes=file.info(p)$size,sha256=digest(file=p,algo="sha256"))
orth_file<-file.path(root,"mapping/babelgene22.9_mouse_human_one_to_one.csv");record_input(orth_file)
orth<-read.csv(orth_file,stringsAsFactors=FALSE)
hfile<-file.path(project,"data/raw/GSE155237/GSE155237_gene_counts.csv.gz");record_input(hfile)
human<-read.csv(gzfile(hfile),check.names=FALSE,stringsAsFactors=FALSE)[[1]]
human<-clean(human)
hkeys<-intersect(human,keys(org.Hs.eg.db,keytype="ENSEMBL"))
ha<-select(org.Hs.eg.db,keys=hkeys,keytype="ENSEMBL",columns=c("SYMBOL","ENTREZID"))
write.csv(ha,file.path(root,"mapping/human_Ensembl_annotation_all_links.csv"),row.names=FALSE,na="")
ha_by<-split(ha,ha$ENSEMBL)
hmap<-data.frame(feature_id=human,human_ensembl=human,display_symbol=vapply(human,function(id)if(id%in%names(ha_by))one_value(ha_by[[id]]$SYMBOL) else "",character(1)),entrez_id=vapply(human,function(id)if(id%in%names(ha_by))one_value(ha_by[[id]]$ENTREZID) else "",character(1)),mapping_basis="Direct human Ensembl feature; annotation never gates membership",stringsAsFactors=FALSE)
write.csv(hmap,file.path(root,"mapping/human_feature_map.csv"),row.names=FALSE,na="")
mfile<-file.path(project,"data/processed/next_batch/GSE197868_mouse_lung_counts/GSE197868_mouse_lung_counts_gene_annotation.csv");record_input(mfile)
ma<-read.csv(mfile,stringsAsFactors=FALSE)
mm<-merge(data.frame(feature_id=ma$feature_id,ensembl=clean(ma$feature_id)),orth,by="ensembl",all.x=TRUE,sort=FALSE)
mm$eligible_mapping<-!is.na(mm$human_ensembl)&nzchar(mm$human_ensembl)
mm$mapping_basis<-"Mouse Ensembl -> frozen babelgene high-support one-to-one human Ensembl"
write.csv(mm,file.path(root,"mapping/mouse_feature_map.csv"),row.names=FALSE,na="")
rfile<-file.path(project,"data/raw/GSE297903/GSE297903_Niewiesk_36CottonRat_rnaseq_featurecounts_anno.txt.gz");record_input(rfile)
rat<-read.delim(gzfile(rfile),check.names=FALSE,stringsAsFactors=FALSE)[,c("GeneID","ENSID","GeneName")]
rat$mouse_protein<-clean(rat$ENSID)
valid<-intersect(unique(rat$mouse_protein[grepl("^ENSMUSP[0-9]+$",rat$mouse_protein)]),keys(org.Mm.eg.db,keytype="ENSEMBLPROT"))
pa<-select(org.Mm.eg.db,keys=valid,keytype="ENSEMBLPROT",columns=c("ENSEMBL","SYMBOL","ENTREZID"))
write.csv(pa,file.path(root,"mapping/cotton_rat_source_protein_to_mouse_all_links.csv"),row.names=FALSE,na="")
pg<-unique(pa[,c("ENSEMBLPROT","ENSEMBL")]);pg<-pg[!is.na(pg$ENSEMBL)&nzchar(pg$ENSEMBL),]
pd<-table(pg$ENSEMBLPROT);pg$protein_gene_degree<-as.integer(pd[pg$ENSEMBLPROT])
pg<-pg[pg$protein_gene_degree==1,]
rm<-merge(rat,pg,by.x="mouse_protein",by.y="ENSEMBLPROT",all.x=TRUE,sort=FALSE)
rm<-merge(rm,orth,by.x="ENSEMBL",by.y="ensembl",all.x=TRUE,sort=FALSE)
rm$provisional_mapping<-!is.na(rm$human_ensembl)&nzchar(rm$human_ensembl)
edges<-unique(rm[rm$provisional_mapping,c("GeneID","human_ensembl")])
source_degree<-table(edges$GeneID);target_degree<-table(edges$human_ensembl)
rm$source_feature_degree<-as.integer(source_degree[rm$GeneID]);rm$target_feature_degree<-as.integer(target_degree[rm$human_ensembl])
rm$eligible_mapping<-rm$provisional_mapping & !is.na(rm$source_feature_degree)&rm$source_feature_degree==1 & !is.na(rm$target_feature_degree)&rm$target_feature_degree==1
rm$feature_id<-rm$GeneID
rm$mapping_basis<-"Author-annotated cotton-rat mouse-protein bridge -> unambiguous mouse gene -> frozen mouse-human map; projection, not direct cotton-rat orthology"
rm$exclusion_reason<-ifelse(!rm$provisional_mapping,"No unambiguous source mouse-protein/orthology bridge",ifelse(!rm$eligible_mapping,"Nonunique cotton-rat feature <-> human gene relation",""))
write.csv(rm,file.path(root,"mapping/cotton_rat_feature_map.csv"),row.names=FALSE,na="")
org_file<-file.path(project,"data/processed/next_batch/GSE263272_organoid_tpm/GSE263272_organoid_tpm_gene_annotation.csv");record_input(org_file)
oa<-read.csv(org_file,stringsAsFactors=FALSE)
omap<-data.frame(feature_id=clean(oa$feature_id),human_ensembl=clean(oa$feature_id),display_symbol=oa$gene_symbol,mapping_basis="Direct human Ensembl in source TPM matrix",stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(omap$feature_id));write.csv(omap,file.path(root,"mapping/organoid_feature_map.csv"),row.names=FALSE,na="")
summary<-data.frame(system=c("human","mouse","cotton_rat","organoid"),n_source_features=c(nrow(hmap),nrow(mm),nrow(rat),nrow(omap)),n_eligible_mappings=c(nrow(hmap),sum(mm$eligible_mapping),sum(rm$eligible_mapping),nrow(omap)),basis=c("Direct human IDs, symbols for display only","Frozen high-support one-to-one mouse/human relation","Author annotation bridge; projection-qualified evidence","Direct human IDs; replicate eligibility separate"))
write.csv(summary,file.path(root,"mapping/mapping_summary.csv"),row.names=FALSE)
write.csv(do.call(rbind,source_hashes),file.path(root,"evidence/mapping_input_hashes.csv"),row.names=FALSE)
write_json(list(org_Hs_eg_db=as.character(packageVersion("org.Hs.eg.db")),org_Mm_eg_db=as.character(packageVersion("org.Mm.eg.db")),mouse_to_human="babelgene22.9 min_support=3 top=FALSE; one-to-one degree filtering",cotton_rat_caveat="Source mouse-protein annotation is an author-supplied bridge. One-to-one within this projection does not establish genome-wide evolutionary orthology."),file.path(root,"mapping/mapping_provenance.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(root,"evidence/R42_mapping_sessionInfo.txt"))
print(summary,row.names=FALSE)
