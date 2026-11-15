#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(edgeR);library(limma);library(Biobase);library(jsonlite);library(digest)})
args<-commandArgs(trailingOnly=TRUE);project<-args[1];root<-args[2]
cfg<-fromJSON(file.path(root,"parameters/experimental_consensus_v1.json"))
file<-file.path(project,"data/raw/GSE231788/GSE231788_SR92_genes_n_regRegions_raw_counts_for_degust_GEO.csv.gz")
raw<-read.csv(gzfile(file),check.names=FALSE,stringsAsFactors=FALSE)
pd<-pData(readRDS(file.path(project,"data/metadata/GSE231788_gse.rds"))[[1]])
meta<-read.csv(file.path(project,"data/processed/phase1/GSE231788_sample_metadata.csv"),stringsAsFactors=FALSE)
meta$condition<-factor(meta$condition,levels=c("mock","rsv"))
stopifnot(identical(as.integer(table(meta$condition)),c(3L,3L)))
meta$gsm<-vapply(seq_len(nrow(meta)),function(i){
  title<-paste0("A549 ",if(meta$condition[i]=="mock")"Mock" else "RSV","-infected, total RNA ",sub(".*_","",meta$sample_id[i]))
  hits<-pd$geo_accession[pd$title==title];stopifnot(length(hits)==1);hits
},character(1))
stopifnot(all(meta$sample_id%in%names(raw)))
ids<-sub("[.][0-9]+$","",as.character(raw$gene_id))
isgene<-grepl("^ENSG[0-9]+$",ids)
primary_assembly<-grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$",raw$chr)
canonical_ids<-ids[isgene&primary_assembly]
ambiguous<-names(which(table(canonical_ids)>1))
used<-isgene&primary_assembly&!ids%in%ambiguous
write.csv(data.frame(source_feature=raw$gene_symbol,gene_id=raw$gene_id,biotype=raw$biotype,chromosome=raw$chr,gene_feature=isgene,primary_assembly=primary_assembly,ambiguous_primary_id=ids%in%ambiguous,used_in_fit=used),file.path(root,"mapping/A549_feature_type_audit.csv"),row.names=FALSE)
gene<-raw[used,];gene$feature_id<-ids[used]
stopifnot(!anyDuplicated(gene$feature_id))
counts<-as.matrix(gene[,meta$sample_id,drop=FALSE]);storage.mode(counts)<-"double";rownames(counts)<-gene$feature_id
stopifnot(all(is.finite(counts)),all(counts>=0))
design<-model.matrix(~condition,meta)
y<-DGEList(counts=counts);keep<-filterByExpr(y,group=meta$condition)
y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y)
write.csv(data.frame(feature_id=rownames(counts),retained=keep),file.path(root,"qc/A549_feature_filter.csv"),row.names=FALSE)
write.csv(meta,file.path(root,"qc/A549_sample_manifest.csv"),row.names=FALSE)
write.csv(data.frame(feature_id=gene$feature_id,human_ensembl=gene$feature_id,display_symbol=gene$gene_symbol),file.path(root,"mapping/A549_feature_map.csv"),row.names=FALSE)
y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
r<-as.data.frame(topTags(glmQLFTest(fit,coef="conditionrsv"),n=Inf,sort.by="none")$table)
r<-data.frame(feature_id=rownames(r),r,check.names=FALSE)
r$primary_significant<-!is.na(r$FDR)&r$FDR<cfg$per_family_q_lt&abs(r$logFC)>cfg$per_family_abs_log2FC_gt
r<-r[order(r$PValue,r$feature_id),]
write.csv(r,file.path(root,"de/A549_primary__infected_vs_mock__edgeR_QLF.csv"),row.names=FALSE)
v<-voom(y,design,plot=FALSE);vf<-eBayes(lmFit(v,design),robust=TRUE)
t<-topTable(vf,coef="conditionrsv",number=Inf,sort.by="none")
lr<-data.frame(feature_id=rownames(t),logFC=t$logFC,PValue=t$P.Value,FDR=t$adj.P.Val,statistic=t$t)
lr$primary_significant<-!is.na(lr$FDR)&lr$FDR<cfg$per_family_q_lt&abs(lr$logFC)>cfg$per_family_abs_log2FC_gt
write.csv(lr,file.path(root,"de/A549_primary__infected_vs_mock__limma_voom.csv"),row.names=FALSE)
saveRDS(list(analysis_id="A549_primary",y=y,meta=meta,design=design,fit=fit),file.path(root,"models/A549_primary_QL_fit.rds"))
write_json(list(n_source_features=nrow(raw),n_Ensembl_genes=nrow(gene),n_ambiguous_primary_assembly_gene_ids_excluded=length(ambiguous),n_tested=nrow(y),n_fractional_count_entries=sum(abs(counts-round(counts))>1e-8),n_significant=sum(r$primary_significant),n_voom_significant=sum(lr$primary_significant),source_sha256=digest(file=file,algo="sha256"),status="Exploratory discovery; three culture profiles per group, not donors; alternate loci and residual ambiguous primary gene IDs excluded"),file.path(root,"de/A549_run_summary.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(root,"evidence/R42_A549_sessionInfo.txt"))
cat("A549 tested:",nrow(y),"significant:",sum(r$primary_significant),"\n")
