#!/usr/bin/env Rscript
# Actual corrected bulk expression modeling. Does not read validation/scRNA data.
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite);library(digest)})
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
args<-commandArgs(trailingOnly=TRUE)
project<-args[1]; phase1<-args[2]; root<-args[3]
for(d in c("de","qc","models","logs","evidence")) dir.create(file.path(root,d),recursive=TRUE,showWarnings=FALSE)
config_file<-file.path(root,"parameters/pre_analysis_v1.json")
cfg<-fromJSON(config_file)
config_sha<-digest(file=config_file,algo="sha256")
timestamp<-function()format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z",tz="Asia/Shanghai")
write_status<-function(stage,status,detail="")write_json(list(stage=stage,status=status,detail=detail,time=timestamp(),config_sha256=config_sha),file.path(root,"logs/bulk_status.json"),pretty=TRUE,auto_unbox=TRUE)
read_counts<-function(file){
  t<-read.csv(if(grepl("[.]gz$",file))gzfile(file) else file,check.names=FALSE,stringsAsFactors=FALSE)
  ids<-sub("[.][0-9]+$","",as.character(t[[1]])); stopifnot(!anyDuplicated(ids))
  m<-as.matrix(t[,-1,drop=FALSE]);storage.mode(m)<-"double";rownames(m)<-ids
  stopifnot(all(is.finite(m)),all(m>=0),all(colSums(m)>0));m
}
read_meta<-function(file)read.csv(file,check.names=FALSE,stringsAsFactors=FALSE)
flag<-function(x)tolower(as.character(x))%in%c("true","1")
metrics<-list();inputs<-list()
add_input<-function(p){inputs[[length(inputs)+1]]<<-data.frame(path=p,bytes=file.info(p)$size,sha256=digest(file=p,algo="sha256"))}
write_table<-function(tab,file)write.csv(tab,file,row.names=FALSE,na="",fileEncoding="UTF-8")
contrast_vector<-function(weights,design){
  stopifnot(all(names(weights)%in%colnames(design)))
  c<-setNames(rep(0,ncol(design)),colnames(design));c[names(weights)]<-weights;c
}
run_count<-function(id,counts,meta,design,groups,contrasts,role){
  message(timestamp()," START ",id," n=",ncol(counts)," p=",ncol(design))
  write_status(id,"running")
  stopifnot(nrow(meta)==ncol(counts),nrow(design)==ncol(counts),qr(design)$rank==ncol(design))
  y<-DGEList(counts=counts)
  keep<-filterByExpr(y,group=groups)
  y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method="TMM")
  stopifnot(nrow(y)>10)
  qct<-data.frame(count_matrix_column=colnames(y),library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors,detected_features=colSums(counts>0),meta,check.names=FALSE)
  write_table(qct,file.path(root,"qc",paste0(id,"_samples.csv")))
  write_table(data.frame(feature_id=rownames(counts),retained=keep),file.path(root,"qc",paste0(id,"_feature_filter.csv")))
  write.csv(design,file.path(root,"qc",paste0(id,"_design.csv")),row.names=TRUE)
  y<-estimateDisp(y,design,robust=TRUE)
  fit<-glmQLFit(y,design,robust=TRUE)
  v<-voom(y,design,plot=FALSE)
  vf<-lmFit(v,design)
  for(label in names(contrasts)){
    c<-contrast_vector(contrasts[[label]],design)
    tested<-glmQLFTest(fit,contrast=c)
    r<-as.data.frame(topTags(tested,n=Inf,sort.by="none")$table)
    r<-data.frame(feature_id=rownames(r),r,check.names=FALSE)
    r$primary_significant<-!is.na(r$FDR) & r$FDR<cfg$primary_q_threshold & abs(r$logFC)>cfg$absolute_log2FC_threshold
    r$contrast<-label;r$analysis_id<-id;r$method<-"edgeR_QLF";r$role<-role
    r<-r[order(r$PValue,r$feature_id),]
    write_table(r,file.path(root,"de",paste0(id,"__",label,"__edgeR_QLF.csv")))
    lf<-contrasts.fit(vf,matrix(c,ncol=1,dimnames=list(names(c),label)))
    lf<-eBayes(lf,robust=TRUE)
    lr<-topTable(lf,coef=1,number=Inf,sort.by="none")
    lr<-data.frame(feature_id=rownames(lr),logFC=lr$logFC,logCPM=lr$AveExpr,statistic=lr$t,PValue=lr$P.Value,FDR=lr$adj.P.Val,check.names=FALSE)
    lr$primary_significant<-!is.na(lr$FDR) & lr$FDR<cfg$primary_q_threshold & abs(lr$logFC)>cfg$absolute_log2FC_threshold
    lr$contrast<-label;lr$analysis_id<-id;lr$method<-"limma_voom";lr$role<-"Method sensitivity; not used to select minimum FDR"
    lr<-lr[order(lr$PValue,lr$feature_id),]
    write_table(lr,file.path(root,"de",paste0(id,"__",label,"__limma_voom.csv")))
    metrics[[length(metrics)+1]]<<-data.frame(analysis_id=id,contrast=label,role=role,n_samples=ncol(y),n_features_input=nrow(counts),n_features_tested=nrow(y),design_rank=qr(design)$rank,residual_df=nrow(design)-ncol(design),primary_q_lt_005=sum(r$FDR<0.05),primary_q_lfc=sum(r$primary_significant),primary_up=sum(r$primary_significant&r$logFC>0),primary_down=sum(r$primary_significant&r$logFC<0),voom_q_lfc=sum(lr$primary_significant))
    message(timestamp()," RESULT ",id,"/",label," tested=",nrow(y)," significant=",sum(r$primary_significant))
  }
  saveRDS(list(analysis_id=id,y=y,meta=meta,design=design,filter_groups=groups,contrasts=contrasts,fit=fit,config_sha256=config_sha),file.path(root,"models",paste0(id,"_QL_fit.rds")),compress="gzip")
  write_table(do.call(rbind,metrics),file.path(root,"de/DE_run_summary.csv"))
  write_status(id,"complete")
  invisible(NULL)
}
human_meta<-read_meta(file.path(phase1,"tables/06_human_challenge_104_record_manifest.csv"))
add_input(file.path(phase1,"tables/06_human_challenge_104_record_manifest.csv"))
human_p<-file.path(project,"data/raw/GSE155237/GSE155237_gene_counts.csv.gz")
add_input(human_p);human_counts<-read_counts(human_p)
prep_human<-function(meta,day){
  meta$donor<-factor(meta$donor_id_raw)
  meta$outcome<-factor(meta$eventual_infection_outcome,levels=c("Uninfected","Infected"))
  meta$post<-as.numeric(meta$timepoint==day)
  meta$post_infected<-meta$post*as.numeric(meta$outcome=="Infected")
  meta$batch<-factor(meta$batch);meta
}
extract<-function(counts,meta)counts[,meta$expression_column,drop=FALSE]
h<-prep_human(human_meta[flag(human_meta$include_primary_D3),],"D3")
stopifnot(nrow(h)==76,length(unique(h$donor))==38)
ct<-list(infected_change=c(post=1,post_infected=1),outcome_by_time=c(post_infected=1),uninfected_change=c(post=1))
run_count("H_D3_primary",extract(human_counts,h),h,model.matrix(~donor+post+post_infected+batch,h),interaction(h$outcome,h$post),ct,"Primary human family; only infected_change supplies anchor")
run_count("H_D3_no_batch",extract(human_counts,h),h,model.matrix(~donor+post+post_infected,h),interaction(h$outcome,h$post),ct,"Batch sensitivity, same family")
pairs<-read_meta(file.path(phase1,"tables/09_human_longitudinal_pairs.csv"))
stable<-pairs$donor_id[pairs$contrast_id=="human_D3_vs_D0" & !flag(pairs$batch_changes)]
hn<-droplevels(h[h$donor_id_raw%in%stable,]);stopifnot(nrow(hn)==66)
run_count("H_D3_new_new",extract(human_counts,hn),hn,model.matrix(~donor+post+post_infected,hn),interaction(hn$outcome,hn$post),ct,"33 donor-pair stable-batch sensitivity")
d7_p<-file.path(project,"data/raw/GSE166161/GSE166161_gene_counts.csv.gz");add_input(d7_p)
d7_counts<-read_counts(d7_p)
d7<-prep_human(human_meta[human_meta$accession=="GSE166161",],"D7")
stopifnot(nrow(d7)==26,length(unique(d7$donor))==13)
run_count("H_D7_same_family",extract(d7_counts,d7),d7,model.matrix(~donor+post+post_infected,d7),interaction(d7$outcome,d7$post),ct,"Temporal sensitivity, not independent replication")
baseline<-prep_human(human_meta[human_meta$accession=="GSE155237" & human_meta$timepoint=="D0",],"D3")
run_count("H_baseline",extract(human_counts,baseline),baseline,model.matrix(~batch+outcome,baseline),baseline$outcome,list(outcome_association=c(outcomeInfected=1)),"Baseline susceptibility/context only; excluded from construction")
pooled<-human_meta[!flag(human_meta$duplicate_alias)&human_meta$donor_id_raw!="JF_RF017",]
pooled<-prep_human(pooled,"D3")
pooled$day3<-as.numeric(pooled$timepoint=="D3");pooled$day7<-as.numeric(pooled$timepoint=="D7")
pooled$day3_infected<-pooled$day3*as.numeric(pooled$outcome=="Infected")
pooled$day7_infected<-pooled$day7*as.numeric(pooled$outcome=="Infected")
stopifnot(identical(rownames(human_counts),rownames(d7_counts)))
merged<-cbind(human_counts,d7_counts[,!colnames(d7_counts)%in%colnames(human_counts),drop=FALSE])
stopifnot(!anyDuplicated(colnames(merged)),nrow(pooled)==90)
run_count("H_joint_time",extract(merged,pooled),pooled,model.matrix(~donor+day3+day7+day3_infected+day7_infected+batch,pooled),interaction(pooled$outcome,pooled$timepoint),list(D3_infected_change=c(day3=1,day3_infected=1),D7_infected_change=c(day7=1,day7_infected=1)),"Same-family joint temporal sensitivity; baseline aliases used once")
mouse_block<-"GSE197868_mouse_lung_counts"
mouse_p<-file.path(project,"data/processed/next_batch",mouse_block,paste0(mouse_block,"_expression_matrix.csv"))
mouse_m<-file.path(project,"data/processed/next_batch",mouse_block,paste0(mouse_block,"_sample_metadata.csv"))
add_input(mouse_p);add_input(mouse_m)
mouse<-read_meta(mouse_m);mouse$condition<-factor(mouse$condition,levels=c("Control","RSV"))
mouse_counts<-read_counts(mouse_p)[,mouse$sample_id,drop=FALSE]
run_count("M_lung_primary",mouse_counts,mouse,model.matrix(~condition,mouse),mouse$condition,list(infected_vs_control=c(conditionRSV=1)),"Primary mouse family, lung only")
rat_p<-file.path(project,"data/raw/GSE297903/GSE297903_Niewiesk_36CottonRat_rnaseq_featurecounts_anno.txt.gz")
rat_m<-file.path(project,"data/processed/phase1/GSE297903_sample_metadata.csv")
add_input(rat_p);add_input(rat_m)
rat_raw<-read.delim(gzfile(rat_p),check.names=FALSE,stringsAsFactors=FALSE)
rat_meta<-read_meta(rat_m)
rat_counts<-as.matrix(rat_raw[,rat_meta$sample_id,drop=FALSE]);storage.mode(rat_counts)<-"double";rownames(rat_counts)<-rat_raw$GeneID
stopifnot(!anyDuplicated(rownames(rat_counts)))
for(day in c("day1","day4")){
  r<-rat_meta[rat_meta$time_point%in%c("day0",day),]
  r$group<-factor(paste(r$age_group,ifelse(r$time_point=="day0","PBS","RSV"),sep="_"),levels=c("adult_PBS","adult_RSV","geriatric_PBS","geriatric_RSV"))
  contrasts<-list(age_average=c(groupadult_PBS=-0.5,groupadult_RSV=0.5,groupgeriatric_PBS=-0.5,groupgeriatric_RSV=0.5),adult=c(groupadult_PBS=-1,groupadult_RSV=1),geriatric=c(groupgeriatric_PBS=-1,groupgeriatric_RSV=1))
  run_count(paste0("R_",ifelse(day=="day1","D1_primary","D4_secondary")),rat_counts[,r$sample_id,drop=FALSE],r,model.matrix(~0+group,r),r$group,contrasts,ifelse(day=="day1","Primary cotton-rat family: age_average only; author-annotated projection","Secondary day4, no parallel day4 mock"))
}
org_block<-"GSE263272_organoid_tpm"
org_p<-file.path(project,"data/processed/next_batch",org_block,paste0(org_block,"_expression_matrix.csv"))
org_m<-file.path(project,"data/processed/next_batch",org_block,paste0(org_block,"_sample_metadata.csv"))
add_input(org_p);add_input(org_m)
org<-read_meta(org_m);org$condition<-factor(org$condition,levels=c("Mock","RSV_A","RSV_B"))
tpm<-read_counts(org_p)[,org$sample_id,drop=FALSE]
keep<-rowSums(tpm>=1)>=4
org_design<-model.matrix(~0+condition,org)
ofit<-lmFit(log2(tpm[keep,,drop=FALSE]+1),org_design)
for(label in c("strain_average","RSV_A","RSV_B")){
  weights<-switch(label,strain_average=c(conditionMock=-1,conditionRSV_A=0.5,conditionRSV_B=0.5),RSV_A=c(conditionMock=-1,conditionRSV_A=1),RSV_B=c(conditionMock=-1,conditionRSV_B=1))
  c<-contrast_vector(weights,org_design)
  f<-eBayes(contrasts.fit(ofit,matrix(c,ncol=1)),trend=TRUE,robust=TRUE)
  t<-topTable(f,number=Inf,sort.by="none")
  r<-data.frame(feature_id=rownames(t),logFC=t$logFC,AveExpr=t$AveExpr,statistic=t$t,PValue=t$P.Value,FDR=t$adj.P.Val)
  r$primary_significant<-r$FDR<cfg$primary_q_threshold & abs(r$logFC)>cfg$absolute_log2FC_threshold
  r$contrast<-label;r$analysis_id<-"O_conditional";r$method<-"limma_logTPM";r$role<-"Conditional only; independent experimental units not yet verified; not an eligible conservation vote"
  r<-r[order(r$PValue,r$feature_id),]
  write_table(r,file.path(root,"de",paste0("O_conditional__",label,"__limma_logTPM.csv")))
  metrics[[length(metrics)+1]]<-data.frame(analysis_id="O_conditional",contrast=label,role=r$role[1],n_samples=ncol(tpm),n_features_input=nrow(tpm),n_features_tested=sum(keep),design_rank=3,residual_df=9,primary_q_lt_005=sum(r$FDR<0.05),primary_q_lfc=sum(r$primary_significant),primary_up=sum(r$primary_significant&r$logFC>0),primary_down=sum(r$primary_significant&r$logFC<0),voom_q_lfc=NA_integer_)
}
write_table(do.call(rbind,metrics),file.path(root,"de/DE_run_summary.csv"))
write_table(do.call(rbind,inputs),file.path(root,"evidence/bulk_input_hashes.csv"))
capture.output(sessionInfo(),file=file.path(root,"evidence/R42_fit_sessionInfo.txt"))
write_status("all_corrected_bulk_models","complete","Actual corrected DE completed; no new core has yet been externally validated")
message(timestamp()," ALL CORRECTED BULK MODELS COMPLETE")
