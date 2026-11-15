suppressPackageStartupMessages({library(Biobase);library(AnnotationDbi);library(org.Hs.eg.db);library(jsonlite);library(digest);library(data.table)})
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)]);OUT<-file.path(Sys.getenv('CRRM_WORKSPACE'),'phase22_Result3_incremental_information_20261112')
REV<-file.path(Sys.getenv('CRRM_WORKSPACE'),"");OLD<-file.path(REV,'phase20_Result3_clinical_analysis_20261112');P<-dirname(dirname(REV))
cfg<-fromJSON(file.path(OUT,'parameters/exploration_plan_v1.json'));receipt<-fromJSON(file.path(OUT,'parameters/declaration_receipt.json'));stopifnot(digest(file=file.path(OUT,'parameters/exploration_plan_v1.json'),algo='sha256')==receipt$plan_sha256)
for(i in seq_len(nrow(cfg$inputs)))stopifnot(digest(file=cfg$inputs$path[i],algo='sha256')==cfg$inputs$sha256[i])
write_tab<-function(x,f)write.csv(x,file.path(OUT,'tables',f),row.names=FALSE,na='')
core<-read.csv(file.path(OUT,'parameters/frozen_CRRM388.csv'));branch<-read.csv(file.path(OUT,'parameters/frozen_functional_branches.csv'));branch<-branch[branch$branch=='Neutrophil degranulation',]
am<-fromJSON(file.path(OUT,'parameters/primary_assembly_panel_annotation_amendment.json'));stopifnot(am$declared_before_any_new_scores_or_tests,am$parent_plan_sha256==receipt$plan_sha256);overrides<-vapply(am$overrides,function(x)x$human_ensembl,character(1));stopifnot(!any(overrides%in%core$human_ensembl))
for(n in names(am$overrides))stopifnot(digest(file=file.path(OUT,'data/public',paste0(n,'_Ensembl_primary_lookup.json')),algo='sha256')==am$overrides[[n]]$source_sha256)
official<-select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,keytype='SYMBOL'),keytype='SYMBOL',columns='ENSEMBL');official<-unique(official[!is.na(official$ENSEMBL),c('SYMBOL','ENSEMBL')]);env<-list2env(split(official$ENSEMBL,official$SYMBOL),hash=TRUE,parent=emptyenv())
gmt<-strsplit(readLines(file.path(OUT,'data/public/hallmark_2026_1.gmt')),'\t',fixed=TRUE);ifn<-unique(unlist(lapply(gmt,function(x)if(x[1]%in%c('HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE'))x[-c(1,2)]else NULL)))
pan<-read.csv(file.path(OUT,'parameters/published_panel_members.csv'));symbols<-unique(c(pan$symbol,ifn));unknown<-setdiff(symbols,ls(env));valid<-intersect(unknown,keys(org.Hs.eg.db,keytype='ALIAS'));al<-if(length(valid))select(org.Hs.eg.db,keys=valid,keytype='ALIAS',columns='ENSEMBL')else data.frame(ALIAS=character(),ENSEMBL=character());al<-al[!is.na(al$ENSEMBL),,drop=FALSE];ae<-list2env(split(al$ENSEMBL,al$ALIAS),hash=TRUE,parent=emptyenv())
resolve<-function(s){if(s%in%names(overrides))return(unname(overrides[s]));v<-get0(s,envir=env,inherits=FALSE);if(is.null(v))v<-get0(s,envir=ae,inherits=FALSE);v<-unique(v[!is.na(v)]);if(length(v)==1)v else ''}
sm<-data.frame(source_symbol=symbols,human_ensembl=vapply(symbols,resolve,character(1)),stringsAsFactors=FALSE);write_tab(sm,'00_comparator_symbol_resolution.csv')
gene_sets<-list(CRRM25=branch$human_ensembl,CRRM24=branch$human_ensembl[branch$display_symbol!='LCN2'],CRRM388=core$human_ensembl)
for(n in c('Zhang10','ZhangSeverityPair'))gene_sets[[n]]<-unique(sm$human_ensembl[match(pan$symbol[pan$set_name==n],sm$source_symbol)])
gene_sets$IFNUnion<-unique(sm$human_ensembl[match(ifn,sm$source_symbol)]);gene_sets<-lapply(gene_sets,function(v)v[nzchar(v)])
gene_sets$CRRMDisjoint<-setdiff(gene_sets$CRRM24,gene_sets$IFNUnion)
markers<-read.delim(file.path(OUT,'data/public/MCPcounter_genes.txt'),check.names=FALSE,stringsAsFactors=FALSE);excluded<-unique(unlist(gene_sets[c('CRRM388','Zhang10','IFNUnion')]))
for(z in list(c('Neutrophils','MCPNeutrophils'),c('Monocytic lineage','MCPMonocytic'))){
    g<-unique(markers[['ENSEMBL ID']][markers[['Cell population']]==z[1]]);g<-g[!is.na(g)&nzchar(g)];gene_sets[[z[2]]]<-setdiff(g,excluded)
}
nominal<-c(CRRM25=25,CRRM24=24,CRRM388=388,Zhang10=10,ZhangSeverityPair=2,IFNUnion=length(ifn),CRRMDisjoint=length(gene_sets$CRRMDisjoint),MCPNeutrophils=length(gene_sets$MCPNeutrophils),MCPMonocytic=length(gene_sets$MCPMonocytic))
stopifnot(length(gene_sets$CRRM25)==25,length(gene_sets$CRRM24)==24,length(gene_sets$CRRMDisjoint)>=5,identical(intersect(gene_sets$CRRM25,gene_sets$Zhang10),'ENSG00000148346'))
members<-do.call(rbind,lapply(names(gene_sets),function(n)data.frame(set_name=n,human_ensembl=gene_sets[[n]],stringsAsFactors=FALSE)));members$display_symbol<-mapIds(org.Hs.eg.db,keys=members$human_ensembl,keytype='ENSEMBL',column='SYMBOL',multiVals='first');members$display_symbol[is.na(members$display_symbol)]<-members$human_ensembl[is.na(members$display_symbol)]
write.csv(members,file.path(OUT,'parameters/resolved_score_members.csv'),row.names=FALSE);write_json(as.list(nominal),file.path(OUT,'parameters/nominal_set_sizes.json'),pretty=TRUE,auto_unbox=TRUE)
over<-do.call(rbind,lapply(names(gene_sets),function(n)do.call(rbind,lapply(names(gene_sets),function(k)data.frame(set1=n,set2=k,n_overlap=length(intersect(gene_sets[[n]],gene_sets[[k]])))))));write_tab(over,'01_signature_overlaps.csv')
target<-unique(unlist(gene_sets));message('Nominal sets:');print(nominal);message('Unique target genes:',length(target));all_cov<-list();probe_override_audit<-list()
series<-function(file){con<-gzfile(file);l<-readLines(con,warn=FALSE);close(con);b<-grep('^!series_matrix_table_begin',l);e<-grep('^!series_matrix_table_end',l);stopifnot(length(b)==1,length(e)==1);d<-read.delim(textConnection(paste(l[(b+1):(e-1)],collapse='\n')),check.names=FALSE,quote='"');x<-as.matrix(d[,-1]);storage.mode(x)<-'double';rownames(x)<-d[[1]];x}
maybe_log<-function(x)if(unname(quantile(x[is.finite(x)],.99))>100)log2(pmax(x,0)+1)else x
for(ds in c('GSE246622','GSE246622_RMA','GSE69606','GSE105450')){
    m<-read.csv(file.path(OLD,'data/processed',paste0(ds,'_metadata.csv')),check.names=FALSE,stringsAsFactors=FALSE)
    if(ds=='GSE246622'){
        es<-readRDS(file.path(P,'data/metadata/GSE246622_gse.rds'))[[1]];xp<-exprs(es)[,m$geo_accession,drop=FALSE];colnames(xp)<-m$sample_id
    }else if(ds=='GSE246622_RMA'){
        d<-as.data.frame(fread(file.path(P,'data/processed/phase1_validation/GSE246622_expression_matrix.csv'),check.names=FALSE));xp<-as.matrix(d[,-1]);storage.mode(xp)<-'double';rownames(xp)<-d[[1]];xp<-xp[,m$sample_id,drop=FALSE]
    }else xp<-maybe_log(series(file.path(P,'data/raw',ds,paste0(ds,'_series_matrix.txt.gz')))[,m$sample_id,drop=FALSE])
    ann<-read.csv(gzfile(file.path(OLD,'mapping',paste0(ds,'_full_reference_mapping.csv.gz'))),stringsAsFactors=FALSE)
    for(sym in names(overrides)){
        ix<-which(ann$source_symbol==sym&!ann$annotation_conflict&ann$n_author_ensembl==0&ann$complete_expression)
        if(length(ix)){stopifnot(all(!nzchar(ann$human_ensembl[ix])|ann$human_ensembl[ix]==overrides[[sym]]));probe_override_audit[[length(probe_override_audit)+1]]<-data.frame(dataset=ds,feature_id=ann$feature_id[ix],source_symbol=sym,old_human_ensembl=ann$human_ensembl[ix],verified_primary_ensembl=overrides[[sym]]);ann$human_ensembl[ix]<-overrides[[sym]];ann$unambiguous[ix]<-TRUE}
    }
    ann<-unique(ann[ann$unambiguous&ann$complete_expression&ann$human_ensembl%in%target,c('feature_id','human_ensembl')]);stopifnot(all(ann$feature_id%in%rownames(xp)));ann<-ann[!duplicated(ann),]
    groups<-split(ann$feature_id,ann$human_ensembl);x<-do.call(rbind,lapply(groups,function(ids){v<-xp[unique(ids),,drop=FALSE];if(nrow(v)==1)as.numeric(v[1,])else apply(v,2,median)}));rownames(x)<-names(groups);colnames(x)<-m$sample_id
    stopifnot(all(is.finite(x)),identical(colnames(x),m$sample_id));ref<-m$reference_for_scaling;sd0<-apply(x[,ref,drop=FALSE],1,sd);scored<-names(sd0)[is.finite(sd0)&sd0>1e-8]
    for(n in names(gene_sets))all_cov[[length(all_cov)+1]]<-data.frame(dataset=ds,set_name=n,n_nominal=nominal[[n]],n_resolved=length(gene_sets[[n]]),n_present=sum(gene_sets[[n]]%in%rownames(x)),n_measured=sum(gene_sets[[n]]%in%scored),measured_fraction=sum(gene_sets[[n]]%in%scored)/nominal[[n]])
    write.csv(data.frame(human_ensembl=rownames(x),x,check.names=FALSE),gzfile(file.path(OUT,'data/processed',paste0(ds,'_target_log_expression.csv.gz'))),row.names=FALSE);write.csv(m,file.path(OUT,'data/processed',paste0(ds,'_metadata.csv')),row.names=FALSE,na='');message(ds,': ',nrow(x),'genes; ',ncol(x),'profiles')
    rm(xp,x);invisible(gc())
}
cv<-do.call(rbind,all_cov);write_tab(cv,'02_comparator_coverage.csv');write_tab(do.call(rbind,probe_override_audit),'00b_primary_assembly_override_probes.csv')
counts<-read.table(file.path(OUT,'data/public/RESCEU_datTraits_counts.txt'),header=TRUE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
traits<-read.table(file.path(OUT,'data/public/RESCEU_datTraits.txt'),header=TRUE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(counts)==66,ncol(counts)==257,!anyDuplicated(rownames(counts)),all(as.character(1:256)%in%names(counts)))
write_tab(data.frame(sample_id=rownames(counts),counts,check.names=FALSE),'03_author_CyTOF_cluster_frequencies.csv');write_tab(traits,'04_author_clinical_traits.csv')
dump<-list(prepared=TRUE,nominal_sizes=as.list(nominal),target_genes=length(target),composition='Author processed measured256-cluster frequencies; no rawFCS analysis');write_json(dump,file.path(OUT,'logs/input_preparation_summary.json'),pretty=TRUE,auto_unbox=TRUE);capture.output(sessionInfo(),file=file.path(OUT,'logs/input_R_sessionInfo.txt'));print(cv)
