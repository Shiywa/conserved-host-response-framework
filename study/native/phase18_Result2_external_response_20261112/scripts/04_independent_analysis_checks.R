suppressPackageStartupMessages({library(jsonlite);library(digest)})
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)]);out<-file.path(Sys.getenv('CRRM_WORKSPACE'),'phase18_Result2_external_response_20261112')
core<-read.csv(file.path(out,'parameters/frozen_CRRM388.csv'));br<-read.csv(file.path(out,'parameters/frozen_functional_branches.csv'));objects<-readRDS(file.path(out,'data/processed/external_log_expression.rds'))
all_scores<-read.csv(file.path(out,'tables/05_sample_scores.csv'),check.names=FALSE);primary<-read.csv(file.path(out,'tables/07_primary_response_effects.csv'));effects<-read.csv(file.path(out,'tables/10_external_gene_effects.csv'));functions<-read.csv(file.path(out,'tables/08_primary_functional_effects.csv'));components<-read.csv(file.path(out,'tables/09_directional_component_effects.csv'))
checks<-list();ck<-function(name,ok){checks[[length(checks)+1]]<<-data.frame(check=name,passed=isTRUE(ok));stopifnot(isTRUE(ok))}
near<-function(a,b,tol=1e-10)isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol,check.attributes=FALSE))
receipt<-fromJSON(file.path(out,'parameters/declaration_receipt.json'));ck('Analysis plan hash unchanged',digest(file=file.path(out,'parameters/analysis_plan_v1.json'),algo='sha256')==receipt$plan_sha256)
ck('Frozen388unchanged',digest(file=file.path(out,'parameters/frozen_CRRM388.csv'),algo='sha256')==receipt$core_sha256);ck('Frozen7functionalsets unchanged',digest(file=file.path(out,'parameters/frozen_functional_branches.csv'),algo='sha256')==receipt$branches_sha256)
scores<-list();genes<-list()
for(name in names(objects)){
  o<-objects[[name]];ref<-o$reference;sd0<-apply(o$x[,ref,drop=FALSE],1,sd);keep<-is.finite(sd0)&sd0>1e-8;x<-o$x[keep,,drop=FALSE];z<-sweep(sweep(x,1,rowMeans(x[,ref,drop=FALSE]),'-'),1,sd0[keep],'/');sign<-ifelse(core$frozen_direction[match(rownames(z),core$human_ensembl)]=='up',1,-1);directed<-sweep(z,1,sign,'*')
  score<-colMeans(directed);want<-all_scores[all_scores$dataset==name,];score<-score[match(want$sample_id,names(score))]
  message(name,': maximum score difference = ',max(abs(score-want$score_Total)), '; missing matches = ',sum(is.na(score)))
  ck(paste(name,'independent frozen score reconstruction'),near(score,want$score_Total));ck(paste(name,'component sum equals total'),near(want$score_Total,want[['score_Up contribution']]+want[['score_Down contribution']]))
  for(branch in unique(br$branch)){
    ids<-intersect(br$human_ensembl[br$branch==branch],rownames(directed));eligible<-length(ids)>=5&&length(ids)/sum(br$branch==branch)>=.5
    v<-if(eligible)colMeans(directed[ids,,drop=FALSE])else rep(NA_real_,ncol(directed));ck(paste(name,branch,'membership and scores'),near(v,want[[paste0('score_',branch)]]))
  }
  scores[[name]]<-score;names(scores[[name]])<-want$sample_id;genes[[name]]<-x
}
compare_primary<-function(key,est,ci,p,genecoef){
  r<-primary[primary$contrast_key==key,];ck(paste(key,'effect'),near(est,r$effect));ck(paste(key,'95percentCI'),near(ci,c(r$CI95_low,r$CI95_high)));ck(paste(key,'modelP'),near(p,r$p_value))
  tab<-effects[effects$contrast_key==key,];ck(paste(key,'all gene coefficients'),near(genecoef[match(tab$human_ensembl,names(genecoef))],tab$log_expression_effect))
}
o<-objects$GSE247298;m<-o$meta;donors<-sort(unique(m$donor_id[m$background=='Non-CF']));i0<-vapply(donors,function(d)which(m$donor_id==d&m$time_h==0),integer(1));i1<-vapply(donors,function(d)which(m$donor_id==d&m$time_h==24),integer(1));v<-scores$GSE247298[i1]-scores$GSE247298[i0];tt<-t.test(v);ge<-rowMeans(genes$GSE247298[,i1]-genes$GSE247298[,i0]);compare_primary('Airway24h',mean(v),tt$conf.int,tt$p.value,ge)

o<-objects$GSE112956;m<-o$meta;i<-which(m$time_d==1);d<-data.frame(cell=factor(paste(m$strain[i],m$condition[i],sep='_')),score=scores$GSE112956[i]);fit<-lm(score~0+cell,d);c0<-ifelse(grepl('_RSV$',names(coef(fit))),.2,-.2);est<-sum(c0*coef(fit));se<-sqrt(drop(t(c0)%*%vcov(fit)%*%c0));ci<-est+c(-1,1)*qt(.975,df.residual(fit))*se;p<-2*pt(-abs(est/se),df.residual(fit));X<-model.matrix(fit);ge<-drop(t(c0)%*%solve(crossprod(X),crossprod(X,t(genes$GSE112956[,i]))));names(ge)<-rownames(genes$GSE112956);compare_primary('Mouse_multistrain_D1',est,ci,p,ge);ck('Actual multistrain groups14plus15primary',sum(m$condition=='Control')==14&&sum(m$condition=='RSV'&m$time_d==1)==15)

for(key in c('Mouse119849','Mouse140226_D5','Mouse18170','Hamster_lung')){
  dataset<-primary$dataset[primary$contrast_key==key];o<-objects[[dataset]];m<-o$meta;keep<-if(dataset=='GSE140226')m$time_d==5 else rep(TRUE,nrow(m));a<-which(keep&m$condition=='RSV');b<-which(keep&m$condition=='Control');tt<-t.test(scores[[dataset]][a],scores[[dataset]][b]);ge<-rowMeans(genes[[dataset]][,a,drop=FALSE])-rowMeans(genes[[dataset]][,b,drop=FALSE]);compare_primary(key,mean(scores[[dataset]][a])-mean(scores[[dataset]][b]),tt$conf.int,tt$p.value,ge)
}

# Verify the Python donor-change regression against the full76-profile donor-fixed model.
o<-objects$Human_challenge;m<-o$meta;keep<-tolower(as.character(m$primary))=='true';d<-m[keep,];d$donor_id<-factor(d$donor_id);d$post<-as.numeric(d$time_d==3);d$post_infected<-d$post*as.numeric(d$condition=='Infected');d$batch<-factor(d$batch);d$score<-scores$Human_challenge[keep]
fit<-lm(score~donor_id+post+post_infected+batch,d);cn<-names(coef(fit));c0<-as.numeric(cn%in%c('post','post_infected'));est<-sum(c0*coef(fit));se<-sqrt(drop(t(c0)%*%vcov(fit)%*%c0));ci<-est+c(-1,1)*qt(.975,df.residual(fit))*se;p<-2*pt(-abs(est/se),df.residual(fit));X<-model.matrix(fit);ge<-drop(t(c0)%*%solve(crossprod(X),crossprod(X,t(genes$Human_challenge[,keep]))));names(ge)<-rownames(genes$Human_challenge);compare_primary('Human_D3',est,ci,p,ge)
ck('Human primary38paired donors',nrow(d)==76&&nlevels(d$donor_id)==38);ck('Human full model residualdf35',df.residual(fit)==35)
ck('Holm7matchesR',near(p.adjust(primary$p_value,'holm',n=7),primary$primary_Holm_P));ck('BH49matchesR',near(p.adjust(functions$p_value,'BH',n=49),functions$FDR_BH49));ck('BH14matchesR',near(p.adjust(components$p_value,'BH',n=14),components$FDR_BH14))
for(key in primary$contrast_key){v<-components[components$contrast_key==key,];ck(paste(key,'effect components sum'),near(sum(v$effect),primary$effect[primary$contrast_key==key]))}
ck('No negative/weak primary result deleted',nrow(primary)==7&&all(tolower(as.character(primary$eligible_primary))=='true')&&sum(primary$primary_Holm_P>=.05)==2)
write.csv(do.call(rbind,checks),file.path(out,'logs/independent_analysis_checks.csv'),row.names=FALSE)
write_json(list(status='PASS',checks=length(checks),passed=length(checks),independent_implementation='BaseRscore reconstruction,t.test,lm,andp.adjust;full76-profile challenge model compared with donor-change model',members_unchanged=TRUE),file.path(out,'logs/independent_analysis_checks.json'),pretty=TRUE,auto_unbox=TRUE)
print(paste('PASS',length(checks),'independent checks'))
