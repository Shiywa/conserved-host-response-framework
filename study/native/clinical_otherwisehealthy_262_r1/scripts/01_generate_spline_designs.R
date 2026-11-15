suppressPackageStartupMessages({library(splines);library(jsonlite)})
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)])
O<-file.path(Sys.getenv('CRRM_WORKSPACE'),'clinical_otherwisehealthy_262_r1');P<-file.path(O,'main')
m<-read.csv(file.path(P,'data/processed/GSE246622_metadata.csv'),check.names=FALSE)
m$sex<-factor(m$sex);strict<-m[m$group%in%c('Mild','Severe'),]
strict<-strict[is.finite(strict$age_days)&!is.na(strict$sex),];strict$high<-as.numeric(strict$group=='Severe')
stopifnot(nrow(strict)==115,sum(strict$high)==36)
records<-list()
for(kind in c('unadjusted','adjusted')){
 X<-model.matrix(if(kind=='adjusted')~high+ns(age_days,df=3)+sex else ~high,strict)
 stopifnot(qr(X)$rank==ncol(X))
 write.csv(data.frame(sample_id=strict$sample_id,X,check.names=FALSE),file.path(P,'models',paste0('GSE246622_',kind,'_design.csv')),row.names=FALSE)
}
records$severity<-list(n=nrow(strict),population='115complete-case Severe/Mild profiles',df=3,
 knots=attr(ns(strict$age_days,df=3),'knots'),boundary_knots=attr(ns(strict$age_days,df=3),'Boundary.knots'),
 inherited_between_exploratory_folds=TRUE,foldwise_spline_refit=FALSE)
m$group<-factor(m$group,levels=c('Healthy','Mild','Moderate','Severe'))
X<-model.matrix(~group+ns(age_days,df=3)+sex,m);stopifnot(nrow(X)==262,qr(X)$rank==ncol(X))
write.csv(data.frame(sample_id=m$sample_id,X,check.names=FALSE),file.path(P,'models/GSE246622_gene_comparison_design.csv'),row.names=FALSE)
records$gene_effects<-list(n=nrow(m),population='262retained profiles in four clinical groups',df=3,
 knots=attr(ns(m$age_days,df=3),'knots'),boundary_knots=attr(ns(m$age_days,df=3),'Boundary.knots'))
write_json(records,file.path(O,'parameters/age_spline_parameters.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(O,'audit/spline_R_sessionInfo.txt'))
print(records)
