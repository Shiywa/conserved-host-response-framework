options(digits=17)
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)])
O<-file.path(Sys.getenv('CRRM_WORKSPACE'),'clinical_otherwisehealthy_262_r1');A<-file.path(O,'audit')
values<-read.delim(file.path(A,'affected_P_value_pairs.tsv'),check.names=FALSE)
checks<-lapply(split(values,values$family_id),function(d){
 p<-d$P;ok<-is.finite(p);p[!ok]<-1;q<-p.adjust(p,method=d$method[1],n=d$n[1]);q[!ok]<-NA_real_
 good<-is.finite(q)&is.finite(d$stored_adjusted_P)
 err<-if(any(good))max(abs(q[good]-d$stored_adjusted_P[good]))else 0
 stopifnot(identical(is.finite(q),is.finite(d$stored_adjusted_P)),err<1e-10)
 data.frame(family=d$family_id[1],n=d$n[1],positions=nrow(d),max_error=err)
});checks<-do.call(rbind,checks)
write.table(checks,file.path(A,'R_adjustment_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

d<-read.delim(file.path(A,'HC3_models_for_R.tsv'),check.names=FALSE)
models<-lapply(seq_len(nrow(d)),function(i){
 z<-read.csv(d$file[i],check.names=FALSE);X<-as.matrix(z[,setdiff(names(z),c('sample_id','outcome_score')),drop=FALSE]);y<-z$outcome_score
 fit<-lm(y~0+X);j<-which(colnames(X)=='high');res<-residuals(fit);h<-hatvalues(fit);bread<-solve(crossprod(X));cov<-bread%*%crossprod(X*(res/(1-h)))%*%bread
 e<-unname(coef(fit)[j]);se<-sqrt(cov[j,j]);df<-df.residual(fit);q<-qt(.975,df);vv<-c(e,e-q*se,e+q*se,2*pt(-abs(e/se),df))
 old<-as.numeric(d[i,c('effect','CI95_low','CI95_high','p_value')]);err<-max(abs(vv-old));stopifnot(err<1e-8)
 data.frame(label=d$label[i],n=nrow(X),df=df,max_error=err)
});models<-do.call(rbind,models)
write.table(models,file.path(A,'R_HC3_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
d<-read.delim(file.path(A,'Welch_models_for_R.tsv'),check.names=FALSE)
wt<-lapply(seq_len(nrow(d)),function(i){
 z<-read.csv(d$file[i]);tt<-t.test(z$score[z$high==1],z$score[z$high==0],var.equal=FALSE)
 vv<-c(unname(diff(rev(tt$estimate))),tt$conf.int,tt$p.value)
 old<-as.numeric(d[i,c('effect','CI95_low','CI95_high','p_value')]);err<-max(abs(vv-old));stopifnot(err<1e-8)
 data.frame(label=d$label[i],max_error=err)
});wt<-do.call(rbind,wt)
write.table(wt,file.path(A,'R_Welch_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

# Independent Spearman + all70whole-culture label assignments for the updated bridge.
z<-read.delim(file.path(A,'S7E_eight_culture_expression.tsv'),check.names=FALSE)
v<-read.delim(file.path(A,'S7E_whole_blood_vector.tsv'),check.names=FALSE)
X<-as.matrix(z[,v$human_ensembl,drop=FALSE]);y<-v$clinical_log_effect
eff<-colMeans(X[z$BRD4i==1,,drop=FALSE])-colMeans(X[z$BRD4i==0,,drop=FALSE]);rho<-cor(eff,y,method='spearman')
assign<-combn(8,4);null<-apply(assign,2,function(ix)cor(colMeans(X[ix,,drop=FALSE])-colMeans(X[-ix,,drop=FALSE]),y,method='spearman'))
pp<-mean(abs(null)>=abs(rho)-1e-12)
old<-read.csv(file.path(O,'host/results/clinical_gene_alignment.csv'));old<-old[old$clinical_contrast=='GSE246622_Severe_Mild',]
stopifnot(abs(rho-old$rho)<1e-12,abs(pp-old$permutation_P)<1e-12)
write.table(data.frame(rho=rho,permutation_P=pp,allocations=ncol(assign)),file.path(A,'R_S7E_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(A,'independent_R_sessionInfo.txt'))
cat('PASS:',nrow(checks),'families;',sum(checks$positions),'P positions;',nrow(models),'HC3 fits;',nrow(wt),'Welch contrasts;70culture allocations\n')
cat('Largest errors:adjustment',max(checks$max_error),'HC3',max(models$max_error),'Welch',max(wt$max_error),'\n')
