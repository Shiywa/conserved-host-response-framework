suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)]);OUT<-file.path(Sys.getenv('CRRM_WORKSPACE'),'Result6_analysis')
read<-function(name)read.csv(file.path(OUT,name),check.names=FALSE,stringsAsFactors=FALSE)
e<-read('results/score_effects.csv');u<-read('data/inference_units.csv')
replay<-function(a,b){
 a<-a[is.finite(a)];b<-b[is.finite(b)]
 if(min(length(a),length(b))<3)return(c(estimate=NA,lower=NA,upper=NA,P=NA))
 if(diff(range(a))<=1e-12 && diff(range(b))<=1e-12)return(c(estimate=mean(a)-mean(b),lower=NA,upper=NA,P=1))
 t<-t.test(a,b,var.equal=FALSE)
 c(estimate=unname(t$estimate[1]-t$estimate[2]),lower=t$conf.int[1],upper=t$conf.int[2],P=t$p.value)
}
spl<-split(u,u$test_id);rr<-lapply(seq_len(nrow(e)),function(i){
 z<-spl[[e$test_id[i]]]
 if(is.null(z))return(data.frame(test_id=e$test_id[i],estimate=NA,lower=NA,upper=NA,P=NA))
 q<-replay(z$value[z$role=='case'],z$value[z$role=='control'])
 data.frame(test_id=e$test_id[i],as.list(q),check.names=FALSE)
});r<-do.call(rbind,rr)
stopifnot(identical(e$test_id,r$test_id))
maxerr<-max(abs(as.matrix(e[,c('estimate','lower','upper','P')])-as.matrix(r[,c('estimate','lower','upper','P')])),na.rm=TRUE)
stopifnot(maxerr<1e-9)
qh<-function(p,n){v<-p;v[!is.finite(v)]<-1;q<-p.adjust(v,'BH',n=n);q[!is.finite(p)]<-NA;q}
r$adjusted_P<-NA
for(ix in split(seq_len(nrow(e)),e$family))r$adjusted_P[ix]<-qh(r$P[ix],e$family_size[ix[1]])
qerr<-max(abs(r$adjusted_P-e$adjusted_P),na.rm=TRUE);stopifnot(qerr<1e-9)

g<-read('results/intervention_gene_effects.csv');gu<-read('data/gene_inference_units.csv.gz');sp<-split(gu,gu$test_id)
gr<-lapply(seq_len(nrow(g)),function(i){
 z<-sp[[g$test_id[i]]]
 if(is.null(z))return(c(estimate=NA,lower=NA,upper=NA,P=NA))
 replay(z$value[z$role=='case'],z$value[z$role=='control'])
});gr<-as.data.frame(do.call(rbind,gr));gerr<-max(abs(as.matrix(g[,c('estimate','lower','upper','P')])-as.matrix(gr)),na.rm=TRUE);stopifnot(gerr<1e-8)
gqerr<-max(abs(qh(gr$P,4656)-g$adjusted_P_BH4656),na.rm=TRUE);stopifnot(gqerr<1e-8)

# Reconstruct primary human scores directly from deposited TPM and original signs.
x<-read('data/human_author_TPM.csv.gz');core<-read('parameters/frozen_CRRM388.csv');hm<-read('data/human_sample_manifest.csv');ps<-read('results/profile_scores.csv')
xx<-log2(1+t(as.matrix(x[match(core$human_ensembl,x$gene_id),hm$column])))
sdx<-apply(xx,2,sd);keep<-is.finite(sdx)&sdx>1e-8
z<-sweep(sweep(xx[,keep,drop=FALSE],2,colMeans(xx[,keep,drop=FALSE]),'-'),2,sdx[keep],'/')
z<-sweep(z,2,ifelse(core$frozen_direction[keep]=='up',1,-1),'*');total<-rowMeans(z)
pp<-ps[ps$source=='BRD4_hSAEC' & ps$method=='Primary',];pp<-pp[match(hm$unit,pp$unit),]
serr<-max(abs(total-pp$CRRM388));stopifnot(serr<1e-10)

# Replicate vector-level clinical alignment with whole-treatment-label permutations.
al<-read('results/clinical_gene_alignment.csv');pts<-read('results/clinical_gene_alignment_points.csv')
tx<-which(hm$group%in%c('RSV','RSV + BRD4i'));asgn<-combn(seq_along(tx),4)
alignment<-lapply(seq_len(nrow(al)),function(i){
 q<-pts[pts$clinical_contrast==al$clinical_contrast[i],]
 gx<-log2(1+t(as.matrix(x[match(q$human_ensembl,x$gene_id),hm$column[tx]])))
 observed<-colMeans(gx[hm$group[tx]=='RSV + BRD4i',,drop=FALSE])-colMeans(gx[hm$group[tx]=='RSV',,drop=FALSE])
 rho<-cor(observed,q$clinical_log_effect,method='spearman')
 null<-apply(asgn,2,function(ii)cor(colMeans(gx[ii,,drop=FALSE])-colMeans(gx[-ii,,drop=FALSE]),q$clinical_log_effect,method='spearman'))
 pe<-mean(abs(null)>=abs(rho)-1e-12)
 c(rho=rho,permutation_P=pe)
});ar<-do.call(rbind,alignment);aerr<-max(abs(ar-as.matrix(al[,c('rho','permutation_P')])));stopifnot(aerr<1e-10)

di<-read('results/mouse_identity_composition_effects.csv');du<-read('results/mouse_epithelial_composition_scores.csv')
dr<-lapply(seq_len(nrow(di)),function(i){
 z<-du[du$method==di$method[i],];t<-di$score[i]
 replay(z[z$group=='RSV RelA CKO',t],z[z$group=='RSV vehicle',t])
});dr<-as.data.frame(do.call(rbind,dr));derr<-max(abs(as.matrix(di[,c('estimate','lower','upper','P')])-as.matrix(dr)),na.rm=TRUE)
dqerr<-max(abs(qh(dr$P,9)-di$adjusted_P_BH9));stopifnot(derr<1e-10,dqerr<1e-10)
checks<-list(status='PASS',score_tests_checked=sum(is.finite(r$P)),score_family_count=length(unique(e$family)),max_score_statistic_error=maxerr,max_score_BH_error=qerr,gene_tests_checked=sum(is.finite(gr$P)),gene_family_denominator=4656,max_gene_statistic_error=gerr,max_gene_BH_error=gqerr,human_primary_score_reconstruction_error=serr,clinical_vector_permutations_checked=3,max_alignment_error=aerr,annotation_diagnostic_tests=9,max_diagnostic_error=derr,max_diagnostic_BH_error=dqerr)
write_json(checks,file.path(OUT,'audit/independent_R_statistics.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat(toJSON(checks,auto_unbox=TRUE,pretty=TRUE,digits=8),'\n')
