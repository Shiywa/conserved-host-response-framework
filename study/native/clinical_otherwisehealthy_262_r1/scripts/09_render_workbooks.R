suppressPackageStartupMessages({library(openxlsx);library(jsonlite)})
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)])
O<-file.path(Sys.getenv('CRRM_WORKSPACE'),'clinical_otherwisehealthy_262_r1');D<-file.path(O,'submission_tables')
spec<-fromJSON(file.path(O,'parameters/submission_workbooks.json'),simplifyVector=FALSE)
titleStyle<-createStyle(fontSize=13,textDecoration='bold',fontColour='#173F59',wrapText=TRUE,valign='center')
headStyle<-createStyle(fontSize=10,fgFill='#28556C',fontColour='#FFFFFF',textDecoration='bold',wrapText=TRUE,valign='center')
bodyStyle<-createStyle(fontSize=10,wrapText=TRUE,valign='top')
checks<-list()
for(file in unique(vapply(spec,function(x)x$workbook,character(1)))){
 ss<-Filter(function(x)x$workbook==file,spec)
 wb<-createWorkbook(creator='Liu et al.');modifyBaseFont(wb,fontName='Arial',fontSize=10)
 for(v in ss){
  d<-read.csv(file.path(D,v$file),colClasses='character',check.names=FALSE,na.strings=c('','NA'),stringsAsFactors=FALSE)
  numeric<-names(d)[vapply(d,function(x)all(is.na(x)|grepl('^-?[0-9]+(\\.[0-9]*)?([eE][+-]?[0-9]+)?$',x)),logical(1))]
  for(n in numeric)d[[n]]<-as.numeric(d[[n]])
  addWorksheet(wb,v$sheet,gridLines=FALSE);writeData(wb,v$sheet,v$title,startRow=1,colNames=FALSE)
  mergeCells(wb,v$sheet,rows=1,cols=seq_len(ncol(d)));addStyle(wb,v$sheet,titleStyle,rows=1,cols=1);setRowHeights(wb,v$sheet,1,34)
  writeDataTable(wb,v$sheet,d,startRow=2,tableName=gsub('_','',v$sheet),headerStyle=headStyle,tableStyle='TableStyleLight9',keepNA=FALSE)
  addStyle(wb,v$sheet,bodyStyle,rows=3:(nrow(d)+2),cols=seq_len(ncol(d)),gridExpand=TRUE,stack=TRUE)
  setRowHeights(wb,v$sheet,2,38);setRowHeights(wb,v$sheet,3:(nrow(d)+2),if(grepl('Summary|Key_results',v$sheet))58 else 25)
  setColWidths(wb,v$sheet,seq_len(ncol(d)),18)
  if(grepl('Summary|Key_results',v$sheet)){
   setColWidths(wb,v$sheet,1:6,c(38,23,62,46,45,18));setColWidths(wb,v$sheet,c(10,13,16),c(44,49,86))
   for(n in c('Estimate','Lower bound','Upper bound'))addStyle(wb,v$sheet,createStyle(numFmt='0.0000'),rows=3:(nrow(d)+2),cols=match(n,names(d)),stack=TRUE)
   for(n in c('Model P','Adjusted P','Sensitivity P','Sensitivity adjusted P'))addStyle(wb,v$sheet,createStyle(numFmt='0.000E+00'),rows=3:(nrow(d)+2),cols=match(n,names(d)),stack=TRUE)
  }else{
   setColWidths(wb,v$sheet,1:3,c(19,26,19));setColWidths(wb,v$sheet,4:ncol(d),25)
   addStyle(wb,v$sheet,createStyle(numFmt='0.0000'),rows=3:(nrow(d)+2),cols=4:ncol(d),gridExpand=TRUE,stack=TRUE)
  }
  freezePane(wb,v$sheet,firstActiveRow=3,firstActiveCol=3)
  pageSetup(wb,v$sheet,orientation='landscape',paperSize=9,fitToWidth=1,fitToHeight=0,printTitleRows=1:2)
 }
 dest<-file.path(D,file);saveWorkbook(wb,dest,overwrite=TRUE)
 stopifnot(identical(getSheetNames(dest),vapply(ss,function(x)x$sheet,character(1))))
 for(v in ss){
  expected<-read.csv(file.path(D,v$file),colClasses='character',check.names=FALSE,na.strings=c('','NA'),stringsAsFactors=FALSE)
  actual<-read.xlsx(dest,sheet=v$sheet,startRow=2,check.names=FALSE,sep.names=' ')
  stopifnot(identical(names(actual),names(expected)),nrow(actual)==nrow(expected))
  for(n in names(expected)){
   e<-expected[[n]];aa<-actual[[n]]
   if(is.numeric(aa)){
    stopifnot(isTRUE(all.equal(as.numeric(e),aa,tolerance=1e-12,check.attributes=FALSE)))
   }else{
    stopifnot(identical(ifelse(is.na(e),'',e),ifelse(is.na(aa),'',as.character(aa))))
   }
  }
  checks[[length(checks)+1]]<-list(file=file,sheet=v$sheet,rows=nrow(actual),typed_readback='PASS')
 }
}
# Preserve the established 80-row S7 table; only the clinical-correspondence block changes.
original<-file.path(Sys.getenv('CRRM_DATA_ROOT'),"\u5f53\u524d\u7a3f\u4ef6/\u8865\u5145\u8868/TableS7_\u5bbf\u4e3b\u5e72\u9884.xlsx")
wb<-loadWorkbook(original);sh<-getSheetNames(original);stopifnot(length(sh)==1)
d<-read.xlsx(original,sheet=sh,startRow=2,check.names=FALSE,sep.names=' ')
al<-read.csv(file.path(O,'host/results/clinical_gene_alignment.csv'),check.names=FALSE)
stopifnot(nrow(d)==80)
ix<-which(d$Section=='Clinical alignment');stopifnot(length(ix)==3)
for(i in ix){
 r<-al[al$clinical_contrast==d$Context[i],];stopifnot(nrow(r)==1)
 for(k in c('Estimate','BH P','Exact P')){
  value<-if(k=='Estimate')r$rho else if(k=='BH P')r$adjusted_P_BH3 else r$permutation_P
  writeData(wb,sh,value,startRow=i+2,startCol=match(k,names(d)),colNames=FALSE)
 }
 if(r$clinical_contrast=='GSE246622_Severe_Mild')writeData(wb,sh,paste(r$members,'measured members of the fixed 25-gene annotation; unadjusted Severe–Mild comparison (36/79 expression profiles); 70 culture-label allocations; BH3'),startRow=i+2,startCol=match('Notes',names(d)),colNames=FALSE)
}
dest<-file.path(D,'TableS7_宿主干预.xlsx');saveWorkbook(wb,dest,overwrite=TRUE)
updated<-read.xlsx(dest,sheet=sh,startRow=2,check.names=FALSE,sep.names=' ')
stopifnot(nrow(updated)==80,identical(d[-ix,],updated[-ix,]))
for(i in ix){r<-al[al$clinical_contrast==updated$Context[i],];stopifnot(abs(updated$Estimate[i]-r$rho)<1e-13,abs(updated[['BH P']][i]-r$adjusted_P_BH3)<1e-13)}
checks[[length(checks)+1]]<-list(file='TableS7_宿主干预.xlsx',sheet=sh,rows=80,unchanged_nonclinical_rows=77,typed_readback='PASS')
write_json(list(status='PASS',workbooks=checks),file.path(O,'audit/workbook_readback.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(O,'audit/workbook_R_sessionInfo.txt'))
print(checks)
