suppressPackageStartupMessages({library(openxlsx);library(jsonlite)})
a<-commandArgs(trailingOnly=FALSE);s<-sub('^--file=','',a[grepl('^--file=',a)]);out<-file.path(Sys.getenv('CRRM_WORKSPACE'),'phase18_Result2_external_response_20261112');spec<-fromJSON(file.path(out,'parameters/workbook_spec.json'),simplifyVector=FALSE)
wb<-createWorkbook(creator='Liu et al.',title='Supplementary Table S2: key external CRRM results');modifyBaseFont(wb,fontSize=10,fontName='Arial',fontColour='#253342')
title<-createStyle(fontName='Arial',fontSize=12,textDecoration='bold',fontColour='#173F59',valign='center',wrapText=TRUE);header<-createStyle(fontName='Arial',fontSize=10,textDecoration='bold',fontColour='#FFFFFF',fgFill='#28556C',wrapText=TRUE,valign='center');body<-createStyle(fontName='Arial',fontSize=10,valign='top',wrapText=TRUE);notes<-createStyle(fontName='Arial',fontSize=10,fontColour='#52616B',wrapText=TRUE,valign='top');dec<-createStyle(numFmt='0.000');pv<-createStyle(numFmt='0.00E+00');manifest<-list()
for(s in spec){
 d<-read.csv(file.path(out,'tables',s$file),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA',''));sh<-s$sheet;nc<-ncol(d);nr<-nrow(d);addWorksheet(wb,sh,gridLines=FALSE,zoom=85)
 writeData(wb,sh,s$title,startRow=1,colNames=FALSE);mergeCells(wb,sh,cols=seq_len(nc),rows=1);addStyle(wb,sh,title,rows=1,cols=1);setRowHeights(wb,sh,1,32)
 writeDataTable(wb,sh,d,startRow=2,tableName=sh,tableStyle='TableStyleLight9',withFilter=TRUE,keepNA=TRUE,na.string='NA',headerStyle=header);addStyle(wb,sh,body,rows=3:(nr+2),cols=seq_len(nc),gridExpand=TRUE,stack=TRUE)
 setRowHeights(wb,sh,2,40);setRowHeights(wb,sh,3:(nr+2),if(sh=='S2B_Gene_directions')18 else 32);setColWidths(wb,sh,seq_len(nc),17)
 textcols<-which(vapply(d,is.character,logical(1)));if(length(textcols))setColWidths(wb,sh,textcols,30)
 if(sh=='S2B_Gene_directions')setColWidths(wb,sh,1:3,c(18,24,18))
 if(sh=='S2D_Context_and_components')setColWidths(wb,sh,c(5,12),c(44,64))
 vals<-which(grepl('score change|log-expression effect|95% CI',names(d)));pcols<-which(grepl(' P|^Model P',names(d)));if(length(vals))addStyle(wb,sh,dec,rows=3:(nr+2),cols=vals,gridExpand=TRUE,stack=TRUE);if(length(pcols))addStyle(wb,sh,pv,rows=3:(nr+2),cols=pcols,gridExpand=TRUE,stack=TRUE)
 for(i in seq_along(s$notes)){rr<-nr+3+i;writeData(wb,sh,s$notes[[i]],startRow=rr,colNames=FALSE);mergeCells(wb,sh,cols=seq_len(nc),rows=rr);addStyle(wb,sh,notes,rows=rr,cols=1);setRowHeights(wb,sh,rr,32)}
 freezePane(wb,sh,firstActiveRow=3,firstActiveCol=if(sh=='S2B_Gene_directions')3 else 2);pageSetup(wb,sh,orientation='landscape',paperSize=9,fitToWidth=1,fitToHeight=0);setHeaderFooter(wb,sh,footer=c('Liu et al.','Supplementary Table S2','Page &P of &N'));manifest[[length(manifest)+1]]<-data.frame(sheet=sh,title=s$title,n_rows=nr,n_columns=nc,title_row=1,header_row=2,first_data_row=3)
}
dest<-file.path(out,'reports/SupplementaryTableS2_CRRM_key_results.xlsx');saveWorkbook(wb,dest,overwrite=TRUE);manifest<-do.call(rbind,manifest);write.csv(manifest,file.path(out,'logs/workbook_sheet_manifest.csv'),row.names=FALSE)
for(s in spec){d<-read.csv(file.path(out,'tables',s$file),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA',''));x<-read.xlsx(dest,sheet=s$sheet,rows=2:(nrow(d)+2),check.names=FALSE,sep.names=' ',na.strings='NA');stopifnot(identical(names(d),names(x)),nrow(d)==nrow(x));for(n in names(d)){if(is.numeric(d[[n]]))stopifnot(isTRUE(all.equal(d[[n]],x[[n]],tolerance=1e-12,check.attributes=FALSE)))else stopifnot(identical(as.character(d[[n]]),as.character(x[[n]])))}}
write_json(list(status='PASS',sheets=4,typed_readback=TRUE,row_counts=as.list(setNames(manifest$n_rows,manifest$sheet))),file.path(out,'logs/workbook_readback.json'),pretty=TRUE,auto_unbox=TRUE);print(manifest)
