library(data.table)

dir.create("CHPLOTs",showWarnings=F)

# FUNCTION :
# >>>>>>>>>>
pep.compMIRdt <-
 function(
  RES1,RES2,MRBY="ID",SUFF=c(".res1",".res2"),
  COL=c("gray50","black","black","black","black","black","black"),
  BG =c("gray50","blue","red","limegreen",rgb(0,0,1,0.3,maxColorValue=1),rgb(1,0,0,0.3,maxColorValue=1),rgb(50,205,50,0.3,maxColorValue=250)),
  CEX=c(0.75,1.25,1.50),PCH=c(16,25,24),
  PCL="padj",PTH=0.05,PTP="type",PST="miRNA.STR",PMR="mirtron",xLIM=NULL,yLIM=NULL,xSCALE=c(-0.2,+0.2),ySCALE=c(-0.2,+0.2),xySQ=F,
  TIT=NULL,xLAB=NULL,yLAB=NULL,MIRTRON=T,xTYPE="all",
  pLEG="topleft",pSTAT="topright",tLEG=NULL,no.LEGEND=F,xRpos=NULL,yRpos=NULL,xRadd=+0.00,yRadd=-0.50
 ) {

 cat( "\n miRNA type : \"",xTYPE,"\"\n",sep="")
 cat( " Acceptable inputs : \"miRNA\", \"miRNA*\", \"all\"\n",sep="" )

 if ( is.null(tLEG) ) {
  tLEG <-
   paste(
    c(paste0("only ",deparse(substitute(RES1))),paste0("only ",deparse(substitute(RES2))),"both"),
    c(" sig."," sig."," sig."),
    sep=""
   )
  tLEG <- sprintf("%21s",tLEG)
 }

 plot.dt <-
  merge(
   RES1[,c(MRBY,"log2FoldChange",PCL),with=F],
   RES2[,c(MRBY,"log2FoldChange",PCL,PTP,PST,PMR),with=F],
   by=MRBY,suffixes=SUFF,all=T,sort=F)
 if ( xTYPE!="all" ) {
 plot.dt <- plot.dt[ get(PTP)==xTYPE ]
 }

 if ( is.null(xLIM) ) {
  xLIM <- range(plot.dt[[paste("log2FoldChange",SUFF[1],sep="")]], na.rm=T) ; xLIM <- xLIM + (abs(xLIM)*xSCALE)
 }
 if ( is.null(yLIM) ) {
  yLIM <- range(plot.dt[[paste("log2FoldChange",SUFF[2],sep="")]], na.rm=T) ; yLIM <- yLIM + (abs(yLIM)*ySCALE)
 }
 if ( is.null(xLAB) ) {
  xLAB <- paste("log2FoldChange",SUFF[1],sep="")
 }
 if ( is.null(yLAB) ) {
  yLAB <- paste("log2FoldChange",SUFF[2],sep="")
 }
 if ( is.null(xRpos) ) {
  xRpos <- xLIM[1]*0.75
 }
 if ( is.null(yRpos) ) {
  yRpos <- yLIM[2]*0.75
 }
 if ( no.LEGEND ) { tLEG <- NULL }
 if ( xySQ ) {
  sLIM <- c(min(xLIM[1],yLIM[1]),max(xLIM[2],yLIM[2]))
  xLIM <- sLIM
  yLIM <- sLIM
 }

 cat(" X scale : < ",xLIM[1],";",xLIM[2]," >\n"," Y scale : < ",yLIM[1],";",yLIM[2]," >\n",sep="")

 corFmt <- function(COR) {
  if ( !is.na(COR) ) {
   if ( COR>0 ) { COR<-paste("+",sprintf("%5.3f",COR),sep="") } else { COR<-paste(    sprintf("%5.3f",COR),sep="") }
  }
  return(COR)
 }

#= not-sig.
 plot(
  plot.dt[ ( get(paste0(PCL,SUFF[1]))>PTH ) & ( get(paste0(PCL,SUFF[2]))>PTH ) ,paste("log2FoldChange",SUFF,sep=""),with=F],
  pch=PCH[1],col=COL[1],bg=BG[1],cex=CEX[1],main=TIT,xlab=xLAB,ylab=yLAB,xlim=xLIM,ylim=yLIM,
  panel.first=abline(h=0,v=0,col="black",lty=1,lwd=2),font.lab=2
 )

#= RES1 sig only
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="5p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="5p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[2],cex=CEX[2],bg=BG[5],col=COL[5]
 )
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="3p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="3p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[3],cex=CEX[2],bg=BG[6],col=COL[6]
 )
 if ( MIRTRON ) {
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="5p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="5p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[2],cex=CEX[2],bg=BG[7],col=COL[7]
 )
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="3p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))> PTH ) & ( get(PST)=="3p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[3],cex=CEX[2],bg=BG[7],col=COL[7]
 )
 }
#= RES2 sig only
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[2],cex=CEX[2],bg=BG[5],col=COL[5]
 )
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[3],cex=CEX[2],bg=BG[6],col=COL[6]
 )
 if ( MIRTRON ) {
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[2],cex=CEX[2],bg=BG[7],col=COL[7]
 )
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))> PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[3],cex=CEX[2],bg=BG[7],col=COL[7]
 )
 }
#= both sig
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[2],cex=CEX[3],bg=BG[2],col=COL[2]
 )
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)!=T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[3],cex=CEX[3],bg=BG[3],col=COL[3]
 )
 if ( MIRTRON ) {
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="5p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[2],cex=CEX[3],bg=BG[4],col=COL[4]
 )
 points(
  x=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[1],sep="")]],
  y=plot.dt[ ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) & ( get(PST)=="3p" ) & ( get(PMR)==T ) ][[paste("log2FoldChange",SUFF[2],sep="")]],
  pch=PCH[3],cex=CEX[3],bg=BG[4],col=COL[4]
 )
 }

 cor1 <-
  cor(
   plot.dt[ !is.na(get(paste0(PCL,SUFF[1]))) & !is.na(get(paste0(PCL,SUFF[2])))                                                                         ,paste("log2FoldChange",SUFF,sep=""),with=F ]
  )[1,2]
 cor2 <-
  cor(
   plot.dt[ !is.na(get(paste0(PCL,SUFF[1]))) & !is.na(get(paste0(PCL,SUFF[2]))) & ( get(paste0(PCL,SUFF[1]))<=PTH ) & ( get(paste0(PCL,SUFF[2]))<=PTH ) ,paste("log2FoldChange",SUFF,sep=""),with=F ]
  )[1,2]
 text(
  x=xRpos,      y=yRpos,      labels=paste("R = ",corFmt(COR=cor1),sep=""),col="black"
 )
 text(
  x=xRpos+xRadd,y=yRpos+yRadd,labels=paste("R = ",corFmt(COR=cor2),sep=""),col="red"
 )
 if ( !no.LEGEND ) {
 legend( pLEG,pch=PCH[2:3],pt.bg=BG[2:3],cex=CEX[c(2,2)],bty="n",legend=c("5'","3'") )
 legend(
  pSTAT,bty="n",
  legend=paste0(
   tLEG," [",
   c(
    nrow(plot.dt[(get(paste(PCL,SUFF[1],sep=""))<=PTH)&(get(paste(PCL,SUFF[2],sep=""))> PTH)]),
    nrow(plot.dt[(get(paste(PCL,SUFF[1],sep=""))> PTH)&(get(paste(PCL,SUFF[2],sep=""))<=PTH)]),
    nrow(plot.dt[(get(paste(PCL,SUFF[1],sep=""))<=PTH)&(get(paste(PCL,SUFF[2],sep=""))<=PTH)])
   ),"]"
  )
 )
 }

 return( plot.dt )

}
# <<<<<<<<<<

pdf( file=paste0("CHPLOTs/dYEL-vs-dGRN-",format(Sys.time(), "%Y%m%d",".pdf"),width=15,height=5 )
 par( bg="white",mfrow=c(1,3) )

  pep.compMIRdt(
   RES1=dYEL.dt,RES2=dGRN.dt,
   MIRTRON=T,pLEG="topleft",pSTAT="bottomright",tLEG=c( "only dYEL sig.","only dGRN sig.","both sig." ),
   SUFF=c( "_dYEL","_dGRN" ),
   xTYPE="all",TIT="dYEL vs. dGRN (all miRNAs)"
  )
  pep.compMIRdt(
   RES1=dYEL.dt,RES2=dGRN.dt,
   MIRTRON=T,pLEG="topleft",pSTAT="bottomright",tLEG=c( "only dYEL sig.","only dGRN sig.","both sig." ),
   SUFF=c( "_dYEL","_dGRN" ),
   xTYPE="miRNA",TIT="dYEL vs. dGRN (main strand)"
  )
  pep.compMIRdt(
   RES1=dYEL.dt,RES2=dGRN.dt,
   MIRTRON=T,pLEG="topleft",pSTAT="bottomright",tLEG=c( "only dYEL sig.","only dGRN sig.","both sig." ),
   SUFF=c( "_dYEL","_dGRN" ),
   xTYPE="miRNA*",TIT="dYEL vs. dGRN (passengers)"
  )

dev.off()
