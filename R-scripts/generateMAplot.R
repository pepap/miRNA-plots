library(data.table)

# FUNCTION :
# >>>>>>>>>>
pep.mirMAplot <- function( MIRDT,xTYPE="miRNA",TIT=NULL,lfcLIM=NULL,xLIM=NULL,xLAB=NULL,mirtron=F,COLS=c("blue","red","gray50","limegreen"),PLOT=T ) {

 cat( "\n miRNA type : \"",xTYPE,"\"\n",sep="")
 cat( " Acceptable inputs : \"miRNA\", \"miRNA*\", \"all\"\n",sep="" )

 if ( is.null(TIT) ) {
  TIT <- sub( "[*]","-passenger",x=TIT )
 }

 if ( xTYPE=="all" ) {
  MIRDT          <- MIRDT[ type %in% unique(type) ]
 } else {
  MIRDT          <- MIRDT[ type==xTYPE ]
 }
 MIRDT[["LOG"]] <- log10(MIRDT[["baseMean"]]+1)
 print( range(MIRDT[["LOG"]]) )
 print( range(MIRDT[["log2FoldChange"]]) )
 MIRDT[["BG"]]                    <- COLS[3]
 MIRDT[ miRNA.STR=="5p" ][["BG"]] <- COLS[1]
 MIRDT[ miRNA.STR=="3p" ][["BG"]] <- COLS[2]
 if ( mirtron ) {
  MIRDT[ mirtron==T ][["BG"]] <- COLS[4]
 }
 
 if ( is.null(lfcLIM) ) {
  YMAX           <- max(abs(na.omit(MIRDT[["log2FoldChange"]])))
  YLIM           <- c(YMAX*(-1.01),YMAX*(+1.01))
 } else                 {
  YMAX           <- max(    na.omit(MIRDT[["log2FoldChange"]]) )
  YMIN           <- min(    na.omit(MIRDT[["log2FoldChange"]]) )
  YLIM           <- lfcLIM
 if ( YMIN<YLIM[1] ) { cat(" ! User defined MIN value of log2FoldChange : ",YLIM[1]," !\n ! Dataset      MIN value of log2FoldChange : ",YMIN," !\n",sep="") }
 if ( YMAX>YLIM[2] ) { cat(" ! User defined MAX value of log2FoldChange : ",YLIM[2]," !\n ! Dataset      MAX value of log2FoldChange : ",YMAX," !\n",sep="") }
 }

 if ( is.null(xLIM) ) {
  XLIM           <- range(MIRDT[["LOG"]])*c(0.99,1.01)
 } else               {
  XMAX           <- max(    na.omit(MIRDT[["LOG"]]) )
  XMIN           <- min(    na.omit(MIRDT[["LOG"]]) )
  XLIM           <- xLIM
 if ( XMIN<XLIM[1] ) { cat(" ! User defined MIN value of baseMean : ",      XLIM[1]," !\n ! Dataset      MIN value of baseMean : ",      XMIN," !\n",sep="") }
 if ( XMAX>XLIM[2] ) { cat(" ! User defined MAX value of baseMean : ",      XLIM[2]," !\n ! Dataset      MAX value of baseMean : ",      XMAX," !\n",sep="") }
 }

 if ( is.null(xLAB) ) {
  XLAB           <- seq( round(XLIM[1]),round(XLIM[2]),1 )
 } else {
  XLAB           <- xLAB
 }

 MIRDT[["PCH"]] <- 16
 MIRDT[["CEX"]] <- 1.50
 if ( nrow(MIRDT[ miRNA.STR=="5p" & padj<=0.05 ])!=0 ) {
  MIRDT[ miRNA.STR=="5p" & padj<=0.05 ][["PCH"]] <- 25
  MIRDT[ miRNA.STR=="5p" & padj<=0.05 ][["CEX"]] <- 2.25
 }
 if ( nrow(MIRDT[ miRNA.STR=="3p" & padj<=0.05 ])!=0 ) {
  MIRDT[ miRNA.STR=="3p" & padj<=0.05 ][["PCH"]] <- 24
  MIRDT[ miRNA.STR=="3p" & padj<=0.05 ][["CEX"]] <- 2.25
 }

 if (PLOT) {
  nonsig.PCH=16
  nonsig.CEX=1.50
  nonsig.COL=COLS[3]
  nonsig.BG ="transparent"
  par( bg="white" )
  plot(
   MIRDT[,c("LOG","log2FoldChange"),with=F],
   col=nonsig.COL,bg=nonsig.BG,pch=nonsig.PCH,cex=nonsig.CEX,
   xlab="mean of normalized counts (log10)",ylab="log2FoldChange",main=TIT,xlim=XLIM,ylim=YLIM,
   xaxt="n",
   font.lab=2
  )
  axis( side=1,at=seq(XLIM[1],XLIM[2],length.out=length(XLAB)),labels=XLAB )
  sig5.PCH=25
  sig5.CEX=2.25
  sig5.COL="transparent"
  sig5.BG =COLS[1]
  points(
   MIRDT[ miRNA.STR=="5p" & padj<=0.05 ,c("LOG","log2FoldChange"),with=F],
   pch=sig5.PCH,cex=sig5.CEX,col=sig5.COL,bg=MIRDT[ miRNA.STR=="5p" & padj<=0.05 , BG ]
  )
  sig3.PCH=24
  sig3.CEX=2.25
  sig3.COL="transparent"
  sig3.BG =COLS[2]
  points(
   MIRDT[ miRNA.STR=="3p" & padj<=0.05 ,c("LOG","log2FoldChange"),with=F],
   pch=sig3.PCH,cex=sig3.CEX,col=sig3.COL,bg=MIRDT[ miRNA.STR=="3p" & padj<=0.05 , BG ]
  )
  legend(
   "topleft",
   pch=c(sig5.PCH,sig3.PCH),
   legend=c("5'","3'"),
   bty="n",
   col=c(sig5.COL,sig3.COL),
   pt.bg=c(sig5.BG,sig3.BG),
   pt.cex=c(sig5.CEX,sig3.CEX),
   y.intersp=2
  )
 }

 return( MIRDT )

}
# <<<<<<<<<<

dir.create("MAPLOTs",showWarnings=F)

x.lim=c(0.0,7.0)
x.lab=c("0","10^1","10^2","10^3","10^4","10^5","10^6","10^7")
y.lim=c(-6.50,+4.00)

for ( GNTP in c( "dYEL","dGRN" ) ) {
 pdf( file=paste0("MAPLOTs/MAplot-",GNTP,"-",format(Sys.time(), "%Y%m%d"),".pdf"),  width=50,height=15,pointsize=30 )
  par( bg="white",mfrow=c(1,3) )
  
  otmp_data <-
   pep.mirMAplot(
    MIRDT   = get(paste0(GNTP,".dt")),
    xTYPE   = "all",
    TIT     = paste0(GNTP," : all miRNAs"),
    lfcLIM  = y.lim,
    xLIM    = x.lim,
    xLAB    = x.lab,
    mirtron = T
   )
  otmp_data <-
   pep.mirMAplot(
    MIRDT   = get(paste0(GNTP,".dt")),
    xTYPE   = "miRNA",
    TIT     = paste0(GNTP," : main strand"),
    lfcLIM  = y.lim,
    xLIM    = x.lim,
    xLAB    = x.lab,
    mirtron = T
   )
  otmp_data <-
   pep.mirMAplot(
    MIRDT   = get(paste0(GNTP,".dt")),
    xTYPE   = "miRNA*",
    TIT     = paste0(GNTP," : passengers"),
    lfcLIM  = y.lim,
    xLIM    = x.lim,
    xLAB    = x.lab,
    mirtron = T
   )

 dev.off()
}
