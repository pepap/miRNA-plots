library(data.table)
library(matrixStats)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(scales)

load( "FD01.mirAnnot.pepType-CPs-SA.dt.rda",verbose=T )
load( "FD02.ALL-SA-CPs.21to23nt.rda",verbose=T )

ONAME="21to23nt"

breakList <- seq(-1,+1,0.01)

source("R-scripts/selected_miRNAs.R")

IDR.VB.CON <- c(
   "SOM.r1", "SOM.r2", "SOM.r3",
  "dGRN.r1","dGRN.r2","dGRN.r3",
  "dYEL.r1","dYEL.r2","dYEL.r3"
)

BXP_REPS <-
 c(
   "SOM.r1.CPs.normCent.mat", "SOM.r2.CPs.normCent.mat", "SOM.r3.CPs.normCent.mat",
  "dGRN.r1.CPs.normCent.mat","dGRN.r2.CPs.normCent.mat","dGRN.r3.CPs.normCent.mat",
  "dYEL.r1.CPs.normCent.mat","dYEL.r2.CPs.normCent.mat","dYEL.r3.CPs.normCent.mat"
 )
BXP_CONS <-
 rep.int( c("SOM","dGRN","dYEL"),c(3,3,3) )

# FUNCTION :
# >>>>>>>>>>
# (1) matrix lines => density
normMat <- function( iMAT ) {

 out.mat <- apply(
  X      = iMAT,
  MARGIN = 1,
  FUN    = function(L){
   if ( sum(L)==0 ) {
    return(L)
   } else           {
    return(L/sum(L))
   }
  }
  )

 return( t(out.mat) )

}

# (2) define new CP based on the expression; report multiple CPs & OUTlier CPs
deffCP <- function( iMAT,SEL=as.character(seq(-10,+10)),OUT=5,ANNOT=mirAnnot.pepType.dt ) {

 all.rwn <- rownames( iMAT )
 ADD     <- as.integer(SEL[1]) - as.integer(colnames(iMAT)[1])
 CANON   <- seq_along(SEL)[ SEL=="0" ]

 j=1
 out.mat <-
  apply(
   X      = iMAT,
   MARGIN = 1,
   FUN    = function(L){
    LSEL <- L[ SEL ]
    tmp.rnw <- rownames(iMAT)[j]
    j <<- ( j + 1 )
    if ( sum(LSEL)!=0 ) {
     ZERO <- max( LSEL )
     ZERO <- seq_along(LSEL)[ LSEL==ZERO ]
     if ( length(ZERO)>1 ) {
#      print( c(ZERO,tmp.rnw) )
      return( CANON+ADD )
     }
     if ( abs(CANON-ZERO)>OUT ) {
#      print( c(ZERO,tmp.rnw) )
      return( CANON+ADD )
     } else               {
      return( ZERO +ADD )
     }
    } else              {
     return(  CANON+ADD )
    }
   }
  )

 return( out.mat[all.rwn] )

}

# (3) rename colnames in matrix, center "0" to new CP
deffCPcenter <- function( iMAT,CP.DT,EXPAND=5 ) {

 out.mat <-
  sapply(
   X   = seq(nrow(CP.DT)),
   FUN = function(i){
    cLINE        <- iMAT[ CP.DT[ i , ID ] , seq( from=CP.DT[ i , CP-EXPAND ],to=CP.DT[ i , CP+EXPAND ] ) ]
    names(cLINE) <- as.character(seq((-1)*EXPAND,(+1)*EXPAND))
    return( cLINE )
   },
   simplify=T,
   USE.NAMES=T
  )

 out.mat           <- t(out.mat)
 rownames(out.mat) <- CP.DT[ , ID ]

 return( out.mat )

}

# (4) for each genomic position return "FCE"-value of ( As - Bs )
diffMat <- function( As,Bs,FCE="median" ) {

 i=1
 tmp.ncol <- c()
 for ( a in As ) {
  for ( b in Bs ) {
   if ( i == 1 )  {
    tmp.mat <- ( get(a) - get(b) )
   } else         {
    tmp.mat <- cbind( tmp.mat,( get(a) - get(b) ) )
   }
   tmp.ncol <- append( tmp.ncol,ncol(get(a)) )
   tmp.ncol <- append( tmp.ncol,ncol(get(b)) )
   i <- ( i + 1 )
  }
 }

 tmp.ncol <- unique( tmp.ncol )
 if ( length(tmp.ncol)!=1 ) { stop("\n! Input matrices differ in number of columns !\n\n") }

 rNAME    <- rownames( tmp.mat )
 cNAME    <- c()
 out.vals <- c()
 for ( i in c( seq(1,(tmp.ncol-1)),0 ) ) {
  cNAME    <- append( cNAME,unique(colnames(tmp.mat)[ ( seq_along(tmp.mat[1,]) %% tmp.ncol ) == i ]) )
  out.vals <- append( out.vals,apply( X=tmp.mat,MARGIN=1,FUN=function(L){ return(get(FCE)( L[ ( seq_along(L) %% tmp.ncol ) == i ] )) } ) )
 }

 out.mat <- matrix( data=out.vals,byrow=F,ncol=tmp.ncol,dimnames=list(rNAME,cNAME) )

 return( out.mat )

}

# <<<<<<<<<<

#= I. define miRNA cleavage-points (CPs) from WT reads
cat( "\n ** I. define miRNA cleavage-points (CPs) from WT reads **\n\n",sep="" )
#= (1) IDRs : SOM
 IDR.CPs.dt <- data.table(
  ID=IDR.mirnas
 )
 IDR.CPs.dt <- merge( IDR.CPs.dt,mirAnnot.pepType.dt[,c("ID","gLoc","preID.gLoc","Name","det.STR","type"),with=F],by="ID",all.x=T,sort=F )

 tmp.dt <- deffCP( iMAT=get("SOM.r1.CPs.mat")[ IDR.mirnas,], SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.pepType.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),SOM.r1 =tmp.dt )
 IDR.CPs.dt <- merge( IDR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=get("SOM.r2.CPs.mat")[IDR.mirnas,], SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.pepType.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),SOM.r2=tmp.dt )
 IDR.CPs.dt <- merge( IDR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 tmp.dt <- deffCP( iMAT=get("SOM.r3.CPs.mat")[IDR.mirnas,], SEL=as.character(seq(-10,+10)),OUT=7,ANNOT=mirAnnot.pepType.dt )
 tmp.dt <- data.table( ID=names(tmp.dt),SOM.r3=tmp.dt )
 IDR.CPs.dt <- merge( IDR.CPs.dt,tmp.dt,by="ID",sort=F,all=T )

 IDR.CPs.dt[["CP"]] <- rowMedians(as.matrix(IDR.CPs.dt[,c("SOM.r1","SOM.r2","SOM.r3"),with=F]))

#= II. center CP-matrices to new CPs & normalize
cat( "\n ** II. center CP-matrices to new CPs & normalize **\n\n",sep="" )
 out.objs    <- c()
 for ( iM in paste(IDR.VB.CON,".CPs.mat",sep="") ) {

  cat( " => ",iM,"\n",sep="" )
  tmp.mat <- normMat( iMAT=deffCPcenter( iMAT=get(iM),CP.DT=IDR.CPs.dt,EXPAND=5 ) )
#  print(           dim(tmp.mat)  )
#  print( table(rowSums(tmp.mat)) )
  assign( x=sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ),value=tmp.mat )
  out.objs <- append( out.objs,sub( "[.]CPs[.]mat$",".CPs.normCent.mat",iM ) )

 }
save( list=c("IDR.CPs.dt",out.objs),file=paste("FD03.IDRs-CPs.normCent.",ONAME,".rda",sep="") )

bxp.dt <- data.table()
for ( j in seq_along(BXP_REPS) ) {
 bxp.dt <- rbind( bxp.dt,data.table( ID=rownames(get(BXP_REPS[j])),FRAC0=get(BXP_REPS[j])[,"0"],CON=BXP_CONS[j] ) )
}

iii.dt          <- merge(bxp.dt,mirAnnot.pepType.dt[,c("ID","Name","det.STR"),with=F],by="ID",sort=F,all.x=T)
iii.dt[["CON"]] <- factor( x=iii.dt$CON,levels=c("SOM","dYEL","dGRN") )
iii.dt          <- merge( iii.dt,iii.dt[ CON=="SOM",{ list( INDEX=median(FRAC0) ) },by="ID" ][ order(INDEX) ],by="ID",all.x=T,sort=F )
iii.dt          <- iii.dt[ order(INDEX) ]

bxp.5p <-
 ggboxplot(
  data         = iii.dt[ det.STR=="5p" & ( CON %in% c("SOM","dYEL","dGRN") ) ],
  x            = "Name",
  y            = "FRAC0",
  combine      = T,
  fill         = "Name",
  ylab         = "DROSHA Fraction at CP",
  xlab         = "",
  title        = "5p-miRNAs",
  x.text.angle = 90,
  legend       = "none",
  palette      = hue_pal()(51)
 )
my_comparisons <- list( c("SOM","dYEL"),c("SOM","dGRN") )
bxp.5p <-
 bxp.5p +
 geom_point( data=iii.dt[ det.STR=="5p" & ( CON %in% c("SOM","dYEL","dGRN") ) ],mapping=aes(x=Name,y=FRAC0,fill=Name),size=5,shape=21,position=position_jitterdodge() ) +
 facet_wrap( facets=~CON,scales="free",ncol=1 ) +
 coord_cartesian( ylim=c(-0.2,1.2) ) +
 stat_compare_means( comparisons=my_comparisons,method="wilcox.test",paired=F ) +
 theme(
  panel.grid.major.y = element_line( color="grey80" ),
  plot.title   = element_text( size=30,hjust=0.5,face="bold" ), # Adjust main title size
  axis.title.x = element_text( size=20 ),                       # Adjust x-axis label size
  axis.title.y = element_text( size=20,face="bold" ),           # Adjust y-axis label size
  strip.text   = element_text( size=20,face="bold" ),
  axis.text.x  = element_text( face="bold" )
 )

bxp.3p <-
 ggboxplot(
  data         = iii.dt[ det.STR=="3p" & ( CON %in% c("SOM","dYEL","dGRN") ) ],
  x            = "Name",
  y            = "FRAC0",
  combine      = T,
  fill         = "Name",
  ylab         = "DICER1 Fraction at CP",
  xlab         = "",
  title        = "3p-miRNAs",
  x.text.angle = 90,
  legend       = "none",
  palette      = hue_pal()(51)
 )
my_comparisons <- list( c("SOM","dYEL"),c("SOM","dGRN") )
bxp.3p <-
 bxp.3p +
 geom_point( data=iii.dt[ det.STR=="3p" & ( CON %in% c("SOM","dYEL","dGRN") ) ],mapping=aes(x=Name,y=FRAC0,fill=Name),size=5,shape=21,position=position_jitterdodge() ) +
 facet_wrap( facets=~CON,scales="free",ncol=1 ) +
 coord_cartesian( ylim=c(-0.2,1.2) ) +
 stat_compare_means( comparisons=my_comparisons,method="wilcox.test",paired=F ) +
 theme(
  panel.grid.major.y = element_line( color="grey80" ),
  plot.title   = element_text( size=30,hjust=0.5,face="bold" ), # Adjust main title size
  axis.title.x = element_text( size=20 ),                       # Adjust x-axis label size
  axis.title.y = element_text( size=20,face="bold" ),           # Adjust y-axis label size
  strip.text   = element_text( size=20,face="bold" ),
  axis.text.x  = element_text( face="bold" )
 )

dir.create("FDBOXPLOT",showWarnings=F)

pdf( file=paste0("FDBOXPLOT/FDboxplots-",format(Sys.time(), "%Y%m%d"),".pdf"),width=20,height=15 )
 par( bg="white" )

 print(bxp.5p)
 print(bxp.3p)

dev.off()

