library(data.table)
library(matrixStats)
library(pheatmap)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)

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
save( list=c("IDR.CPs.dt",out.objs),file=paste("FD04.IDRs-CPs.normCent.",ONAME,".rda",sep="") )

#= III. unify KO/WT differences as a median-value for each genomic position
cat( "\n ** III. unify KO/WT differences as a median-value for each genomic position **\n\n",sep="" )

 dYEL.diff <-
  diffMat(
   As=paste(c("dYEL.r1","dYEL.r2","dYEL.r3"),".CPs.normCent.mat",sep=""),
   Bs=paste(c("SOM.r1", "SOM.r2", "SOM.r3"), ".CPs.normCent.mat",sep=""),
   FCE="median"
  )
# >> order from -1 => +1 at "0" position ( == increasing )
 dYEL.diff   <- dYEL.diff[ order(     (dYEL.diff[,"0"]),   decreasing=F) , ]
 tmp.dt      <- data.table(            dYEL.diff,keep.rownames="ID" )
 dYEL.CPs.dt <- merge( IDR.CPs.dt,tmp.dt,by="ID",all=T )
# >> order from biggest changes (abs.values) to smallest at "0" position ( == decreasing )
 dYEL.CPs.dt <- dYEL.CPs.dt[ order(abs(dYEL.CPs.dt[["0"]]),decreasing=T) ]

 dGRN.diff <-
  diffMat(
   As=paste(c("dGRN.r1","dGRN.r2","dGRN.r3"),".CPs.normCent.mat",sep=""),
   Bs=paste(c("SOM.r1", "SOM.r2", "SOM.r3"), ".CPs.normCent.mat",sep=""),
   FCE="median"
  )
# >> order from -1 => +1 at "0" position ( == increasing )
 dGRN.diff   <- dGRN.diff[   order(   (dGRN.diff[,"0"]),   decreasing=F) , ]
 tmp.dt     <- data.table(             dGRN.diff,keep.rownames="ID" )
 dGRN.CPs.dt <- merge( IDR.CPs.dt,tmp.dt,by="ID",all=T )
# >> order from biggest changes (abs.values) to smallest at "0" position ( == decreasing )
 dGRN.CPs.dt <- dGRN.CPs.dt[ order(abs(dGRN.CPs.dt[["0"]]),decreasing=T) ]

#= IV. generate heatmaps
cat( "\n ** IV. generate heatmaps **\n\n",sep="" )

dir.create("FDHEATMAP/")
dir.create("FDHEATMAP/p5/")
dir.create("FDHEATMAP/p3/")
dir.create("FDHEATMAP/all/")
rLAB    <- mirAnnot.pepType.dt[ , paste0(" ",Name," ") ]; names(rLAB) <- mirAnnot.pepType.dt[,ID]

idr3.5p <- IDR.mirnas[ IDR.mirnas %in% mirAnnot.pepType.dt[ det.STR=="5p" , ID ] ]
idr3.3p <- IDR.mirnas[ IDR.mirnas %in% mirAnnot.pepType.dt[ det.STR=="3p" , ID ] ]
for ( deCP in c( "dGRN","dYEL" ) ) {

 cat( " >> ",deCP,"\n",sep="" )

 if ( length(idr3.5p[ !( idr3.5p %in% rownames(get(paste0(deCP,".diff")))[ rownames(get(paste0(deCP,".diff"))) %in% idr3.5p ] ) ])==0 ) {
  ph5.mat <- get(paste0(deCP,".diff"))[ idr3.5p, ]
 } else {
  ph5.mat <- get(paste0(deCP,".diff"))[ rownames(get(paste0(deCP,".diff"))) %in% idr3.5p, ]
  nfm.5p  <- idr3.5p[ !( idr3.5p %in% rownames(ph5.mat) ) ]
  ph5.mat <-
   rbind(
    ph5.mat,
    matrix( NA,nrow=length(nfm.5p),ncol=ncol(get(paste0(deCP,".diff"))),dimnames=list( nfm.5p,colnames(ph5.mat) ) )
   )
  ph5.mat <- ph5.mat[ idr3.5p, ]
 }
 if ( length(idr3.3p[ !( idr3.3p %in% rownames(get(paste0(deCP,".diff")))[ rownames(get(paste0(deCP,".diff"))) %in% idr3.3p ] ) ])==0 ) {
  ph3.mat <- get(paste0(deCP,".diff"))[ idr3.3p, ]
 } else {
  ph3.mat <- get(paste0(deCP,".diff"))[ rownames(get(paste0(deCP,".diff"))) %in% idr3.3p, ]
  nfm.3p  <- idr3.3p[ !( idr3.3p %in% rownames(ph3.mat) ) ]
  ph3.mat <-
   rbind(
    ph3.mat,
    matrix( NA,nrow=length(nfm.3p),ncol=ncol(get(paste0(deCP,".diff"))),dimnames=list( nfm.3p,colnames(ph3.mat) ) )
   )
  ph3.mat <- ph3.mat[ idr3.3p, ]
 }

 pdf( file=paste("FDHEATMAP/p5/FDheatmap-5p-",deCP,"-",format(Sys.time(), "%Y%m%d"),".pdf",sep=""),width=20,height=60 )
  par( bg="white",mar=c(1,1,1,1) )

  h5p <-
  pheatmap(
   mat=ph5.mat,cluster_rows=F,cluster_cols=F,
   main=paste0(deCP,"-miRNA-5p","\n\n"),border_color="black",cellwidth=50,cellheight=50,
   color=colorRampPalette(c("blue","white","red"))( length(breakList) ),na_col="gray50",breaks=breakList,
   heatmap_legend_param = list(
     title = "",
     labels_gp = gpar(fontsize = 40),
     legend_height = unit(2, "cm"),
     grid_width = unit(2, "cm"),
     margin = unit(c(0, 1, 0, 4), "cm")
   ),legend=F,
   labels_row=rLAB[idr3.5p],labels_col=rep.int("",ncol(ph5.mat)),fontsize=50,silent=T
  )

  draw( h5p,padding=unit( c( 2,100,2,180 ),"mm" ) )

  my_colors <- colorRamp2(
   c( -1, 0, +1 ),
   c("blue", "white", "red")   # Adjust to match your preferred scale
  )
  lgd <- Legend(
   col_fun = my_colors,
   labels_gp = gpar(fontsize = 40),
   legend_height = unit( 5, "cm"),
   grid_width = unit(2, "cm")
  )

  draw(
   lgd,
   x = unit(1, "cm"),     # Exactly 5cm from the left edge of the PDF
   y = unit(30, "cm"),    # Exactly 30cm from the bottom edge of the PDF
   just = c("left", "bottom") # Anchors the bottom-left corner of the legend box to that exact point
  )

 dev.off()
 pdf( file=paste("FDHEATMAP/p3/FDheatmap-3p-",deCP,"-",format(Sys.time(), "%Y%m%d"),".pdf",sep=""),width=20,height=60 )
  par( bg="white",mar=c(1,1,1,1) )

  h3p <-
  pheatmap(
   mat=ph3.mat,cluster_rows=F,cluster_cols=F,
   main=paste0(deCP,"-miRNA-3p","\n\n"),border_color="black",cellwidth=50,cellheight=50,
   color=colorRampPalette(c("blue","white","red"))( length(breakList) ),na_col="gray50",breaks=breakList,
   heatmap_legend_param = list(
     title = "",
     labels_gp = gpar(fontsize = 40),
     legend_height = unit(2, "cm"),
     grid_width = unit(2, "cm"),
     margin = unit(c(0, 1, 0, 4), "cm")
   ),legend=F,
   labels_row=rLAB[idr3.3p],labels_col=rep.int("",ncol(ph3.mat)),fontsize=50,silent=T
  )

  draw( h3p,padding=unit( c( 2,100,2,180 ),"mm" ) )

  my_colors <- colorRamp2(
   c( -1, 0, +1 ),
   c("blue", "white", "red")   # Adjust to match your preferred scale
  )
  lgd <- Legend(
   col_fun = my_colors,
   labels_gp = gpar(fontsize = 40),
   legend_height = unit( 5, "cm"),
   grid_width = unit(2, "cm")
  )

  draw(
   lgd,
   x = unit(1, "cm"),     # Exactly 5cm from the left edge of the PDF
   y = unit(30, "cm"),    # Exactly 30cm from the bottom edge of the PDF
   just = c("left", "bottom") # Anchors the bottom-left corner of the legend box to that exact point
  )

 dev.off()
 pdf( file=paste("FDHEATMAP/all/FDheatmap-5p3p-",deCP,"-",format(Sys.time(), "%Y%m%d"),".pdf",sep=""),width=20,height=120 )
  par( bg="white",mar=c(1,1,1,1) )

  h53 <- h5p %v% h3p

  draw( h53,ht_gap=unit( 1.5,"cm" ),padding=unit( c( 2,100,2,180 ),"mm" ) )

  my_colors <- colorRamp2(
   c( -1, 0, +1 ),
   c("blue", "white", "red")   # Adjust to match your preferred scale
  )
  lgd <- Legend(
   col_fun = my_colors,
   labels_gp = gpar(fontsize = 40),
   legend_height = unit( 5, "cm"),
   grid_width = unit(2, "cm")
  )

  draw(
   lgd,
   x = unit(1, "cm"),     # Exactly 5cm from the left edge of the PDF
   y = unit(70, "cm"),    # Exactly 30cm from the bottom edge of the PDF
   just = c("left", "bottom") # Anchors the bottom-left corner of the legend box to that exact point
  )

 dev.off()

}

