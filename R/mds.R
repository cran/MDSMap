#' Estimate marker positions using spherically constrained weighted MDS
#'
#' Reads a text file of pairwise recombination fractions and LOD scores, 
#' estimates marker positions using spherically constrained weighted MDS
#'
#' @param fname Character string specifying the base name of the file 
#' fname.txt which contains the data to be analysed. The file should be 
#' white space or tab separated.
#' @param p Integer - the penalty for deviations from the sphere - higher p 
#' forces points more closely onto a sphere.
#' @param n Vector of integers or strings containing markers to be omitted from 
#' the analysis.
#' @param weightfn Character string specifying the values to use for the weight 
#' matrix in the MDS 'lod2' or 'lod'.
#' @param mapfn Character string specifying the map function to use on the 
#' recombination fractions 'haldane' is default, 'kosambi' or 'none'.
#'
#' @details
#' This can be very slow with large sets of markers, in which case it may be 
#' better to consider \code{\link{calc.maps.pc}}.
#'
#' Reads a file of the form described below and casts the data into matrices of 
#' pairwise recombination fractions and weights determined by the \code{weightfn} 
#' parameter (\code{LOD} or \code{LOD^2^}) calculates a distance matrix from the map 
#' function. Haldane is the default map function, None just uses recombination 
#' fractions and the other alternative is Kosambi (see \code{link{dmap}} 
#' for details). 
#'
#' Performs both an unconstrained and dual spherically constrained weighted MDS 
#' on the distance matrix using \code{\link[smacof]{smacofSym}} and 
#' \code{\link[smacof]{smacofSphere}} (\cite{de Leeuw & Mair 2009}) 
#' and maps this to an interval (see \code{\link{map.to.interval}} for details).
#'
#' Inevitably the constrained MDS has higher stress than the unconstrained MDS and 
#' a good rule of thumb is that this should not be more than about 10% higher.
#'
#' File names should be of the form \code{fname.txt} and it is assumed that they are in 
#' a tab or space separated file of the format displayed below. The first entry on 
#' the first row is the number of markers to be analysed. Underneath this is a 
#' table in which the first two columns contain marker names, the third column 
#' contains the pairwise recombination fractions between the markers and the 
#' fourth column the associated LOD score.  Note that marker names in the first 
#' column vary more slowly than in the second column. Missing recombination pairs 
#' are acceptable. Recombination fractions greater than 0.499999 are set to that 
#' value.
#'
#' \tabular{llll}{
#' \cr
#' \code{nmarkers} \tab \tab \tab \cr
#' \code{marker_1} \tab \code{marker_2} \tab \code{recombination fraction} \tab \code{LOD}\cr
#' \code{1} \tab \code{2} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' }
#'
#' @return A list (S3 class 'spheremap') with the following elements:
#' \item{smacofsym}{The unconstrained wMDS results.}
#' \item{smacofsphere}{The spherically constrained wMDS results.}
#' \item{mapsphere}{Map of the markers onto an interval containing order-the 
#' rank of each marker.}
#' \item{distmap}{A symmetric matrix of pairwise distances between markers where 
#' the columns are in the estimated order.}
#' \item{lodmap}{A symmetric matrix of lod scores associated with the distances in distmap.}
#' \item{locimap}{A data frame of the markers containing the name of each marker, 
#' the number in the configuration plot if that is being used, the position of 
#' each marker in order of increasing distance and the nearest neighbour fit of 
#' the marker.}
#' \item{length}{Integer giving the total length of the segment.}
#' \item{removed}{A vector of the names of markers removed from the analysis.}
#' \item{locikey}{A data frame showing the number associated with each marker 
#' name for interpreting the wMDS configuration plots.}
#' \item{stressratio}{The ratio of the constrained to unconstrained stress.}
#' \item{ssphere}{The stress per point of the spherically constrained wMDS.}
#' \item{ssym}{Stress per point of the unconstrained wMDS.}
#' \item{meannnfit}{The mean across all markers of the nearest neighbour fits.}
#'
#' @references
#' \cite{de Leeuw J, Mair P (2009) Multidimensional scaling using majorization: SMACOF 
#' in R. J Stat Softw 31: 1-30} \url{http://www.jstatsoft.org/v31/i03/}
#' @seealso
#' \code{\link{calc.maps.pc}}, \code{\link{calc.pair.rf.lod}}, \code{\link[smacof]{smacofSym}}, \code{\link[smacof]{smacofSphere}}, \code{\link{map.to.interval}}, \code{\link{dmap}}, \code{\link{calc.nnfit}}
#'
#' @examples
#' smap<-calc.maps.sphere(system.file("extdata", "lgI.txt", package="MDSMap"),
#' weightfn='lod',mapfn='kosambi')
#' @export
calc.maps.sphere<-function(fname,p=100,n=NULL,weightfn='lod2',mapfn='haldane'){
  lodrf<-calc.pair.rf.lod(fname,weightfn)
  confplotno<-1:lodrf$nloci
  if(!is.null(n)){
    if(!is.numeric(n))n<-which(lodrf$locinames%in%n)    
    r<-lodrf$rf[-n,-n]
    lod<-lodrf$lod[-n,-n]
    confplotno<-confplotno[-n]
  } else {
    r<-lodrf$rf
    lod<-lodrf$lod
  }
  M<-dmap(r,mapfn)
  nloci=length(confplotno)
  smacofsym<-smacof::smacofSym(M,ndim=2,weightmat=lod,itmax=100000)
  smacofsphere<-smacof::smacofSphere(M,ndim=2,algorithm="dual",weightmat=lod,penalty=p,itmax=1000000,mod=10,verbose=FALSE)
  mapsphere<-map.to.interval(smacofsphere,nloci)
  length<-mapsphere$chromlength[nloci]
  distmap<-outer(mapsphere$maporder,mapsphere$maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(mapsphere$maporder,mapsphere$maporder,Vectorize(function(i,j)lod[i,j]))
  #stressratio=smacofsphere$stress/smacofsym$stress
  if(!is.null(n)) {
	locikey<-data.frame(locus=lodrf$locinames[-n],confplotno=confplotno)
  } else {
	locikey<-data.frame(locus=lodrf$locinames,confplotno=confplotno)
  }
  sr=smacofsphere$stress/smacofsym$stress
  ssphere=smacofsphere$stress
  ssym=smacofsym$stress
  nnfit<-calc.nnfit(distmap,lodmap,mapsphere$chromlength)
  locimap<-data.frame(confplotno=confplotno[mapsphere$maporder],
	locus=locikey$locus[mapsphere$maporder],position=mapsphere$chromlength,
	nnfit=nnfit$pointfits,row.names=1:nloci)

  if(!is.null(n)) {
    removedloci<-data.frame(n,lodrf$locinames[n],row.names=NULL) 
  } else {
    removedloci<-n
  }
  retlist<-list(smacofsym=smacofsym,smacofsphere=smacofsphere,mapsphere=mapsphere,distmap=distmap,
	lodmap=lodmap,locimap=locimap,length=length,removed=n,locikey=locikey,stressratio=sr,
	ssphere=ssphere,ssym=ssym,meannnfit=nnfit$meanfit)
  class(retlist) <- "spheremap"
  retlist
}


#' Estimate marker positions using Principal Curves
#' 
#' Reads a text file of pairwise recombination fractions and LOD scores, reduces 
#' to 2 or 3 dimensions using wMDS and projects onto a single dimension using 
#' principal curves to estimate marker positions.
#'
#' @param fname Character string the name of the file of recombination fractions 
#' and scores it should not contain any suffices (the file should be a .txt file 
#' as described below).
#' @param spar Integer - the smoothing parameter for the principal curve. If 
#' NULL this will be done using leave one out cross validation.
#' @param n	Vector of integers or character strings containing markers to be 
#' omitted from the analysis.
#' @param ndim Number of dimensions in which to perform the wMDS and fit the 
#' curve - can be 2 or 3.
#' @param weightfn Character string specifying the values to use for the weight 
#' matrix in the MDS \code{'lod2'} or \code{'lod'}.
#' @param mapfn Character string specifying the map function to use on the 
#' recombination fractions \code{'haldane'} is default, \code{'kosambi'} or \code{'none'}.
#'
#' @details
#' Reads a file of the form described below and casts the data into matrices of 
#' pairwise recombination fractions and weights determined by the \code{weightfn} 
#' parameter (\code{LOD} or \code{LOD^2^}) calculates a distance matrix from the map function. 
#' Haldane is the default map function, none just uses recombination fractions 
#' and the other alternative is Kosambi (see \code{\link{dmap}} for details). 
#'
#' Performs both an weighted MDS on the distance matrix using \code{\link[smacof]{smacofSym}} and 
#' \code{\link[smacof]{smacofSphere}} (\cite{de Leeuw & Mair 2009}) and fits a 
#' principal curve to map this to an interval (\code{\link[princurve]{principal.curve}} for details).
#'
#' File names should be of the form \code{fname.txt} and it is assumed that they are in 
#' a tab or space separated file of the format displayed below. The first entry on 
#' the first row is the number of markers to be analysed. Underneath this is a 
#' table in which the first two columns contain marker names, the third column 
#' contains the pairwise recombination fractions between the markers and the 
#' fourth column the associated lod score.  Note that marker names in the first 
#' column vary more slowly than in the second column. Missing recombination pairs 
#' are acceptable. Recombination fractions greater than 0.499999 are set to that 
#' value.
#'
#' \tabular{llll}{
#' \cr
#' \code{nmarkers} \tab \tab \tab \cr
#' \code{marker_1} \tab \code{marker_2} \tab \code{recombination fraction} \tab \code{LOD}\cr
#' \code{1} \tab \code{2} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' }
#'
#' @return A list (S3 class pcmap or pcmap3d depending on ndim) with the following elements:
#' \item{smacofsym}{The unconstrained wMDS results.}
#' \item{pc}{The results from the principal curve fit.}
#' \item{distmap}{A symmetric matrix of pairwise distances between markers where 
#' the columns are in the estimated order.}
#' \item{lodmap}{A symmetric matrix of lod scores associated with the distances 
#' in distmap.}
#' \item{locimap}{A data frame of the markers containing the name of each marker, 
#' the number in the configuration plot if that is being used, the position of each 
#' marker in order of increasing distance and the nearest neighbour fit of the marker.}
#' \item{length}{Integer giving the total length of the segment.}
#' \item{removed}{A vector of the names of markers removed from the analysis.}
#' \item{locikey}{A data frame showing the number associated with each marker name 
#' for interpreting the wMDS configuration plots.}
#' \item{meannnfit}{The mean across all markers of the nearest neighbour fits.}
#'
#' @references
#' \cite{de Leeuw J, Mair P (2009) Multidimensional scaling using majorization: SMACOF in R. J Stat Softw 31: 1-30} \url{http://www.jstatsoft.org/v31/i03/}
#'
#' \cite{Hastie T, Weingessel A (2013) princurve: Fits a Principal Curve in Arbitrary Dimension. ) R package version 1.1-12.} \url{https://CRAN.R-project.org/package=princurve}
#'
#' @seealso
#' \code{\link{calc.maps.sphere}}, \code{\link{calc.pair.rf.lod}}, \code{\link[smacof]{smacofSym}}, \code{\link[smacof]{smacofSphere}}, \code{\link{map.to.interval}}, \code{\link{dmap}}
#'
#' @examples
#' map<-calc.maps.pc(system.file("extdata", "lgV.txt", package="MDSMap"),
#' ndim=2,weightfn='lod2',mapfn='kosambi')
#' plot(map)
#' @export

calc.maps.pc<-function(fname,spar=NULL,n=NULL,ndim=2,weightfn='lod2',mapfn='haldane'){
  lodrf<-calc.pair.rf.lod(fname,weightfn)
  confplotno<-1:lodrf$nloci
  if(!is.null(n)){
    if(!is.numeric(n))n<-which(lodrf$locinames%in%n)    
    r<-lodrf$rf[-n,-n]
    lod<-lodrf$lod[-n,-n]
    confplotno<-confplotno[-n]
  } else {
    r<-lodrf$rf
    lod<-lodrf$lod
  }
  M<-dmap(r,mapfn)
  nloci=length(confplotno)

  smacofsym<-smacof::smacofSym(M,ndim=ndim,weightmat=lod,itmax=100000)
  pc1<-princurve::principal.curve(smacofsym$conf,maxit=150,spar=spar)
  scale<-sum(smacofsym$delta)/sum(smacofsym$dhat) 
  # Configuration dissim are based on the normalized observed diss - dhat. 
  # True observed dissimilarities are delta
  maporder<-pc1$tag
  estpos<-pc1$lambda[maporder]*scale*100
  # gives the estimated length from the beginning of the line
  rownames<-lodrf$locinames[maporder]
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder, Vectorize(function(i,j)lod[i,j]))
  rownames(distmap)<-rownames;colnames(distmap)<-rownames
  rownames(lodmap)<-rownames;colnames(lodmap)<-rownames
  if(!is.null(n))  {
    locikey<-data.frame(locus=lodrf$locinames[-n],confplotno=confplotno)
  } else {
    locikey<-data.frame(locus=lodrf$locinames,confplotno=confplotno)
  }
  nnfit<-calc.nnfit(distmap,lodmap,estpos)
  locimap<-data.frame(confplotno=confplotno[maporder],locus=locikey$locus[maporder],position=estpos,nnfit=nnfit$pointfits,row.names=1:nloci)
  if(!is.null(n)) {
    removedloci<-data.frame(n,lodrf$locinames[n],row.names=NULL)
  } else {
    removedloci<-n
  }
  
  retlist<-list(smacofsym=smacofsym,pc=pc1,distmap=distmap,lodmap=lodmap,locimap=locimap,length=max(estpos),removed=n,locikey=locikey,meannnfit=nnfit$meanfit)
  if(ndim == 2) {
	class(retlist) <- "pcmap"
  } else {
    class(retlist) <- "pcmap3d"
  }
  retlist
}


#' Create recombination matrix from pairwise data file.
#' 
#' Reads a text file of pairwise recombination fractions and LOD scores and casts 
#' it into a matrix of recombination fractions and weights.
#'
#' @param fname Character string specifying the base name of the file \code{fname.txt} which 
#' contains the data to be analysed the file should be white space or tab separated.
#' @param weightfn Character string specifiying the values to use for the weight 
#' matrix \code{'lod2'} or \code{'lod'}.
#' @param ... \code{\link[utils]{read.table}} arguments.
#'
#' @details
#' File names should be of the form \code{fname.txt} and it is assumed that they are in 
#' a tab or space separated file of the format displayed below. The first entry on 
#' the first row is the number of markers to be analysed. Underneath this is a 
#' table in which the first two columns contain marker names, the third column 
#' contains the pairwise recombination fractions between the markers and the 
#' fourth column the associated LOD score.  Note that marker names in the first 
#' column vary more slowly than in the second column. Missing recombination pairs 
#' are acceptable. Recombination fractions greater than 0.499999 are set to that 
#' value
#'
#' \tabular{llll}{
#' \cr
#' \code{nmarkers} \tab \tab \tab \cr
#' \code{marker_1} \tab \code{marker_2} \tab \code{recombination fraction} \tab \code{LOD}\cr
#' \code{1} \tab \code{2} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' }
#'
#' @return A list with the following elements:
#' \item{rf}{A symmetric matrix of recombination fractions.}
#' \item{nloci}{The number of markers in the analysis.}
#' \item{locinames}{The names of the markers in the analysis.}
#'
#' @examples
#' lodrf<-calc.pair.rf.lod(system.file("extdata", "lgV.txt", package="MDSMap"), 
#' "lod2")
#' @export
calc.pair.rf.lod<-function(fname,weightfn='lod',...){
  if(!file.exists(fname)) {
    fname2 <- paste(fname,'.txt',sep="")
	if(file.exists(fname2)) {
		fname <- fname2
	}
  }  
  nloci<-scan(fname,what=integer(),nmax=1)
  d<-utils::read.table(fname,skip=1,header=FALSE)
  names(d)<-c("name1","name2","rfreq","lodscore")
  if(weightfn=='lod2') d$lodscore<-d$lodscore^2
  dd<-d
  
  missing1<-with(d,unique(as.character(name1[!name1%in%name2])))
  missing2<-with(d,as.character(unique(name2[!name2%in%name1])))
  
  if(length(missing1)>1){
    dd$name1<-as.character(dd$name1);dd$name2<-as.character(dd$name2)
    for(i in 2:length(missing1))dd<-rbind(dd,list(missing1[1],missing1[i],0,0))
  }
  if(length(missing2)>1){
    dd$name1<-as.character(dd$name1);dd$name2<-as.character(dd$name2)
    for(i in 2:length(missing2))dd<-rbind(dd,list(missing2[i],missing2[1],0,0))
  }
  dd$name1<-factor(dd$name1,unique(as.character(dd$name1)))
  dd$name1<-stats::relevel(dd$name1,missing1[1])
  dd$name2<-factor(dd$name2,levels=c(as.character(levels(dd$name1)[2:length(levels(dd$name1))]),as.character(missing2[1])))

  d<-dd
  b<-matrix(0,ncol=nloci,nrow=nloci)
  temp<-reshape::cast(name1~name2,data=d,value="rfreq",add.missing=TRUE,fill=0)
  tt<-as.matrix(temp[,2:(nloci)])
  colnames(tt)<-names(temp)[2:nloci]
  rownames(tt)<-temp$name1
  b[upper.tri(b)]<-tt[upper.tri(tt,diag=TRUE)]
  rfmat<-b+t(b)
  colnames(rfmat)<-c(rownames(tt)[1],colnames(tt))
  rownames(rfmat)<-c(rownames(tt),colnames(tt)[nloci-1])
  rm(temp,tt,b)
  b<-matrix(0,ncol=nloci,nrow=nloci)  
  temp<-reshape::cast(name1~name2,data=d,value="lodscore",add.missing=TRUE,fill=0)
  tt<-as.matrix(temp[,2:(nloci)])  
  colnames(tt)<-names(temp)[2:nloci]
  rownames(tt)<-temp$name1
  b[upper.tri(b)]<-tt[upper.tri(tt,diag=TRUE)]
  lodmat<-b+t(b)  
  colnames(lodmat)<-c(rownames(tt)[1],colnames(tt))
  rownames(lodmat)<-c(rownames(tt),colnames(tt)[nloci-1])
  rfmat[rfmat>0.499999]<-0.499999
  rm(temp,tt,b)
  diag(rfmat)<-NA
  diag(lodmat)<-NA
  lodmat[lodmat<0]<-0  
  list(rf=rfmat,lod=lodmat,nloci=nloci,locinames=rownames(rfmat))
}

#' Load data, estimate a linkage map and plot diagnostics for the fit.
#'
#' Load data, estimate a linkage map and plot diagnostics for the fit.
#'
#' @param fname Character string containing the base file from which the data 
#' should be read - should contain the complete file name excluding the suffix 
#' which should be \code{.txt}
#' @param p Smoothing parameter.
#' @param n Vector of integers or character strings containing the name or 
#' position in the input list of loci to be excluded from the analysis.
#' @param ispc Logical determining the method to be used to estimate the map. By 
#' default this is \code{TRUE} and the method of principal curves will be used. If 
#' \code{FALSE} then the constrained MDS method will be used.
#' @param ndim Integer the number of dimensions to use if the Principal curves 
#' method is used. By default this is 2, but it can also be 3.
#' @param weightfn Character string specifying the values to use for the weight 
#' matrix in the MDS \code{lod2} or \code{lod}.
#' @param mapfn Character string specifying the map function to use on the 
#' recombination fractions \code{'haldane'} is default, \code{'kosambi'} or \code{'none'}.
#' @param D1lim Numeric vector specifying the limits of the axis relating to 
#' dimension 1 of the wMDS used to estimate the map.
#' @param D2lim Numeric vector specifying the limits of the axis relating to 
#' dimension 1 of the wMDS used to estimate the map.
#' @param D3lim Numeric vector specifying the limits of the axis relating to 
#' dimension 1 of the wMDS used to estimate the map.
#' @param displaytext Logical argument determining how markers should be labelled 
#' in the wMDS configuration plot. If \code{TRUE} then marker names are used. If 
#' \code{FALSE} then numbers are used.
#'
#' @details 
#' Data is read from a text file which should be of the form described below. 
#' By default, \code{ispc=TRUE}, in which case maps are estimated using unconstrained 
#' weighted MDS followed by fitting a principal curve. Details can be found in 
#' the description of the function \code{\link{calc.maps.pc}}. If \code{ispc=FALSE} 
#' maps are estimated using spherically constrained weighted MDS. Details can be 
#' found in the description of the function \code{\link{calc.maps.sphere}}.
#'
#' \code{ndim} is only relevant if \code{ispc=TRUE}, in which case it specifies the number of 
#' dimensions to be used, the default is 2 but it can also be 3 dimensions. 
#'
#' Diagnostic plots are then produced using \code{\link{plot.pcmap}} for the method of 
#' principal curves in 2 dimensions, \code{\link{plot.pcmap3d}} for the method of principal 
#' curves in 3 dimensions and \code{\link{plot.spheremap}} for the method using spherically 
#' constrained MDS.
#'
#' \code{n} specifies markers to be omitted from the analysis. It can be a vector of 
#' character strings specifying makers to be omitted, or a vector of integers 
#' specifying the markers to omit. The latter method is likely to be useful when 
#' removing outliers after inspection of the diagnostic plot, because the output 
#' contains a dataframe, locikey, which associates each marker with its 
#' identifying number. By default this is NULL and all markers in the file will 
#' be analysed.
#'
#' \code{p} is a smoothing parameter which operates quite differently depending on 
#' whether map estimation is performed using Principal Curves or Constrained 
#' MDS. If the PC method is used, \code{p} determines the smoothing parameter spar in 
#' the function \code{\link[princurve]{principal.curve}} from the package 
#' \pkg{princurve}. If \code{NULL} then the most appropriate value will be determined 
#' using leave one out cross validation. 
#' If Constrained MDS is used then \code{p} must be set to a number which specifies the 
#' penalty for deviations from the sphere in the function \code{\link[smacof]{smacofSphere}} from the 
#' \pkg{smacof} package. Something between 50 and 100 is generally appropriate and this 
#' penalty can be decreased if stress from the constrained analysis is more than 
#' about 10% higher than that from the unconstrained MDS (see \code{\link{calc.maps.sphere}} 
#' for details)
#'
#' File names should be of the form \code{fname.txt} and it is assumed that they are in 
#' a tab or space separated file of the format displayed below. The first entry on 
#' the first row is the number of markers to be analysed. Underneath this is a 
#' table in which the first two columns contain marker names, the third column 
#' contains the pairwise recombination fractions between the markers and the 
#' fourth column the associated LOD score.  Note that marker names in the first 
#' column vary more slowly than in the second column. Missing recombination pairs 
#' are acceptable. Recombination fractions greater than 0.499999 are set to that 
#' value.
#'
#' \tabular{llll}{
#' \cr
#' \code{nmarkers} \tab \tab \tab \cr
#' \code{marker_1} \tab \code{marker_2} \tab \code{recombination fraction} \tab \code{LOD}\cr
#' \code{1} \tab \code{2} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' }
#'
#' @return map (s3 class pcmap, pcmap3d or spheremap) from \code{\link{calc.maps.pc}} if \code{ispc=TRUE} or 
#' \code{\link{calc.maps.sphere}} if \code{ispc=FALSE}.
#'
#' @references
#' \cite{de Leeuw J, Mair P (2009) Multidimensional scaling using majorization: SMACOF in R. J Stat Softw 31: 1-30} \url{http://www.jstatsoft.org/v31/i03/}
#'
#' \cite{Hastie T, Weingessel A (2013) princurve: Fits a Principal Curve in Arbitrary Dimension. ) R package version 1.1-12.} \url{https://CRAN.R-project.org/package=princurve}
#'
#' @seealso \code{\link[smacof]{smacofSphere}}, \code{\link[princurve]{principal.curve}}, \code{\link{calc.maps.pc}}, \code{\link{calc.maps.sphere}}, \code{\link{plot.pcmap}}, \code{\link{plot.pcmap3d}}, \code{\link{plot.spheremap}}
#'
#' @examples
#' estimate.map(system.file("extdata", "lgI.txt", package="MDSMap"),
#' ndim=3)
#' @export
estimate.map<-function(fname,p=NULL,n=NULL,ispc=TRUE,ndim=2,weightfn='lod2',mapfn='haldane',D1lim=NULL,D2lim=NULL,D3lim=NULL,displaytext=TRUE){
  if(!is.null(n)){
    t<-table(n)[(table(n)>1)]
    if(length(t)>0){
      write("Warning: the list of ommitted markers is not unique, please remove duplicates of the loci listed below","")
      write(names(t),"")
      return(0)
    }
  }
  if(ispc==FALSE){
    map<-calc.maps.sphere(fname,p,n,weightfn=weightfn,mapfn=mapfn)
    graphics::plot(map,displaytext=displaytext)
  } else {
    map<-calc.maps.pc(fname,spar=p,n,ndim=ndim,weightfn=weightfn,mapfn=mapfn)
    if(ndim==2) {
	  graphics::plot(map,displaytext=displaytext)
	} else { 
	  graphics::plot(map,D1lim,D2lim,D3lim,displaytext=displaytext)
	}
  }
  write(paste('Stress:',map$smacofsym$stress),"")
  write(paste('Mean Nearest Neighbour Fit:',map$meannnfit),"")
  write('Markers omitted:',"")
  write(n,"")
  map
}

#' Nearest neighbour fit from estimated map and file of pairwise recombination fractions.
#' 
#' Calculates a nearest neighbour fit based from an estimated map and a file 
#' containing pairwise recombination fractions and LOD scores.
#'
#' @param estmap A character string indicating the name of a comma separated 
#' value file with  the first column containing marker names in the order 
#' of their estimated position.
#' @param fname A character string specifying the base name of the file 
#' \code{fname.txt} which contains the data to be analysed the file should be 
#' white space or tab separated.
#' @param mapfn Character string, \code{'haldane'}, \code{'kosambi'} or \code{'none'}
#' specifying the values to use to estimate the map distance from the 
#' recombination fractions. Default is \code{'haldane'}.
#' @param n Vector of character strings or numbers specifying the markers to be 
#' omitted from the analysis. Default is \code{NULL}.
#' @param header Logical argument indicating whether the .csv file estmap contains headers - default is TRUE
#'
#' @details
#' Reads in two files \code{fname.txt} and \code{estmap}. 
#'
#' The data is cast the data into symmetric matrices of pairwise recombination 
#' fractions and LOD scores with the order of columns and rows in the matrix 
#' determined by the order specified in \code{estmap}. A distance matrix is 
#' calculated according to the method specified by \code{mapfn}. Haldane is the 
#' default map function, None just uses recombination fractions and the other 
#' alternative is Kosambi (see \code{\link{dmap}} for details). The nearest 
#' neighbour fit is then calculated (see \code{\link{calc.nnfit}} for details)
#'
#' \code{estmap} should contain marker names in the first column in the order 
#' of the estimated map.
#'
#' \code{fname} should be of the form \code{fname.txt} and it is assumed that they are in 
#' a tab or space separated file of the format displayed below. The first entry on 
#' the first row is the number of markers to be analysed. Underneath this is a 
#' table in which the first two columns contain marker names, the third column 
#' contains the pairwise recombination fractions between the markers and the 
#' fourth column the associated LOD score.  Note that marker names in the first 
#' column vary more slowly than in the second column. Missing recombination pairs 
#' are acceptable. Recombination fractions greater than 0.499999 are set to that 
#' value.
#'
#' \tabular{llll}{
#' \cr
#' \code{nmarkers} \tab \tab \tab \cr
#' \code{marker_1} \tab \code{marker_2} \tab \code{recombination fraction} \tab \code{LOD}\cr
#' \code{1} \tab \code{2} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' }
#'
#' @return A list with the following elements:
#' \item{fit}{Sum over all markers of the nearest neighbour fits.}
#' \item{pointfits}{The nearest neighbour fit for each marker.}
#' \item{meanfit}{Mean of the nearest neighbour fits over all markers.}
#' @seealso
#' \code{\link{dmap}}, \code{\link{calc.nnfit}}, \code{\link{calc.pair.rf.lod}}
#' @export
calc.nnfit.from.file<-function(estmap,fname,mapfn='haldane',n=NULL,header=FALSE){
  estmap<-utils::read.csv(estmap,header=header)
  lodrf<-calc.pair.rf.lod(fname)
  if(!is.null(n)){
    if(!is.numeric(n))n<-which(lodrf$locinames%in%n)    
    r<-lodrf$rf[-n,-n]
    lod<-lodrf$lod[-n,-n]
  } else {
    r<-lodrf$rf
    lod<-lodrf$lod
  }
  M<-dmap(r,mapfn)
  lnames<-colnames(M)
  names<-estmap[,1]
  maporder<- sapply(1:length(names),function(i)which(lnames==names[i]))
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder,Vectorize(function(i,j)lod[i,j]))
  nnfit<-calc.nnfit(distmap,lodmap,estmap[,2])
  newmap<-data.frame(name=estmap[,1],position=estmap[,2],nnfit=nnfit$pointfits)
  if(!is.null(fname)) utils::write.table(newmap,file=fname,sep=',')
  nnfit
}


#' Calculates pairwise map distances from the recombination fraction.
#' 
#' Calculates pairwise map distances from the recombination fraction.
#'
#' @param rf A symmetric matrix of pairwise recombination fractions.
#' @param mapfn A character string specifying the map function to be 
#' used in calculated the distance \code{'haldane'}, \code{'kosambi'}, \code{'none'}.
#' 
#' @details
#' The default is the \code{'haldane'} map function \eqn{0.5\ln(1-2rf)}{0.5ln(1-2rf)}, \code{'kosambi'} returns 
#' \eqn{0.25\ln((1+2rf)/(1-2rf))}{0.25ln((1+2rf)/(1-2rf))} and \code{'none'} returns \eqn{rf}, the recombination fraction.
#'
#' @return a symmetric matrix of pairwise map distances in the same format as the recombination matrix supplied.
#'
#' @examples
#' lodrf<-calc.pair.rf.lod(system.file("extdata", "lgV.txt", package="MDSMap"))
#' mdist=dmap(lodrf$rf,mapfn="haldane")
#' @export
dmap<-function(rf,mapfn="haldane"){
  if (mapfn=="haldane") return(-0.5*log(1-2*rf)) 
  if (mapfn=="kosambi") return(0.25*log((1+2*rf)/(1-2*rf)))
  if (mapfn=="none") return (rf)
}

#' Convert Cartesian coordinates from wMDS coordinates to polar coordinates.
#' 
#' Converts the coordinates of points in the final configuration of a spherically 
#' constrained wMDS from Cartesian to polar coordinates. 
#'
#' @param mdsobject Output from \code{\link[smacof]{smacofSphere}} using the dual method in the \pkg{smacof} package.
#' @param nloci The number of markers in the configuration.
#'
#' @details
#' Centres the circle on zero if necessary, finds a the most natural break in 
#' the points to start as 0, then calculates the angle of each point relative 
#' to this. The radius is the median distance of points from the centre.
#'
#' @return
#' \item{theta}{A vector of angles one for each point.}
#' \item{radius}{A scalar the radius of sphere.}
#'
#' @references 
#' \cite{de Leeuw J, Mair P (2009) Multidimensional scaling using majorization: SMACOF in R. J Stat Softw 31:1-30} \url{http://www.jstatsoft.org/v31/i03/}
#' @seealso
#' \code{\link[smacof]{smacofSphere}}
#' @examples
#' #M and lod should be n x n symmetric matrices of the same dimensions where n 
#' #is the number markers to be analysed
#' \dontrun{
#' mds1<-smacofSphere(M,ndim=2,algorithm="dual",weightmat=lod,penalty=100)
#' pol<-convert.polar(mds1,n)
#' }
#' @export
convert.polar<-function(mdsobject,nloci){
  conf=mdsobject$conf
  l<-dim(conf)[1]
  start<-l+1-nloci
  if(start>1){
    x<-conf[start:l,1]-conf[1,1]
    y<-conf[start:l,2]-conf[1,2]
  } else {
    x<-conf[start:l,1]
    y<-conf[start:l,2]
  }
  yadd<-ifelse(y<=0,2*pi,0)
  xadd<-ifelse(x<=0,pi,yadd)
  theta<-atan(y/x)+xadd
  newtheta<-sort(theta)
  diff=newtheta[2:(length(newtheta))]-newtheta[1:(length(newtheta)-1)]
  maxd<-max(diff)
  rotation<-ifelse(maxd>(pi/3),-min(newtheta[(which(diff==maxd)+1):length(theta)]),0)
  
  rtheta<-(theta+rotation)%%(2*pi)
  radius<-sqrt(x^2+y^2)
  list(theta=(rtheta-min(rtheta))%%(2*pi),radius=radius)  
}


#' Map points from MDS final configuration to interval starting at 0.
#' 
#' Maps points from the final configuration of a 2-dimensional spherically 
#' constrained wMDS to an interval with a starting point at 0.
#'
#' @param mdsobject The output from \code{\link[smacof]{smacofSphere}}.
#' @param nloci The number of markers in the configuration.
#'
#' @details
#' Centres the configuration on zero and calculates the median distance of the points from the origin.
#' Finds the largest gap in the spherical configuration and assigns the marker on the right hand side of it angle 0. 
#' Converts Cartesian coordinates to polar coordinates and projects points onto 
#' the arc centred on 0 with radius the median distance from the origin. 
#'
#' @return A list with the elements:
#' \item{chromlength}{A named vector giving the position of each marker.}
#' \item{order}{A named vector giving the rank order of the markers.}
#' \item{locilength}{A named vector giving the position of each marker in order 
#' of increasing distance along the segment.}
#' \item{maporder}{A named vector of the position in the input list of each 
#' marker in order of increasing distance along the segment.}
#'
#' @seealso
#' \code{\link{convert.polar}}
#' @examples
#' # M and lod should be n x n symmetric matrices of the same dimensions where
#' # n is the number markers to be analysed
#' \dontrun{
#' mds1<-smacofSphere(M,ndim=2,algorithm=dual,weightmat=lod,penalty=100)
#' pol<-map.to.interval (m1,n)
#' }
#' @export
map.to.interval<-function(mdsobject,nloci){
  pol<-convert.polar(mdsobject,nloci) #detrend
  lin<-pol$theta
  radmed<-stats::median(pol$radius)
  scale<-sum(mdsobject$delta[lower.tri(mdsobject$delta)])/sum(mdsobject$dhat) 
  # configuration dissim are based on the normalized observed diss - dhat. 
  # True observed dissimilarities are delta
  rlin<-rank(lin,ties.method="random")
  path<-sapply(1:nloci,function(i)return(lin[which(rlin==i)]))
  maporder<-sapply(1:nloci,function(i)return(which(rlin==i)))
  thetalength<-path-path[1]
  chromlength<-scale*radmed*thetalength*100
  locilength<-chromlength[rlin]
  list(chromlength=chromlength,order=rlin,locilength=locilength,maporder=maporder)
}


#' Reorders a distance map by a new marker order.
#' 
#' Reorders a distance map by a new marker order.
#' 
#' @param distmap A symmetric matrix of pairwise inter-marker distances.
#' @param newrank A vector of scalars giving the new rank of each marker, markers 
#' should appear in the same order as in the distmap.
#'
#' @details
#' The rows and columns in distmap are reordered such that if entry i in \code{newrank} 
#' has value j then row j and column j in the new matrix are row i and column i 
#' from distmap.
#'
#' @return Matrix of pairwise inter-marker distances.
#'
#' @examples
#' s<-matrix(1:25,nrow=5)
#' s<-0.5*(s+t(s))
#' rank<-c(1,3,4,2,5)
#' dmap.check(s,rank)
#' @export
dmap.check<-function(distmap,newrank){
  return(outer(newrank,newrank,Vectorize(function(i,j)distmap[i,j]) ))
}


#' For a given marker finds the nearest neighbours with LOD scores > 0.
#' 
#' Finds the nearest neighbours of a marker with LOD scores > 0. 
#'
#' @param loci Scalar indicating a marker number
#' @param lodmap Symmetric matrix of pairwise LOD scores
#'
#' @details
#' The columns and rows of the matrix should be in the order corresponding to 
#' the estimated map order. The function then returns the ranks of first markers 
#' to the left and right of the marker of interest with non-zero lod scores.
#'
#' @return A vector of length 1 or 2 containing the rank of the nearest informative markers.
#' @export
get.nearest.informative<-function(loci,lodmap){
  #split matrix by loci
  neighbours<-NULL
  
  if(loci>1) {
    locileft<-lodmap[loci,(loci-1):1]
    if(length(which(locileft!=0))>0)    neighbours<-loci-min(which(locileft!=0))
  }
  if(loci<dim(lodmap)[2]){
    lociright<-lodmap[loci,(loci+1):dim(lodmap)[2]]
    if(length(which(lociright!=0))>0)  neighbours<-c(neighbours,loci+min(which(lociright!=0)))
  }
  neighbours
}

#' Calculates the nearest neighbour fit for an individual marker.
#' 
#' Calculates the nearest neighbour fit for an individual marker.
#'
#' @param loci Scalar indicating the estimated rank position of the marker.
#' @param distmap Symmetric matrix of pairwise inter-marker distances with columns
#' and rows corresponding to the estimated map order.
#' @param lodmap Symmetric matrix of pairwise lod scores with columns and rows
#' corresponding to the estimated map order.
#' @param estmap Vector of estimated marker positions.
#'
#' @details
#' The nearest neighbour fit for a marker is the sum of the difference between 
#' the observed and estimated distances between the marker and its nearest 
#' informative neighbour. A neighbour is informative if the LOD score for the 
#' inter-marker distance is non zero. This function finds the nearest markers 
#' with a non-zero LOD score (this may be one or two markers). Calculates the 
#' estimated distances between these markers and the marker of interest and 
#' returns the sum of the absolute values of the difference between the observed 
#' and estimated distances.
#'
#' @return Scalar corresponding to the difference between the observed and 
#' estimated intermarker differences.
#' @export
calc.nnfit.loci<-function(loci,distmap,lodmap,estmap){
  nns<-get.nearest.informative(loci,lodmap)
  obs<-distmap[loci,nns]
  est<-estmap[loci]-estmap[nns]
  nn.fit<-sum(abs(obs-est))
  nn.fit
}


#' Calculate the nearest neighbour fit.
#' 
#' Calculates the total, mean and individual differences between the observed and 
#' estimated distances from all loci and their nearest neighbours with non-zero 
#' LOD scores.
#'
#' @param distmap Symmetric matrix of pairwise inter-marker distances with 
#' columns and rows corresponding to the estimated map order.
#' @param lodmap Symmetric matrix of pairwise lod scores with columns and
#' rows corresponding to the estimated map order.
#' @param estmap Vector of estimated marker positions.
#'
#' @details
#' The nearest neighbour fit for a marker is the sum of the difference between
#' the observed and estimated distances between the marker and its nearest
#' informative neighbour. A neighbour is informative if the LOD score for the 
#' inter-marker distance is greater than zero. This function calculates the nearest 
#' neighbour fit for each marker and returns the fit for each point and the sum 
#' of all the fits.
#'
#' @return A list with the elements:
#' \item{fit}{Sum of the nearest neighbour fits over all markers.}
#' \item{pointfits}{Vector of nearest neighbour fits for each marker.}
#' \item{meanfit}{Mean of the nearest neighbour fits over all markers.}
#'
#' @seealso
#' \code{\link{calc.nnfit.loci}}
#' @export
calc.nnfit<-function(distmap,lodmap,estmap){
  pointfits<-unlist(lapply(1:dim(distmap)[2],calc.nnfit.loci,distmap=distmap,lodmap=lodmap,estmap=estmap))
  fit<-sum(pointfits)
  list(fit=fit,pointfits=pointfits,meanfit=mean(pointfits))
}

#' Calculates the number of swaps required to move from one order to another.
#' 
#' Calculates the number of swaps required to move from one order to another.
#'
#' @param map1 Vector of marker positions or ranks.
#' @param map2 Vector of marker positions or ranks.
#'
#' @details
#' This is intended to be used when comparing an estimated marker ordering to 
#' some perceived "truth". It is most likely to be useful when dealing with 
#' simulated data where the concept of truth makes most sense. It calculates 
#' the minimum number of single place swaps that would be needed to move from 
#' \code{map1} to \code{map2} and it does this by reverse engineering kendall's tau b 
#' correlation coefficient 
#' \deqn{\tau=\frac{2(C-D)}{N}}{\tau=2(C-D)/N}
#' where \eqn{N} is the total number of pairs of markers, \eqn{C} the number of concordant pairs and \eqn{D} the number 
#' of discordant pairs. If there are \eqn{n} markers then the total number of pairs 
#' \eqn{N={{n}\choose{2}}}{N=nC_2} and \eqn{C=N-D} so \eqn{D=0.5 {{n}\choose{2}}(1-\tau)}{D=0.5 nC_2(1-\tau)} and the minimum 
#' number of swaps is the minimum of \eqn{D} and \eqn{N-D}
#'
#' @return Scalar giving the number of swaps.
#' @export
calc.nswaps<-function(map1, map2){
  n<-length(map1)
  N<-n*(n-1)/2
  nswaps<-0.5*N*(1-stats::cor(map1,map2,method="kendall"))
  return(round(min(nswaps,N-nswaps)))
}




#' Diagnostic plots for the map estimation using calc.maps.pc with 2 dimensions.
#'
#' Diagnostic plots for the map estimation using calc.maps.pc with 2 dimensions.
#'
#' @param x Map object from \code{\link[=calc.maps.pc]{calc.maps.pc()}} with 2 dimensions.
#' @param D1lim Numeric vector specifying the limits of the horizontal axis.
#' @param D2lim Numeric vector specifying the limits of the vertical axis.
#' @param displaytext Logical argument determining how markers should be labelled 
#' in the wMDS configuration plot. If TRUE then marker names are used. If FALSE 
#' then numbers are used.
#' @param ... Further arguments are ignored. (accepted for compatibility with generic plot)
#'
#' @details
#' Plots 2 panels:
#'
#' Panel 1 the final MDS configuration and the fitted principal curve from the 
#' \code{\link[=calc.maps.pc]{calc.maps.pc()}} in 2 dimensions. If D1lim or D2lim 
#' is not specified, then limits are defined by \code{\link[smacof]{plot.smacof}}. 
#'
#' Panel 2 the pointwise nearest neighbour fits in order of the position in the 
#' estimated map.
#'
#' Markers are assigned numbers according to the order in which they occur in 
#' the input file. The locikey output of the map object is a data frame 
#' associating marker names with their numbers. This can be accessed using 
#' \code{pcmap$locikey}. If \code{displaytext=FALSE} then markers will be labelled 
#' by these numbers. By default displaytext=TRUE and markers are labelled by 
#' marker name.
#'
#' @references
#' \cite{de Leeuw J, Mair P (2009) Multidimensional scaling using majorization: SMACOF in R. J Stat Softw 31: 1-30} \url{http://www.jstatsoft.org/v31/i03/}
#'
#' @seealso
#' \code{\link{plot.pcmap3d}}, \code{\link{plot.spheremap}},\code{\link[smacof]{plot.smacof}}, \code{\link{calc.maps.pc}}
#'
#' @examples
#' map<-calc.maps.pc(system.file("extdata", "lgV.txt", package="MDSMap"),
#' ndim=2,weightfn='lod2',mapfn='haldane')
#' plot(map)
#' @export
plot.pcmap <- function (x,D1lim=NULL,D2lim=NULL,displaytext=TRUE,...){
  graphics::par(mfrow=c(1,2))
    with(x,{
      if (displaytext==TRUE) {
	    labels=locikey$locus
	  } else {
	    labels=locikey$confplotno
	  }
      graphics::plot(smacofsym$conf,type="n",main='MDS with principal curve',xlim=D1lim,ylim=D2lim,xlab='Dim 1',ylab='Dim 2')
      text(smacofsym$conf,labels=labels,cex=0.8)
      lines(pc)
      if (displaytext==TRUE)  {
	    labels1=locimap$locus
	  } else  {
	    labels1=locimap$confplotno
	  }
      graphics::plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit',main='nearest neighbour fits')
      text(locimap$position,locimap$nnfit,labels1)
    })
}


#' Diagnostic plots for the map estimation using calc.maps.pc with 3 dimensions.
#' 
#' Diagnostic plots for the map estimation using calc.maps.pc with 3 dimensions.
#'
#' @param x Map object from calc.maps.pc() with 3 dimensions.
#' @param D1lim Numeric vector specifying the limits of the axis relating to dimension 1 of the wMDS used to obtain pcmap3d.
#' @param D2lim Numeric vector specifying the limits of the axis relating to dimension 1 of the wMDS used to obtain pcmap3d.
#' @param D3lim Numeric vector specifying the limits of the axis relating to dimension 1 of the wMDS used to obtain pcmap3d.
#' @param displaytext Logical argument determining how markers should be labelled 
#' in the wMDS configuration plot. If TRUE then marker names are used. If FALSE then numbers are used.
#' @param ... Further arguments are ignored. (accepted for compatibility with generic plot)
#'
#' @details
#' Plots 4 panels
#'
#' Panels 1-3 show the final MDS configuration and the fitted principal curve from 
#' the \code{\link[=calc.maps.pc]{calc.maps.pc()}} in 3 dimensions. plots \code{D1} vs 
#' \code{D2}, \code{D1} vs \code{D3} and \code{D2} vs \code{D3}.  If \code{D1lim}, 
#' \code{D2lim} or \code{D3lim} is not specified, then limits are defined by \code{\link[smacof]{plot.smacof}}.
#'
#' Panel 4 shows the pointwise nearest neighbour fits in order of the position 
#' in the estimated map.
#'
#' Also plots a 3 dimensional scatterplot of the final MDS configuration and the 
#' fitted principal curve in a new window using \code{\link[rgl]{plot3d}} from the 
#' \code{\link[=rgl]{rgl package}}. 
#'
#' Markers are assigned numbers according to the order in which they occur in the 
#' input file. The locikey output of the map object is a data frame associating 
#' marker names with their numbers. This can be accessed using \code{pcmap3d$locikey}. 
#' If \code{displaytext=FALSE} then markers will be labelled by these numbers. 
#' By default \code{displaytext=TRUE} and markers are labelled by marker name.
#'
#' @references
#' \cite{de Leeuw J, Mair P (2009) Multidimensional scaling using majorization: SMACOF in R. J Stat Softw 31: 1-30} \url{http://www.jstatsoft.org/v31/i03/}
#' @seealso 
#' \code{\link{plot.pcmap}}, \code{\link{plot.spheremap}},\code{\link[smacof]{plot.smacof}}, \code{\link{calc.maps.pc}}, \code{\link[rgl]{plot3d}}
#'
#' @examples
#' map<-calc.maps.pc(system.file("extdata", "lgV.txt", package="MDSMap"),
#' ndim=3,weightfn='lod2',mapfn='haldane')
#' plot(map)
#' @export
plot.pcmap3d <- function (x,D1lim=NULL,D2lim=NULL,D3lim=NULL,displaytext=TRUE,...) {
  graphics::par(mfrow=c(2,2))
  with(x,{
    if (displaytext==TRUE) {
	  labels=locikey$locus 
	} else {
	  labels=locikey$confplotno
	}
    graphics::par(mfrow=c(2,2))
    graphics::plot(smacofsym$conf[,'D1'],smacofsym$conf[,'D2'],type="n",main='MDS with principal curve',xlab='Dimension 1',ylab='Dimension 2',xlim=D1lim,ylim=D2lim)
    text(smacofsym$conf[,'D1'],smacofsym$conf[,'D2'],labels=labels,cex=0.8)
    lines(pc$s[,'D1'][pc$tag],pc$s[,'D2'][pc$tag])
    graphics::plot(smacofsym$conf[,'D1'],smacofsym$conf[,'D3'],type="n",main='MDS with principal curve',xlab='Dimension 1',ylab='Dimension 3',xlim=D1lim,ylim=D3lim)
    text(smacofsym$conf[,'D1'],smacofsym$conf[,'D3'],labels=labels,cex=0.8)
    lines(pc$s[,'D1'][pc$tag],pc$s[,'D3'][pc$tag])
    graphics::plot(smacofsym$conf[,'D2'],smacofsym$conf[,'D3'],type="n",main='MDS with principal curve',xlab='Dimension 2',ylab='Dimension 3',xlim=D2lim,ylim=D3lim)
    text(smacofsym$conf[,'D2'],smacofsym$conf[,'D3'],labels=labels,cex=0.8)
    lines(pc$s[,'D2'][pc$tag],pc$s[,'D3'][pc$tag])
    if (displaytext==TRUE) {
	  labels1=locimap$locus
	} else {
	  labels1=locimap$confplotno
	}
    graphics::plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit',main='nearest neighbour fits')
    text(locimap$position,locimap$nnfit,labels1)
    rgl::plot3d(smacofsym$conf,type="n")
    rgl::text3d(smacofsym$conf,text=labels)
    lines3d(pc$s[pc$tag,])
  })
}

#' Produces diagnostic plots for the estimated map using \code{\link{calc.maps.sphere}}.
#' 
#' Produces diagnostic plots for the estimated map using \code{\link{calc.maps.sphere}}.
#'
#' @param x Map object from \code{\link{calc.maps.sphere}}
#' @param displaytext Logical argument determining how markers should be labelled in the MDS 
#' configuration plot. If TRUE then marker names are used. If FALSE then numbers are used.
#' @param ... Further arguments are ignored. (accepted for compatibility with generic plot)
#'
#' @details 
#' Produces a figure with 3 panels from a map object produced by \code{\link{calc.maps.sphere}}.
#'
#' Panel one shows the stress of the unconstrained MDS, the stress of the constrained 
#' MDS and the ratio of the two. A good rule of thumb is the stress from the constrained 
#' MDS should not be more than 10% higher than the stress of the unconstrained MDS.
#'
#' Panel 2 shows the final configuration of the unconstrained MDS which can be used to 
#' identify outliers.
#'
#' Panel 3 shows the final configuration of the constrained MDS in black and the 
#' unconstrained MDS in red. This can be used to check that the constrained fit 
#' is not distorting the data - large changes in the rank of a point in either 
#' dimension 1 or dimension 2 are indications of a problem with the fit.
#'
#' Panel 4 shows the pointwise nearest neighbour fits in order of the position 
#' in the estimated map.
#'
#' If \code{D1lim} or \code{D2lim} is not specified, then limits of panels 2 and 
#' 3 are defined by \code{\link[smacof]{plot.smacof}}. 
#'
#' Markers are assigned numbers according to the order in which they occur in the 
#' input file. The locikey output of the map object is a data frame associating 
#' marker names with their numbers. This can be accessed using \code{pcmap3d$locikey}. 
#' If \code{displaytext=FALSE} then in panels 2 and 3 markers will be labelled by 
#' these numbers. By default \code{displaytext=TRUE} and markers are labelled by 
#' marker name.
#'
#' @references
#' \cite{de Leeuw J, Mair P (2009) Multidimensional scaling using majorization: SMACOF in R. J Stat Softw 31: 1-30} \url{http://www.jstatsoft.org/v31/i03/}
#'
#' @seealso
#' \code{\link{plot.pcmap}}, \code{\link{plot.pcmap3d}}, \code{\link[smacof]{plot.smacof}}, \code{\link{calc.maps.sphere}}
#'
#' @examples
#' map<-calc.maps.sphere(system.file("extdata", "lgI.txt", package="MDSMap"),
#' weightfn='lod', mapfn='kosambi')
#' plot(map)
#' @export
plot.spheremap <- function (x,displaytext=TRUE,...) {
  
  graphics::par(mfrow=c(2,2))
  with(x,{
    if (displaytext==TRUE) {
	  labels=locikey$locus 
	} else {
	  labels=locikey$confplotno
	}
    graphics::plot(c(0,1),c(0,1),type='n',axes=F,xlab="",ylab="")
    text(0.5,0.7,paste('Sym Stress =',round(ssym,digits=4)))
    text(0.5,0.55,paste('Sphere Stress/Sym Stress =',round(stressratio,digits=4)))
    text(0.5,0.4,paste('Sphere Stress =',round(ssphere,digits=4)))

    graphics::plot(smacofsym,plot.type="confplot",type="n",main='Unconstrained',label.conf=list(label=FALSE,pos=1,col=1))
    text(smacofsym$conf,labels=labels,cex=0.8)
	xlower=min(smacofsym$conf[,1],smacofsphere$conf[,1])
    xupper=max(smacofsym$conf[,1],smacofsphere$conf[,1])
    ylower=min(smacofsym$conf[,2],smacofsphere$conf[,2])
    yupper=max(smacofsym$conf[,2],smacofsphere$conf[,2])
    graphics::plot(smacofsym,plot.type="confplot",type="n",main='Unconstrained + Spherical',label.conf=list(label=FALSE,pos=1,col=1),xlim=c(xlower,xupper),ylim=c(ylower,yupper))
    text(smacofsym$conf,labels=labels,cex=0.8)
    l=dim(smacofsphere$conf)[1]-1
    text(smacofsphere$conf[1:l+1,],labels=labels,cex=0.8,col="red") 
    if (displaytext==TRUE)  {
	  labels1=locimap$locus 
	} else {
	  labels1=locimap$confplotno
	}
    graphics::plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit',main='nearest neighbour fits')
    text(locimap$position,locimap$nnfit,labels1)
  })
}





#' Calculate a nearest neighbour fit from an estimated map object.
#' 
#' Calculates a new nearest neighbour fit based on a new order from a map object 
#' generated by \code{\link{calc.maps.pc}}, \code{\link{calc.maps.sphere}} or 
#' \code{\link{estimate.map}}
#'
#' @param estmap A character string indicating the name of a comma separated value 
#' file with  the first column containing marker names in the order of their estimated position.
#' @param mapobject A map object generated by \code{\link{calc.maps.pc}}, 
#' \code{\link{calc.maps.sphere}} or \code{\link{estimate.map}}.
#' @param header Logical argument indicating whether the .csv file \code{estmap} 
#' contains headers - default is \code{TRUE}
#'
#' @details
#' Reads in a new estimated order, reorders the distance map and LOD scores by 
#' the new order and recalculates the nearest neighbour fit.
#'
#' @return A list with the elements:
#' \item{fit}{Sum over all markers of the nearest neighbour fits.}
#' \item{pointfits}{The nearest neighbour fit for each marker.}
#' \item{meanfit}{Meanv of the nearest neighbour fits over all markers.}
#'
#' @seealso
#' \code{\link{calc.maps.pc}}, \code{\link{calc.maps.sphere}}, \code{\link{estimate.map}}, \code{\link{calc.nnfit}}
#' @export
recalc.nnfit.from.map<-function(estmap,mapobject,header=TRUE){
  estmap<-utils::read.csv(estmap,header=header)
  M<-mapobject$distmap
  lod<-mapobject$lodmap
  lnames<-colnames(M)
  names<-estmap[,1]
  maporder<- sapply(1:length(names),function(i)which(lnames==names[i]))
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder,Vectorize(function(i,j)lod[i,j]))
  nnfit<-calc.nnfit(distmap,lodmap,estmap[,2])
  nnfit
}

#' Calculates the distance of a marker from some objective "truth".
#' 
#' Calculates the distance of a marker from some objective "truth".
#'
#' @param loci Character string or number specifying the marker name.
#' @param estmap Data frame in which has a column called "names" containing marker 
#' names and a column called "position" containing marker positions.
#' @param realmap Data frame in which the first column contains marker names and 
#' the second column marker positions. Column names are not necessary.
#'
#' @details
#' Both the first column of \code{realmap} and \code{estmap$name} must contain 
#' \code{markername}, but aside from this they do not have to have identical entries. 
#'
#' @return  position in estmap-position in realmap
#' @export
get.dist.loci<-function(loci,estmap,realmap){
  l<-estmap$name[loci]
  dist<-estmap[estmap$name==l,]$position-realmap[which(realmap[,1]==l),2]
  dist
}


#' Calculates mean square distance between marker positions in two different maps.
#' 
#' Calculates mean square distance of markers in the analysis from some objective 
#' "truth".
#'
#' @param estmap Estimated map with 2 columns, \code{name} and \code{position} 
#' contain marker names and positions.
#' @param realmap Map in which the first column contains marker names and the 
#' second contains marker positions. Column names are not necessary.
#'
#' @details
#' The first column of \code{realmap} must contain identical entries to 
#'  \code{estmap$name}. However, the order of entries can be different.
#'
#' For every marker the difference between the position stated in \code{estmap} 
#' and in \code{realmap} is calculated (see \code{\link{get.dist.loci}}).
#'
#' Every difference is squared and the mean of the square differences is returned.
#'
#' Note that where different weights are used in estimating maps, it is valid to 
#' compare the mean distance from the truth. However, if different map functions 
#' are used then the distances are not comparable. 
#'
#' Therefore, if there is some knowledge of markers on a chromosome and data is 
#' simulated so that there is some objective knowledge of the truth then this 
#' function could be used to decide whether to use \code{lod} or \code{lod2} 
#' weightings to estimate maps attempting to locate additional markers. However, 
#' it is not suitable for deciding on the map functions used to calculated the 
#' pairwise marker distances. 
#'
#' @return A list with the following elements:
#' \item{pointdist}{Data frame containing marker names and the distance between 
#' the estimated position and the "real" position.}
#' \item{meansquaredist}{mean square distance between the estimated real position 
#' of markers.}
#' @export
meandist.from.truth<-function(estmap,realmap){
  dist<-unlist(lapply(1:dim(estmap)[1],get.dist.loci,estmap=estmap,realmap=realmap))
  meansquaredist<-mean(dist^2)
  pointdist<-cbind(name=estmap$name,dist=dist)
  list(pointdist=pointdist, meansquaredist=meansquaredist)
}



#' Invert the order of locimap from an estimated map.
#' 
#' Takes the locimap from estimate.maps a dataframe containing names and positions 
#' and any other information in increasing order of distance and inverts the order.
#' 
#' @param locimap a data frame containing the markers names and positions
#' 
#' @details
#' The map should be a data frame with a column called 'position'. It should have a starting marker a position zero. 
#' The function then inverts the distances from so that the marker at maximum distance from the starting marker (the end marker) is at distance 0 
#' and the original starting marker is now at the maximum distance. It also inverts the order of the rows in the data frame. 
#' Thus if the markers were originally in order of increasing distance from the starting marker they will now be in order of increasing distance 
#' from the end marker.
#' @return The original data frame in inverted order with the distances inverted so that the end marker is now the starting marker. 
#' @export
invert.map<-function(locimap){
  neworder<-dim(locimap)[1]:1
  newposition=max(locimap$position)-locimap$position
  locimap$position<-newposition
  return(locimap[neworder,])
}



 
#' Simulate a backcross population from homozygous parents.
#'
#' Simulates a backcross population from homozygous parents and writes a file 
#' containing the number of markers and observed pairwise distances, the pairwise 
#' recombination fractions and LOD scores in a text file suitable for analysis by 
#' other functions in the package.
#'
#' @param fname	a character string specifying the base name of the file fname.txt to which the data should be written
#'
#' @details
#'This function simply generates data for use with the vignette. The R/qtl package is used to simulate a backcross #'population of 200 individuals from homozygous parents with 200 markers in a single linkage group of length #'100cM. The recombination fractions and LOD scores are calculated. The data is written to a text file in the #'format of output from JoinMap 4.  In particular, the data is cast into a data frame with marker names in the #'first two columns, pairwise recombination fractions in the third column and associated LOD scores in the fourth #'column. The data is written to a text file 'fname.txt' where the first row contains two entries - the number of #'markers and the number of pairwise observations. Below this the data frame containing the distance data is #'appended with no column headings. 
#' \tabular{llll}{
#' \cr
#' \code{nmarkers} \tab \tab \tab \cr
#' \code{marker_1} \tab \code{marker_2} \tab \code{recombination fraction} \tab \code{lod}\cr
#' \code{1} \tab \code{2} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{1} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{3} \tab \code{.} \tab \code{.} \cr
#' \code{2} \tab \code{4} \tab \code{.} \tab \code{.} \cr
#' \code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
#' }
#'
#' @return No output - just the text file as above
#'
#' @references 
#' \cite{Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics. 189: 889-890}
#' \cite{Van Ooijen JW (2006) JoinMap 4; Software for the calculation of genetic linkage maps in experimental populations. Wageningen; Netherlands: Kyazma B.V}
#' @export
sim.bc.rflod.file<-function(fname){
  mymap<-qtl::sim.map(100,n.mar=200,eq.spacing=FALSE) 
  crossmap<-qtl::sim.cross(mymap,n.ind=200,type="bc",error.prob=0.01,map.function="haldane")
  popdata<-crossmap$geno$X$data
  nloci<-dim(popdata)[2]
  locinames<-colnames(popdata)
  l<-dim(popdata)[1]
  d<-dim(popdata)[2]
  nrmat<-outer(1:d,1:d,FUN=Vectorize( function(i,j) sum(popdata[,i]==popdata[,j])))
  colnames(nrmat)<-locinames
  rownames(nrmat)<-locinames
  temp<-nrmat
  temp[upper.tri(temp,diag=TRUE)]<-NA
  temp<-as.data.frame(as.table(temp))
  temp<-temp[!is.na(temp[,3]),]
  rf.frame<-data.frame(marker1=temp[,2],marker2=temp[,1],rf=(l-temp[,3])/l)
  rf.frame$lod=log((1-rf.frame$rf)^temp[,3]*rf.frame$rf^(l-temp[,3])/(0.5^l),10)
  fn<-paste(fname,".txt",sep="")
  utils::write.table(t(c(d,dim(rf.frame)[1])),file=fn,row.names=FALSE,col.names=FALSE,sep="\t")
  utils::write.table(rf.frame,file=fn,append=TRUE, quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")}


