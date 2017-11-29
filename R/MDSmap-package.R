#' High density Genetic Linkage Mapping using Multidimensional Scaling
#'
#' MDSmap provides functions for estimating genetic linkage maps for
#' markers from a single linkage group from pairwise intermarker map
#' distances using the Haldane or Kosambi map function; or recombination 
#' fractions. It either uses constrained weighted metric multidimensional 
#' scaling (cMDS) in 2 dimensions or unconstrained weighted metric 
#' multidimensional scaling (MDS) followed by fitting a principal curve 
#' (PC) in either 2 or 3 dimensions. Pairwise distances can be weighted 
#' either by the LOD score or LOD2. There are functions for diagnostic 
#' plots, estimating the difference between the observed and estimated 
#' difference between points and their nearest informative neighbour, 
#' which may be useful in deciding which weights to use and also for 
#' testing estimated maps against a map estimated externally.
#'
#' The main top level functions to use: \code{\link{calc.maps.pc}} and 
#' \code{\link{calc.maps.sphere}}, and use \code{\link{plot.pcmap}}, 
#' \code{\link{plot.spheremap}} or \code{\link{plot.pcmap3d}} to visualize 
#' the result.
#'
#' @name MDSMap-package
#' @author Katharine F. Preedy <Katharine.preedy@bioss.ac.uk>
#' @import smacof princurve rgl reshape
#' @examples
#' map<-calc.maps.pc(system.file("extdata", "lgV.txt", package="MDSMap"),
#' ndim=2,weightfn='lod2',mapfn='kosambi')
#' plot(map)
#' 
#'
NULL
