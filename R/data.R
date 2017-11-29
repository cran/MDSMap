#' Dataset lgV.txt: pairwise recombination fractions for 238 markers.
#'
#' A dataset containing the pairwise recombinations fractions for 238 SNP markers 
#' from linkage group V of potato. These are derived from the genotypes of 190 
#' offspring from a cross between potato cultivar Stirling and the breeding line 
#' 12601ab1. Further details are available in Hackett et al. (2013).
#'
#' @name lgV.txt
#'
#' @section lgV.txt
#' 
#' @format An ascii text file in the format described in \code{\link{calc.maps.pc}} 
#' and \code{\link{calc.maps.sphere}}. The first line contains the number of markers
#' and the number of combinations. Then follow the space-separated combinations with
#' their recombination fractions and LOD scores:
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
#' @examples 
#' system.file("extdata", "lgV.txt", package="MDSMap")
#'
#' @source \cite{Hackett, C.A., McLean, K. and Bryan, G.J. (2013). Linkage analysis and QTL mapping using SNP dosage data in a tetraploid potato mapping population. PLoS ONE 8, e63939}
NULL

#' Dataset lgI.txt: pairwise recombination fractions for 143 markers.
#'
#' A dataset containing the pairwise recombinations fractions for 143 SNP markers 
#' from linkage group I of potato. These are derived from the genotypes of 190 
#' offspring from a cross between potato cultivar Stirling and the breeding line 
#' 12601ab1. Further details are available in Hackett et al. (2013).
#'
#' @name lgI.txt
#'
#' @section lgI.txt
#' 
#' @format An ascii text file in the format described in \code{\link{calc.maps.pc}} 
#' and \code{\link{calc.maps.sphere}}. The first line contains the number of markers
#' and the number of combinations. Then follow the space-separated combinations with
#' their recombination fractions and LOD scores:
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
#' @examples 
#' system.file("extdata", "lgI.txt", package="MDSMap")
#'
#' @source \cite{Hackett, C.A., McLean, K. and Bryan, G.J. (2013). Linkage analysis and QTL mapping using SNP dosage data in a tetraploid potato mapping population. PLoS ONE 8, e63939}
NULL
