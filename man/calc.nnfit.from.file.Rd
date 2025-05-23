% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mds.R
\name{calc.nnfit.from.file}
\alias{calc.nnfit.from.file}
\title{Nearest neighbour fit from estimated map and file of pairwise recombination fractions.}
\usage{
calc.nnfit.from.file(
  estmap,
  fname,
  mapfn = "haldane",
  n = NULL,
  header = FALSE
)
}
\arguments{
\item{estmap}{A character string indicating the name of a comma separated 
value file with  the first column containing marker names in the order 
of their estimated position.}

\item{fname}{A character string specifying the base name of the file 
\code{fname.txt} which contains the data to be analysed the file should be 
white space or tab separated.}

\item{mapfn}{Character string, \code{'haldane'}, \code{'kosambi'} or \code{'none'}
specifying the values to use to estimate the map distance from the 
recombination fractions. Default is \code{'haldane'}.}

\item{n}{Vector of character strings or numbers specifying the markers to be 
omitted from the analysis. Default is \code{NULL}.}

\item{header}{Logical argument indicating whether the .csv file estmap contains headers - default is TRUE}
}
\value{
A list with the following elements:
\item{fit}{Sum over all markers of the nearest neighbour fits.}
\item{pointfits}{The nearest neighbour fit for each marker.}
\item{meanfit}{Mean of the nearest neighbour fits over all markers.}
}
\description{
Calculates a nearest neighbour fit based from an estimated map and a file 
containing pairwise recombination fractions and LOD scores.
}
\details{
Reads in two files \code{fname.txt} and \code{estmap}. 

The data is cast the data into symmetric matrices of pairwise recombination 
fractions and LOD scores with the order of columns and rows in the matrix 
determined by the order specified in \code{estmap}. A distance matrix is 
calculated according to the method specified by \code{mapfn}. Haldane is the 
default map function, None just uses recombination fractions and the other 
alternative is Kosambi (see \code{\link{dmap}} for details). The nearest 
neighbour fit is then calculated (see \code{\link{calc.nnfit}} for details)

\code{estmap} should contain marker names in the first column in the order 
of the estimated map.

\code{fname} should be of the form \code{fname.txt} and it is assumed that they are in 
a tab or space separated file of the format displayed below. The first entry on 
the first row is the number of markers to be analysed. Underneath this is a 
table in which the first two columns contain marker names, the third column 
contains the pairwise recombination fractions between the markers and the 
fourth column the associated LOD score.  Note that marker names in the first 
column vary more slowly than in the second column. Missing recombination pairs 
are acceptable. Recombination fractions greater than 0.499999 are set to that 
value.

\tabular{llll}{
\cr
\code{nmarkers} \tab \tab \tab \cr
\code{marker_1} \tab \code{marker_2} \tab \code{recombination fraction} \tab \code{LOD}\cr
\code{1} \tab \code{2} \tab \code{.} \tab \code{.} \cr
\code{1} \tab \code{3} \tab \code{.} \tab \code{.} \cr
\code{1} \tab \code{4} \tab \code{.} \tab \code{.} \cr
\code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
\code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
\code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
\code{2} \tab \code{3} \tab \code{.} \tab \code{.} \cr
\code{2} \tab \code{4} \tab \code{.} \tab \code{.} \cr
\code{.} \tab \code{.} \tab \code{.} \tab \code{.} \cr
}
}
\seealso{
\code{\link{dmap}}, \code{\link{calc.nnfit}}, \code{\link{calc.pair.rf.lod}}
}
