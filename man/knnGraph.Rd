\name{knnGraph}
\alias{knnGraph}
\alias{epsGraph}
\alias{cosineDist}
\alias{euclidDist}
\alias{kgraph.predict}
\alias{kgraph.symm}
\alias{eps.predict}
\alias{eps.symm}
\title{
  Convert a \sQuote{data.frame} or \sQuote{matrix} into a
  k-NN graph or epsilon graph.
}
\description{
  These functions are designed to convert data sets or distance matrices
  to graphs. In addition, the Cosine and Euclidean distance functions are also available. 
  The kgraph functions build a knn graph off of a distance matrix (either square (symm) or non-square (predict))
}
\usage{
knnGraph(x,y,k=6L,nok=FALSE,metric="cosine")
epsGraph(x,y,eps=0.2,noeps=FALSE,metric="cosine")
cosineDist(x,y=NULL)
euclidDist(x,y=NULL)
}
\arguments{
 \item{x}{ an object of type \sQuote{matrix}.}
  \item{k}{the lower bound to the number of neighbors for a k-NN graph.}
  \item{eps}{neighborhood radius for an epsilon graph.}
  \item{nok}{the parameter k is ignored.}
  \item{noeps}{the parameter eps is ignored.}
  \item{metric}{the metric to used to compute distances.}
  \item{y}{ a m x p \sQuote{matrix}. If NULL distances
  	between argument x and x are computed.}
}
\author{Mark Vere Culp}
\seealso{ \code{\link{s4pm}}, \code{\link{jtharm}}}
\keyword{classes}
\keyword{models}
\keyword{methods}
