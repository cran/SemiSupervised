\name{SemiSupervised.control}
\alias{SemiSupervised.control}
\title{Control Parameters for the S4-generic generic functions
	derived from virtual class \code{SemiSupervised}}
\description{Controls various aspects of fitting any \sQuote{SemiSupervised} object.}
\usage{             
SemiSupervised.control(normalize=TRUE,stability=NULL,k=NULL,nok=FALSE,
dissimilar=TRUE,l.eps=1e-5,l.thresh=25L,h.thresh=1e-5,U.as.anchor.thresh=600L,
U.as.anchor=TRUE,sig.est=TRUE,sig.frac=0.5,iter.max=1000L,
anchor.seed=100,sfac=5L,cn=4L,LAE.thresh=100L,LAE.eps=1e-4,
cv.fold=3L,cv.seed=100L,cv.cl=TRUE,cv.type="scv",cv.adjust=0.001)
}
\arguments{
\item{normalize}{flags whether or not the normalized or combinatorial graph
	operator should be used. This is not used by \sQuote{jtharm}. Further, it
	has no effect on anchor graph based approaches.}
\item{stability}{stabilization parameter for necessary inverses. A NULL value
	allows it to be set internally (recommended) but can be set manually.}
\item{k}{in fitting a distance graph with the \sQuote{formula} as y~.
	this will set the default k-NN graph parameter. It is not used otherwise.
	In the case of an anchor graph the k parameter is the number of anchors used
	by k-means.}
\item{nok}{flags the y~. \sQuote{formula} call to either fit or not fit a k-NN parameter
	which overrides the k argument.}
\item{dissimilar}{flags whether the graph is similar or dissimilar. This is not
	necessary in the \sQuote{formula} call since a call to \code{\link{sG}} in the formula
	automatically flags the graph as similarity or a call to \code{\link{dG}} to dissimilar.
	This is necessary in all non-formula calls with similarity graphs.}
	\item{sfac}{the \sQuote{s} parameter for the LAE method.}
\item{cn}{the \sQuote{cn} parameter for the LAE method.}
\item{LAE.thresh}{thresh hold for LAE algorithm.}
\item{LAE.eps}{convergence tolerance for LAE algorithm.}
\item{iter.max}{maximum number of iterations for \code{\link{kmeans}}.}
\item{anchor.seed}{sets the seed for the \code{\link{kmeans}} algorithm.}
\item{sig.est}{use an internal estimation scheme to estimate the parameter for 
	an rbf kernel. }
\item{sig.frac}{when \sQuote{sig.est} is true, the fraction of training data
	            necessary for the computaton is used. }
\item{l.thresh}{max iteration parameter for the underlying logistic 
  		regression algorithm when fitting classification. }
\item{l.eps}{threshold for stopping the underlying logistic regression
		 algorithm when fitting classification. }
\item{U.as.anchor}{if n>U.as.anchor.thresh then the anchor points are fitted as 
	the unlabeled cases to speed up the approach. This only works in the formula call
	where the graph is unspecified, i.e., y~.}
\item{U.as.anchor.thresh}{threshold for determining when the unlabeled are to be anchors.}
\item{h.thresh}{minimum allowable choice for h in the grid used by CV.}
\item{cv.fold}{the number of folds used for K-fold CV.}
\item{cv.seed}{the seed to generate the folds for K-fold CV.}
\item{cv.cl}{forces in classification at least two distinct responses
    in each fold and should be set to TRUE.}
\item{cv.type}{use "scv" or "cv".  The "scv" is faster and should be
    used.}
\item{cv.adjust}{forces an inverse stabilization in "scv".}
}
\author{Mark Vere Culp}

\keyword{classes}
\keyword{models}
\keyword{methods}
