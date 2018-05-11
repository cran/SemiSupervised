\name{x.scaleL}
\alias{x.scaleL}
\title{Scale a data set in accordance to a labeled index set}
\description{
	This method centers and scales the columns of a numeric matrix using the means and columns derived from a subset of
	its indices.	
 }
\usage{
	x.scaleL(x,L,sanity=TRUE,sanity.only=FALSE)
}
\arguments{
  \item{x}{a numeric object with n rows and p columns.}
  \item{L}{the index set to define the means and variances.}
  \item{sanity}{a flag to indicate whether \code{x} should be converted to a \sQuote{matrix} form.}
  \item{sanity.only}{this flag bypasses the scaling but does coerce \code{x} to a \sQuote{matrix}.}
}
\note{
This function is primarily an internal function used to make sure that the feature data are entered correctly and has the form expected by the underlying procedures. However, there are circumstances where one should call this directly for data preprocessing. This is why it is released externally. 
}
\author{Mark Vere Culp}

\keyword{classes}
\keyword{models}
\keyword{methods}
