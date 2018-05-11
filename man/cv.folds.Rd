\name{cv.folds}
\alias{cv.folds}
\alias{extract.a}
\alias{extract.xg}
\alias{extract.g}
\alias{h.est}
\alias{inter.predict}
\title{Exported Internal Functions}
\description{
  These functions are internal functions that are exported only so other
  packages can use them to avoid the ::: scope operator. Users are
  strongly encouraged not to use them directly or modify them.
  
  The documentation is minimal to discourage people from further using
  them.
}
\usage{
cv.folds(n, folds = 3L)
h.est(Dist,sym=TRUE,thresh)
extract.a(form,env)
extract.xg(form,env)
extract.g(form,env)
inter.predict(y,W,x=NULL,bet=NULL,pow,fhat=FALSE)
}
\arguments{
  \item{n}{a number greater than folds.}
  \item{folds}{the number of folds.}
  \item{Dist}{a distance matrix.}
  \item{sym}{remove symmetric matrix component below diagonal.}
  \item{thresh}{an internal thresh-hold parameter.}
  \item{form}{a \code{formula}.}
  \item{env}{an \code{environment}.}
  \item{y}{a response for interpolation.}
  \item{W}{a matrix.}
  \item{x}{corresponding feature data.}
  \item{bet}{a default coefficient vector.}
  \item{pow}{the power of the interpolant function. Should
    be ignored.}
  \item{fhat}{whether or not to compute the sum of the linear
    term and the graph term.}
}
\note{
  All of these functions should be treated as internal exported
  functions. The h.est function was primarly taken from the \code{kernlab}
  package. 
}
\author{Mark Vere Culp}

\keyword{classes}
\keyword{models}
\keyword{methods}
