\name{agraph-class}
\docType{class}
\alias{agraph-class}
\alias{lparm,agraph-method}
\alias{gparm,agraph-method}
\alias{lparm<-,agraph-method}
\alias{gparm<-,agraph-method}

\title{Class \sQuote{agraph}}
\description{Anchor Graph Functions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("agraph", ...)}.
or by calling the \code{\link{agraph}} function.
}
\section{Slots}{
  \describe{
    \item{lparm}{Object of class \sQuote{parm} contains
      lambda parameter values used.}
     \item{gparm}{Object of class \sQuote{parm} contains the
      gamma parameter value used.}
  }
}
\section{Methods}{
  \describe{
    \item{lparm}{\code{signature(object = "agraph")}: returns the lambda 
    	parameter values used. }
   \item{gparm}{\code{signature(object = "agraph")}: returns the
      gamma parameter value used.}
  }
}

\author{Mark Vere Culp}


\seealso{
  \code{\link{agraph}}
}
\keyword{classes}
