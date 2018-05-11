\name{jtharm-class}
\docType{class}
\alias{jtharm-class}
\alias{lparm,jtharm-method}
\alias{hparm,jtharm-method}
\alias{gparm,jtharm-method}
\alias{lparm<-,jtharm-method}
\alias{hparm<-,jtharm-method}
\alias{gparm<-,jtharm-method}

\title{Class \sQuote{jtharm}}
\description{Joint Harmonic Functions}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("jtharm", ...)}.
or by calling the \code{\link{jtharm}} function.
}
\section{Slots}{
  \describe{
   \item{lparm}{Object of class \sQuote{parm} contains
      lambda parameter values used.}
    \item{hparm}{Object of class \sQuote{parm} contains
      h parameter for the kernel function used. }
    \item{gparm}{Object of class \sQuote{parm} contains the
      gamma parameter value used.}
  }
}
\section{Methods}{
  \describe{
    \item{lparm}{\code{signature(object = "jtharm")}: returns the lambda 
    	parameter values used. }
    \item{hparm}{\code{signature(object = "jtharm")}:returns the h parameter for the 
    	kernel function used}
    \item{gparm}{\code{signature(object = "jtharm")}: returns the
      gamma parameter value used.}
  }
}

\author{Mark Vere Culp}


\seealso{
  \code{\link{jtharm}} 
}
\keyword{classes}
