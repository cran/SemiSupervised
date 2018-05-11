\name{s4pm-class}
\docType{class}
\alias{s4pm-class}
\alias{lparm}
\alias{hparm}
\alias{gparm}
\alias{lparm<-}
\alias{hparm<-}
\alias{gparm<-}
\alias{lparm,s4pm-method}
\alias{hparm,s4pm-method}
\alias{gparm,s4pm-method}
\alias{lparm<-,s4pm-method}
\alias{hparm<-,s4pm-method}
\alias{gparm<-,s4pm-method}

\title{Class \sQuote{s4pm}}
\description{Safe Semi-Supervised Semi-Parametric Model Class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("s4pm", ...)}.
or by calling the S4 generic \code{\link{s4pm}} function.
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
    \item{lparm}{\code{signature(object = "s4pm")}: returns the lambda 
    	parameter values used.}
    \item{hparm}{\code{signature(object = "s4pm")}: returns the h parameter for the 
    	kernel function used.}
    \item{gparm}{\code{signature(object = "s4pm")}: returns the
      gamma parameter value used.}
  }
}

\author{Mark Vere Culp}


\seealso{
  \code{\link{s4pm}}
}
\keyword{classes}
