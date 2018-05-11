\name{SemiSupervised-class}
\docType{class}
\alias{SemiSupervised-class}
\alias{type}
\alias{gmatrix}
\alias{xmatrix}
\alias{ymatrix}
\alias{type}
\alias{measures}
\alias{type<-}
\alias{gmatrix<-}
\alias{xmatrix<-}
\alias{ymatrix<-}
\alias{measures<-}
\alias{dim}
\alias{fitted}
\alias{type,SemiSupervised-method}
\alias{gmatrix,SemiSupervised-method}
\alias{xmatrix,SemiSupervised-method}
\alias{ymatrix,SemiSupervised-method}
\alias{measures,SemiSupervised-method}
\alias{gmatrix<-,SemiSupervised-method}
\alias{xmatrix<-,SemiSupervised-method}
\alias{ymatrix<-,SemiSupervised-method}
\alias{measures<-,SemiSupervised-method}
\alias{type<-,SemiSupervised-method}
\alias{dim,SemiSupervised-method}
\alias{fitted,SemiSupervised-method}
\alias{parm}
\alias{parm-class}

\title{Class \sQuote{SemiSupervised}}
\description{S4 VIRTUAL class used as a base for the various classes in the \pkg{SemiSupervised} and \pkg{spa} packages.
	
Object \sQuote{parm} is the class union of \sQuote{vector},\sQuote{numeric}, \sQuote{integer}, and \sQuote{NULL} for parmeters and has no methods or slots.
}
\section{Objects from the Class}{  
  Objects from the class cannot be created directly but only contained
  in other classes.
}
\section{Slots}{
  \describe{
    \item{type}{Object of class \sQuote{parm} contains
      lambda parameter values used.}
    \item{gmatrix}{Object of class \sQuote{input} for the data
      graph used during computations.}
    \item{xmatrix}{Object of class \sQuote{output} for the data
      matrix used during computations.}
    \item{ymatrix}{Object of class \sQuote{output} for the response 
    	matrix/vector.}
    \item{measures}{Object of class \sQuote{output} that provides 
    	various measures for the fit.}
    \item{fit}{The raw fitted values for the n observations.}
  \item{.control}{Internal storage of control values. }
  \item{.terminfo}{Internal storage of terms object. }
	\item{.call}{Internal call storage.
		 S4 generics do not display the call exactly as called by the user, 
		 which is well known.}
	\item{.fitinfo}{Internal storage of fit information. }
	\item{.cv_str}{Internal storage of CV information. }
		\item{.respinfo}{Internal storage of response information.}
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(object = "SemiSupervised")}: returns the dimensions of 
    	the feature data for an \sQuote{SemiSupervised} object. }
   \item{fitted}{\code{signature(object = "SemiSupervised")}: returns the fitted values
   	of any \sQuote{SemiSupervised} object. }
   \item{type}{\code{signature(object = "SemiSupervised")}: returns the type of fit
    	 called. }
    \item{measures}{\code{signature(object = "SemiSupervised")}: returns fit statistics. }
    \item{gmatrix}{\code{signature(object = "SemiSupervised")}: returns the adjacency
      matrix for the fit.}
    \item{xmatrix}{\code{signature(object = "SemiSupervised")}: returns the data
      matrix if used or NULL if x data was not used.}
    \item{ymatrix}{\code{signature(object = "SemiSupervised")}: returns the
      response vector.}
  }
}

\author{Mark Vere Culp}


\seealso{
  \code{\link{s4pm}}, \code{\link{jtharm}}, \code{\link{agraph}} 
}
\keyword{classes}
