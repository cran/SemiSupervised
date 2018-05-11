\name{impute.median}
\alias{impute.median}
\title{Median imputation for NA's'}
\description{This function deals with NA for any object inherited by 
	virtual class \dQuote{SemiSupervised} by inputting
  the median of the non-missing values for each missing value.}
 
\usage{
impute.median(object,\dots)
}
\arguments{
  \item{object}{ an \R object for imputation of feature data for
   any \dQuote{SemiSupervised} based objects.}
  \item{\dots}{ignored by this function.}
}
\details{
  This function is primarily designed to be used with an \code{SemiSupervised}
  sub-object to filter missing values for the feature data. This is
  challenging since missing values in the response flag unlabeled
  cases. So missing data in the response is treated differently than
  missing data in the feature data. In order to do this, \code{na.pass} is
  forced as the na.action and cannot be changed. After an initial
  \code{na.pass} for the full data (response and feature) to handle missing
  values in the feature data, we force an additional execution of
  \code{impute.median} on the feature data only. In the case of \code{formula} y~., we
  do this prior to constructing the graph, and all these steps are
  internalized.

  The na.action for the \dQuote{SemiSupervised} based S4 generic functions is fixed for usage
  to \code{na.pass} and \code{impute.median}. So, any other imputation other than the
  median should be done outside the function directly on the data, and the response should 
  not be imputed.  
}
\keyword{classes}
\keyword{models}
\keyword{methods}
