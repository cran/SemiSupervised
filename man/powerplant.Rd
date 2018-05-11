\name{PowerPlant}
\alias{powerplant}
\docType{data}
\title{ The Combined Cycle Power Plant Data Set}
\description{
This dataset has 9568 records collected over 6 years (2006 - 2011). The data come from a Combined Cycle Power Plant (CCPP), which has operated under a full load. The data are available on the UCI repository, but not in an easy way to read into R. This provides a convenient form to access the data. A CCPP is composed of gas turbines (GT), steam turbine (ST), and heat recovery STs. The vacuum has effect on ST, while the others affect GT performance.

The authors of the orginal data performed a 5x2 fold statistical test on 5 shuffles (random Permutations). The first shuffle is in the order that this information is presented. The remaining shuffle indices are provided as the shuffle attribute. These shuffles preserve the original data generators 2-fold CV format for 10 measuresments which they used for statistical testing.
}
\usage{
data(powerplant)
}
\format{
 \code{AT}: Temperatures in the range of 1.81 to 37.11 Celsius.
 
 \code{V}: Exhuast vacuum in the range of 25.36 to 81.56 cm Hg.
 
 \code{AP}: Ambient presure in the range of 992.89 to 1033.30 milibar.
 
 \code{RH}: Relative humidity in the range of 25.56\% to 100.16\%.
 
 \code{PE}: Net Hourly electorial energy output 40.26 to 4895.76 MW.
 
 \code{shuffle2}: CV shuffle 2 indices (refer below for an example).
 
 \code{shuffle3}: CV shuffle 3 indices.
 
 \code{shuffle4}: CV shuffle 4 indices.
 
 \code{shuffle5}: CV shuffle 5 indices.
}
\source{The data was taken direclty from the Excel files on the UCI repository. This allows for a convienent access to these data.}
\note{The averages are taken from various sensors located around the plant that record the ambient variables every second. The variables are given without any normalization.}
\examples{
\donttest{
	data(powerplant)
	x<-powerplant[attr(powerplant,"shuffle")[,1],]
}
}
\keyword{datasets}

