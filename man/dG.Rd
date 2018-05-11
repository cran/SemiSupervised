\name{dG}
\alias{dG}
\alias{sG}
\alias{aG}
\title{Specify graph terms for \sQuote{formula} instances of objects contained in \dQuote{SemiSupervised}}
\description{
  A symbolic wrapper to indicate a graph term in a \code{formula} call to either
  the \code{s4pm}, \code{jtharm}, \code{agraph}, or a call to the \code{spa} from 
  the \pkg{spa} package.
  
  The \code{dG} stands for a dissimilarity graph matrix, \code{sG} stands for
  a similarity graph matrix, and \code{aG} stands for an \sQuote{anchor} graph 
  object. }
\usage{
	dG(x, k = 6L, nok = FALSE, metric = NULL)
	sG(x)
	aG(x)
}
\arguments{
  \item{x}{a symbolic representation depending on the context of call. For \code{dG} or
  	\code{sG} this is either a symmetric \sQuote{data.frame} or \sQuote{matrix}. For \code{aG} this must
  	be an \sQuote{anchor} object. Refer to details below.}
  \item{k}{the k-NN graph parameter for the \code{\link{knnGraph}} call within the
  \code{dG} function.}
  \item{metric}{the metric used to compute distances.}
  \item{nok}{a parameter to treat the graph input as the final adjacency matrix, 
  	i.e., no k-NN graph is computed. This bypasses the k parameter.}
}
\details{
	The y~. \sQuote{formula} case automatically accounts for graph terms
	based on their respective function calls, but there are cases where one must 
	specify the graph terms directly.

	In a \sQuote{formula} call with direct graph input (i.e., not y~.), the
	graph must be specified based on the original function call and the 
	type of graph desired. There are three cases.  
	
	1) For the anchor graph, \code{\link{agraph}}, the \code{aG} function must be invoked
	   and the graph must be of class \sQuote{anchor}. The only way to create this object
	   is to use the \code{\link{AnchorGraph}} function. Refer to \code{\link{AnchorGraph}}
	   for an example.
	
	2) The \code{dG} function passes or creates a dissimilarity graph (i.e., edges
	   correspond to dissimilarity with \sQuote{0} as close and \sQuote{Inf} as far).
	
	3) The \code{sG} passes a similarity graph  (i.e., edges corresponding to dissimilarity
	   with \sQuote{1} (typically) as close and \sQuote{0} as far). 
	
  Examples of each case are provided below.
	
  These commands are designed to work in the \sQuote{formula} instance of the following S4 generics:
  \code{\link{agraph}}, \code{\link{s4pm}}, \code{\link{jtharm}}, and also the \code{spa} from the
  \pkg{spa} package.
}
\author{Mark Vere Culp}
\examples{
#######
## Equivalent uses of the formula and default s4pm call.
#######

## Set up Sonar data with 20\% labeled
library(mlbench)
data(Sonar)

n=dim(Sonar)[1]
p=dim(Sonar)[2]

nu=0.2
set.seed(100)
L=sort(sample(1:n,ceiling(nu*n)))
U=setdiff(1:n,L)

y.true<-Sonar$Class
Sonar$Class[U]=NA
control=SemiSupervised.control(stability=0.0) ## turn off adjustment parameters for comparison

(g.s4pm<-s4pm(Class~.,data=Sonar,control=control)) ### Fit s4pm to Sonar

##The following give equivalent output to the g.s4pm<-s4pm(Class~.,data=Sonar) call.
\donttest{
D11=as.matrix(cosineDist(x.scaleL(Sonar[,-p],L)))
(g.s4pm1<-s4pm(Class~.+dG(D11),data=Sonar,control=control))

#######
## Equivalent uses of the formula, Class ~ ., and default jtharm call.
#######
control=SemiSupervised.control(stability=0.0) 
(g.jtharm1<-jtharm(Class~.,data=Sonar,control=control))
(g.jtharm2<-jtharm(Class~dG(D11),data=Sonar,control=control))

}
}
\keyword{classes}
\keyword{models}
\keyword{methods}
