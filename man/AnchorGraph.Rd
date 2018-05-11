\name{AnchorGraph}
\alias{AnchorGraph}
\alias{getAnchor}
\title{
  Generate an Anchor Graph from an n x p data matrix
}
\description{
  The AnchorGraph function is designed for:
    
  1) Graph generation: creates an \sQuote{anchor} graph object.
  
  2) Graph modification: modify an existing \sQuote{anchor} graph using the fit.g argument, which
  is used primarily for predicting new observations.

  The getAnchor function is a wrapper to kmeans clustering.
}
\usage{
AnchorGraph(x,metric="cosine",anchor=NULL,fit.g=NULL,control=SemiSupervised.control())
getAnchor(x,control)
}
\arguments{
  \item{x}{the n x p \sQuote{vector}, \sQuote{matrix} or \sQuote{data.frame}.}
  \item{fit.g}{an existing \sQuote{anchor} object to be updated
  	   with new observations given by argument x. }
  \item{anchor}{an optional matrix of user provided anchor points. If NULL then \code{link{kmeans}} 
  	is used to generate the anchors.}
  \item{metric}{the metric either cosine or sqDist for graph construction.}
  \item{control}{the \code{\link{SemiSupervised.control}} is used to initialize various LAE settings.
  	 This control should also be used when executing the \code{\link{agraph}} function but is
  	 not enforced. }
}
\value{
\item{Z}{the n x k Z-matrix where each row has at most sfrac (refer to \code{\link{SemiSupervised.control}}) non-zero entries.}
\item{rL}{reduced Laplacian matrix.}
\item{g.scaling}{the scaling used to scale the x data prior to graph construction.}
\item{anchor}{the anchor points.}
\item{metric}{the metric.}	
}
\details{
  The AnchorGraph function creates an informal S3-object of class \sQuote{anchor}. This is required as
  input to the S4 generic function \code{\link{agraph}} when using the \sQuote{anchor} version or \code{\link{agraph.default}}
  version.
  
  This call is performed in the \sQuote{formula}, \sQuote{data.frame}, \sQuote{matrix}, or \sQuote{vector} S4 generic instances
  of \code{\link{agraph}}.
}
\note{
  An S3-object was sufficient for our purposes. There is no need for the overhead or flexibility of a S4-class which
  is why it is programmed this way.
}

\examples{
## Set up Sonar data with 20% labeled
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

g.agraph1<-agraph(Class~.,data=Sonar) 

##The following gives an equivalent output to the g.agraph1<-agraph(Class~.,data=Sonar) call.
\donttest{
ctrl<-SemiSupervised.control()
g<-AnchorGraph(x.scaleL(Sonar[,-p],L),control=ctrl)
g.agraph2<-agraph(Class~.+aG(g),data=Sonar,control=ctrl) 

## For performance comparison check against agraph with graph only
tab=table(fitted(g.agraph2)[U],y.true[U])
1-sum(diag(tab))/sum(tab)

## Fit agraph to Sonar but graph only 
g.agraph3<-agraph(Class~aG(g),data=Sonar,control=ctrl) 
g.agraph3
tab=table(fitted(g.agraph3)[U],y.true[U])
1-sum(diag(tab))/sum(tab)
}
}
\author{Mark Vere Culp}
\seealso{ \code{\link{agraph}}, \code{\link{aG}}}
\keyword{classes}
\keyword{models}
\keyword{methods}
