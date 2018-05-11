\name{agraph}
\alias{agraph}
\alias{agraph-methods}
\alias{agraph,formula-method}
\alias{agraph,data.frame-method}
\alias{agraph,matrix-method}
\alias{agraph,vector-method}
\alias{agraph,anchor-method}
\alias{show,agraph-method}
\alias{agraph.default}
\title{Anchor Graph Functions (\code{agraph})}
\description{
  This performs the \code{agraph} machine learning method in \R as
  described in the references below. This is documented in the user manual and references.
}
\usage{
\S4method{agraph}{formula}(x,data,metric= c("cosine","euclidean"),...,
  est.only=FALSE,control=SemiSupervised.control())

\S4method{agraph}{matrix}(x,y,\dots,metric="cosine",est.only=FALSE, control=SemiSupervised.control())

\S4method{agraph}{vector}(x,\dots)

\S4method{agraph}{data.frame}(x,\dots)

\S4method{agraph}{anchor}(x,\dots)


## Default Method (not meant to be run direclty).
agraph.default(graph,x,y,weights,lams,gams,type=c("r","c"),
               est.only=FALSE,control=SemiSupervised.control())
}
\arguments{
  \item{x}{a symbolic description of the model to be fit. A formula can
  	be used directly or the feature data/graph are inputted
  	directly. Refer to the details below for more on this.} 
  \item{data}{the response and or feature data are a \sQuote{data.frame}. The 
  	missing or unlabeled responses must be NA.}
  \item{graph}{an \sQuote{anchor} object that resulted from a prior \code{\link{AnchorGraph}} call. 
  	This is only required for the default method.}
   \item{metric}{the metric to be fit. It is either Cosine dissimilarity
	or the Euclidean distance (using the \code{\link{daisy}} function
	in the \pkg{cluster} library. Only used with the \sQuote{formula} call.
	Default is \dQuote{cosine}.}
   \item{y}{response could be either of length m=|L|<=n with the first
	m observations labeled, or be of length n with the unlabeled
	cases flagged as NA. Refer to details below.}
  \item{weights}{an optional vector of weights to be used in the fitting
          process. If missing the unweighted algorithm is fit.}
  \item{lams}{a \sQuote{vector} or \sQuote{matrix} of Lagrangian parameters for the graph 
	and ridge penalties depending on the call. If missing then the
	procedure provides the grid by default.}
\item{gams}{a vector of Lagrangian parameters corresponding to the 
	latent unlabeled response penalty. If missing then the
	procedure provides the grid by default.}	
 \item{type}{use \sQuote{r} for regression, which is the default
	unless y is a factor. The \sQuote{c} option performs
	classification with logistic loss.}
\item{est.only}{returns only the fitted vector (sign is the class in
	classification) with no agraph object. Designed for quickly
	fitting.}
\item{control}{ control parameters (refer to \code{\link{SemiSupervised.control}} for
    more information).}
\item{\dots}{mop up additional inputs and checks against the \code{agraph.default}
  	arguments.}
}
\details{
Details on Response: 
  
  The response is set where NA's denote the missing responses (unlabeled set), and 
  values (either \sQuote{factor} or \sQuote{numeric}) for the known cases (labeled set). 
  This type of input must be used for the formula interface. In some circumstances, 
  the response can be inputted with the first m cases defined. So the
  dimension of argument y and graph/x are different. 
  
Details on Formula:
  
  1) y~.: The most common case. One starts with a typical data set x and 
    wishes to fit the model xb+f(G[x]) but approximates f(G[x]) with linear term
    Za using anchors. The predict function and all other aspects are done internally
    just like a typical \R function, i.e., the algorithm inputs are the same as if a linear 
    model or random forest were fit. 
   
  2) y~.+aG(g): The graph g corresponds an object resulting from \code{\link{AnchorGraph}}. 
     Refer to \code{\link{aG}} for some examples. An additional feature data set whose 
     is also provided. This could be the data set used to construct g (safe mode) or an 
     additional data set. 
   
  3) y~aG(g): This bypasses the safe semi-parametric component of 
     the model. A more traditional semi-supervised graph-only method is fit. This tends 
     to perform worse but is faster to fit. 
  
  4) Non-formula call: The \sQuote{formula} is a much simpler interface, but the inputs
     of the function can be manipulated through the \sQuote{data.frame}, \sQuote{vector},
     \sQuote{matrix}, and \sQuote{NULL} interfaces of \code{agraph}. The NULL interface fits 
     version (3). The \code{agraph.default} could also be fit, but it is not recommended. 

  Other details:
  
  The \code{predict} generic for \code{agraph} is optimized for formula
  input. This is the simplest  way to input the function and build the
  graph. However, one can use these functions with direct inputs, but 
  later functions also need to adjust accordingly. 
  
  The code provides semi-supervised graph-based support for \R.
}
\value{
	An object of class \code{agraph}.
}
\note{
  K-fold cross-validation has now been implement in C++ to perform parameter
  estimation for the graph penalty, the safe ridge penalty, and 
  the latent unlabeled response penalty.
}
\author{Mark Vere Culp}

\references{
MV Culp and KJ Ryan (2016). SemiSupervised: Scalable Semi-Supervised Routines for Real Data Problems.
}


\examples{
## Set up Sonar data with 20% labeled (comparing randomForest and glmnet)
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

## Fit agraph to Sonar
g.agraph<-agraph(Class~.,data=Sonar) 
g.agraph
tab=table(fitted(g.agraph)[U],y.true[U])
1-sum(diag(tab))/sum(tab)

}
\keyword{classes}
\keyword{models}
\keyword{methods}
