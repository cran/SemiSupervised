\name{s4pm}
\alias{s4pm}
\alias{s4pm-methods}
\alias{s4pm,formula-method}
\alias{s4pm,NULL-method}
\alias{s4pm,lapGraph-method}
\alias{s4pm,vector-method}
\alias{s4pm,data.frame-method}
\alias{s4pm,matrix-method}
\alias{show,s4pm-method}
\alias{s4pm.default}
\title{Safe Semi-Supervised Semi-Parametric Model (\code{s4pm})}
\description{
	This performs the \code{s4pm} machine learning method in \R as
  described in the references below. It can fit a graph only estimate
  (y=f(G)) and graph-based semi-parametric estimate (y=xb+f(G)). Refer
  to the example below.
  
  This is documented in the user manual and references.
 }
\usage{
\S4method{s4pm}{formula}(x,data,metric= c("cosine","euclidean"),...,
  est.only=FALSE,control=SemiSupervised.control())

\S4method{s4pm}{matrix}(x,y,graph,\dots)

\S4method{s4pm}{data.frame}(x,y,graph,\dots)

\S4method{s4pm}{vector}(x,y,graph,\dots)

\S4method{s4pm}{NULL}(x,y,graph,\dots)

## Default Method (not meant to be run direclty).
s4pm.default(x,y,graph,weights,hs,lams,gams,type=c("r","c"), est.only=FALSE,
	control=SemiSupervised.control())
}
\arguments{
  \item{x}{a symbolic description of the model to be fit. A formula can
  	be used directly or the feature data/graph are inputted
  	directly. Refer to the details below for more on this.} 
  \item{data}{the response and or feature data are a \sQuote{data.frame}. The 
  	missing or unlabeled responses must be NA.}
   \item{graph}{the graph matrix. Typically an n by n dissimilarity matrix but can
  	also be similarity. Refer to the documentation in the details section 
  	below.}
   \item{y}{response could be either of length m=|L|<=n with the first
	m observations labeled, or be of length n with the unlabeled
	cases flagged as NA. Refer to details below.}
  \item{weights}{an optional vector of weights to be used in the fitting
          process. If missing the unweighted algorithm is fit.}
  \item{metric}{the metric to be fit. It is either Cosine dissimilarity
	or the Euclidean distance (using the \code{\link{daisy}} function
	in the \pkg{cluster} library). Only used with the \sQuote{formula} call.
	Default is \dQuote{cosine}.}
  \item{hs}{a sequence of tuning parameters for the kernel function if
	the graph is a dissimilarity matrix. If missing then the
	procedure provides the grid by default.}
  \item{lams}{a \sQuote{vector} or \sQuote{matrix} of Lagrangian parameters for the graph 
	and ridge penalties depending on the call. If missing then the
	procedure provides the grid by default.}
\item{gams}{a vector of Lagrangian parameters corresponding to the 
	latent unlabeled response penalty. If missing then the
	procedure provides the grid by default.}	
    \item{type}{use \dQuote{r} for regression, which is the default
	unless y is a factor. The \dQuote{c} option performs
	classification with logistic loss.}
  \item{est.only}{returns only the fitted vector (sign is the class in
	classification) with no \sQuote{s4pm} object. Designed for quickly
	fitting.}
  \item{control}{control parameters (refer to \code{\link{SemiSupervised.control}} for
    more information).}
  \item{\dots}{mop up additional inputs and checks against the \code{s4pm.default}
  	arguments.}
}
\details{
  Details on Response: 
  
  The response is set where NA's denote the missing responses (unlabeled set) and 
  values (either \sQuote{factor} or \sQuote{numeric}) for the known cases (labeled set). 
  This type of input must be used for the \sQuote{formula} interface. In some circumstances, 
  the response can be inputted with the first m cases defined. So the
  dimension of argument y and graph/x are different. 
  
  Details on Formula:
  
  1) y~.: The most common case. One starts with a typical data set x and 
    wishes to fit the model xb+f(G[x]). All fitting is done internally using a
    k=6 NN graph and cosine distance as default but these can be modified 
    respectively through the \code{\link{SemiSupervised.control}} or the metric argument. 
    The predict function and all other aspects are done internally just like a 
    typical \R function, i.e., the algorithm inputs are the same as if a linear 
    model or random forest were fit. 
   
   2) y~.+dG(G[x]): The graph corresponds to a dissimilarity matrix 
   (\sQuote{0} is similar, \sQuote{Inf} is dissimilar). The graph is 
   constructed outside the function using the \code{\link{knnGraph}} command, 
   and the data argument must be set to the feature data part with one column 
   for the response. This is convenient for benchmarking. 
   
   Note the dG(G[x],k=5L) allows one to modify k and 
   \code{\link{dG}}(G[x],nok=TRUE) bypasses the \code{\link{knnGraph}} treating
   the matrix G[x] as a proper adjacency matrix. Refer to \code{\link{dG}} 
   for some examples.
   
   3) y~.+sG(G): This case corresponds to a similarity (\sQuote{0} is dissimilar and
   \sQuote{1} [typically] is similar) graph and an
      additional feature data set whose rows correspond to the nodes of G. Refer 
      to \code{\link{dG}} help page for some examples.
      
   4) y~dG(G[x]) or y~sG(G): This bypasses the safe semi-parametric component of 
    the model. A more traditional semi-supervised graph-only method is fit. This tends 
      to perform worse but is faster to fit. 
  
  5) Non-formula call: The \sQuote{formula} is a much simpler interface, but the inputs
     of the function can be manipulated through the \sQuote{data.frame}, \sQuote{vector},
     \sQuote{matrix}, and \sQuote{NULL} interfaces of \code{s4pm}. The NULL interface 
     fits version (4) of the \sQuote{formula } call above. The \code{s4pm.default} could 
     also be fit, but it is not recommended. 
  
  Other Details:
  
  The approach only fits the Laplace kernel for the graph. Practically, 
  this is all that is necessary since we optimize the tuning parameters  
  for it. Also, the interest is in a computationally fast algorithm, so 
  determining the grid internally is optimized for performance and time.

  The \code{predict} generic for \code{s4pm} is optimized for the \sQuote{formula}
  input. This is the simplest way to input the function and build the
  graph. 
  
  The code provides semi-supervised graph-based support for \R.
}
\value{
	An object of class \code{\link{s4pm-class}}. 
}
\note{
  K-fold cross-validation was implement in C++ to perform parameter
  estimation for the h, graph penalty, safe ridge, and 
  the latent unlabeled response penalties. Several LAPACK routines are used to fit 
  the underlying functions for both CV and regular fit. 
}
\author{Mark Vere Culp}

\references{
MV Culp, KJ Ryan, and P Banerjee (2015). On Safe Semi-supervised Learning. IEEE Pattern Analysis and Machine Intelligence. Submitted.
  
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

## Fit s4pm to Sonar
g.s4pm<-s4pm(Class~.,data=Sonar) 
g.s4pm
tab=table(fitted(g.s4pm)[U],y.true[U])
1-sum(diag(tab))/sum(tab)

## For comparison
\donttest{
library(caret)
library(randomForest)
     
g.glmnet=train(Class~.,data=Sonar[L,],method="glmnet",preProc = c("center", "scale"))
tab=table(predict(g.glmnet,newdata=Sonar[U,-p]),y.true[U])
1-sum(diag(tab))/sum(tab)
     
g.rf<-randomForest(Class~.,data=Sonar[L,])
tab=table(predict(g.rf,newdata=Sonar[U,-p]),y.true[U])
1-sum(diag(tab))/sum(tab)
}
}
\keyword{classes}
\keyword{models}
\keyword{methods}
