\name{jtharm}
\alias{jtharm}
\alias{jtharm-methods}
\alias{jtharm,formula-method}
\alias{jtharm,data.frame-method}
\alias{jtharm,matrix-method}
\alias{show,jtharm-method}
\alias{jtharm.default}
\title{Joint Harmonic Functions (\code{jtharm})}
\description{
This performs the \code{jtharm} machine learning method in \R as
described in the references below. It can fit a graph only estimate
(y=f(G)). No linear term is used since the safe mechanism works
through regularization of the latent response as described in the 
references below.
}
\usage{
\S4method{jtharm}{formula}(x,data,metric= c("cosine","euclidean"),...,
  est.only=FALSE,control=SemiSupervised.control())

\S4method{jtharm}{matrix}(x,\dots)

## Default Method (not meant to be run direclty).
jtharm.default(graph,y, weights,hs,lams,gams,type=c("r","c"),est.only=FALSE,
             control=SemiSupervised.control())
}
\arguments{
\item{x}{a symbolic description of the model to be fit. A \sQuote{formula} can
  	be used directly or the feature data/graph are inputted
  	directly. Refer to the details below for more on this.} 
\item{data}{the response and or feature data are a \sQuote{data.frame}. The 
  	missing or unlabeled responses must be NA.}
\item{graph}{the graph matrix. Typically an n by n dissimilarity matrix but can
  	also be similarity. Refer to the documentation in the details section 
  	below.}
\item{y}{response could be either of length m=|L|<=n with the first
	m observations labeled or be of length n with the unlabeled
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
\item{lams}{a vector or matrix of Lagrangian parameters for the graph penalty 
	term. If missing then the procedure provides the grid by default.}
\item{gams}{a vector of Lagrangian parameters corresponding to the 
	latent unlabeled response penalty. If missing then the
	procedure provides the grid by default.}	
\item{type}{use \dQuote{r} for regression, which is the default
	unless y is a factor. The \dQuote{c} option performs
	classification.} 
\item{est.only}{returns only the fitted vector (sign is the class in
	classification) with no \sQuote{jtharm} object. Designed for quickly
	fitting.}
\item{control}{control parameters (refer to \code{\link{SemiSupervised.control}} for
    more information).}
  \item{\dots}{mop up additional inputs and checks against the \code{jtharm.default}
  	arguments.}
}
\details{
  Details on Response: 
  
  The response is set where NA's denote the missing responses (unlabeled set) and 
  values (either \sQuote{factor} or \sQuote{numeric}) for the known cases (labeled set). 
  This type of input must be used for the \sQuote{formula} interface. In some circumstances, 
  the response can be inputted with the first m cases defined. So the
  dimension of argument y and graph/x are different in this case. 
  
  Details on Formula:
  
  1) y~.: The most common case. One starts with a typical data set x and 
    wishes to fit the model f(G[x]). All fitting is done internally using a
    k=6 NN graph and cosine distance as default but these can be modified 
    respectively through the \code{\link{SemiSupervised.control}} or the metric argument. 
    The predict function and all other aspects are done internally just like a 
    typical \R function, i.e., the algorithm inputs are the same as if a linear 
    model or random forest were fit. 
   
  2) y~dG(G[x]) or y~sG(G): This fits the semi-supervised graph-only method using a graph fit
     previously outside of the function. The \code{\link{sG}} corresponds to a similarity
     (\sQuote{0} is dissimilar and \sQuote{1} [typically] is similar) graph while \code{\link{dG}}
     corresponds to a dissimilarity matrix (\sQuote{0} is similar, \sQuote{Inf} is dissimilar)
  
  3) Non-formula call: The \sQuote{formula} is a much simpler interface but the inputs
     of the function can be manipulated through the sQuote{matrix} interfaces of \code{jtharm}. 
     The NULL interface fits version (4) of the \sQuote{formula } call above. The 
     \code{jtharm.default} could also be fit, but it is not recommended. 

  
  Other details:
  
  The approach only fits the Laplace kernel for the Graph. Practically, 
  this is all that is necessary since we optimize the tuning parameters  
  for it. Also, the interest is in a computationally fast algorithm, so 
  determining the grid internally is optimized for performance and time.

  The \code{predict} generic for \code{jtharm} is optimized for the \sQuote{formula}
  input. This is the simplest way to input the function and build the
  graph. 
  
  The code provides semi-supervised graph-based support for \R.
}
\value{
	An object of class \code{jtharm}. 
}
\note{
  K-fold cross-validation was implement in C++ to perform parameter
  estimation for the h, graph penalty, and 
  the latent unlabeled response penalties. Several LAPACK routines are used to fit 
  the underlying functions for both CV and regular fit. 
}
\author{Mark Vere Culp}

\references{
MV Culp and KJ Ryan (2013). On Joint Harmonic Functions and Their Supervised Connections. Journal of Machine Learning Research. 14:3721--3752.

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

##Fit jtharm to Sonar
g.jtharm<-jtharm(Class~.,data=Sonar)
tab=table(fitted(g.jtharm)[U],y.true[U])
1-sum(diag(tab))/sum(tab)
}
\keyword{classes}
\keyword{models}
\keyword{methods}
