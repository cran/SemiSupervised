\name{predict.jtharm}
\alias{predict.jtharm}
\alias{predict,jtharm-method}
\title{Out-of-Sample Predict Procedure for \code{jtharm}}
\description{
  This implements the out-of-sample prediction for an
  \sQuote{jtharm} object. 
}
\usage{\S4method{predict}{jtharm}(object,xnew,gnew,type=c("vector","response","prob"),pow=1,\dots)}
\arguments{
  \item{object}{ an existing \sQuote{jtharm} object.}
  \item{xnew}{ an object of class \sQuote{data.frame}, \sQuote{vector}, or \sQuote{matrix}.
  	This is not always necessary and depends on call context (refer to details below).
  	}
  \item{gnew}{ the \sQuote{matrix} of new graph links between the data to predict and 
  	the data used for training. This is not always necessary and depends on call
  	context (refer to details below).}
  \item{type}{ the type of prediction to return.}
  \item{pow}{ tuning parameter for the weighted power in the
    interpolation predictions.}
  \item{\dots}{mop up additional arguments.}
}
\details{
	The prediction inputs are dependent upon how one calls the original \code{\link{jtharm}} generic function. 
	The cases are discussed next:
	
	1) y~.: This is the default and most common case. Set xnew to your new hold-out data
	       set and do not initialize gnew.
	       
	2) y~dG(G): The gnew argument will [most likely] be a  non-symmetric 
	\sQuote{matrix} of adjacencies between some new set of observations and the original x data.
	
	3) y~sG(G): The gnew argument will [most likely] be a  non-symmetric \sQuote{matrix}
	of similarity adjacencies [most likely] observed directly.

    4) Non-formula call: gnew will have to provided in this case but xnew is ignored.
}
\value{
	If \code{type}(object) is \sQuote{r}, a vector of predicted values is
     returned. If \code{type}(object) is \sQuote{c}, the object returned depends
     on the type argument.
}

\author{Mark Vere Culp}

\examples{
## Prediction depends on the nature of the call. Consider some examples.
library(mlbench)
data(Sonar)

n=dim(Sonar)[1]
p=dim(Sonar)[2]

nu=0.2
set.seed(100)
L=sort(sample(1:n,ceiling(nu*n)))
U=setdiff(1:n,L)
U1=sample(U,ceiling(0.5*n))

y.true<-Sonar$Class
Sonar$Class[U]=NA

## Typical, call to jtharm and predict
g.jtharm1<-jtharm(Class~.,data=Sonar[c(L,U1),])
p.jtharm1<-predict(g.jtharm1,xnew=Sonar[U,-p])
tab=table(y.true[U],p.jtharm1)
1-sum(diag(tab))/sum(tab)
\donttest{
## Predict the graph only case Debug later
Dij<-x.scaleL(Sonar[,-p],L)
Dij<-as.matrix(cosineDist(Dij))
Dij1<-Dij[c(L,U1),c(L,U1)]
attr(Dij1,"metric")=attr(Dij,"metric")
attr(Dij1,"distance.graph")=attr(Dij,"distance.graph")
g.jtharm2<-jtharm(Class~dG(Dij1),data=Sonar[c(L,U1),])
p.jtharm2<-predict(g.jtharm2,gnew=Dij[U,c(L,U1)])
tab=table(y.true[U],p.jtharm2)
1-sum(diag(tab))/sum(tab)
}
}


\keyword{classes}
\keyword{models}
\keyword{methods}
