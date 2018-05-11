\name{predict.s4pm}
\alias{predict.s4pm}
\alias{predict,s4pm-method}
\title{Out-of-Sample Predict Procedure for \code{s4pm}}
\description{
  This implements the out-of-sample prediction for an
  \sQuote{s4pm} object. 
}
\usage{\S4method{predict}{s4pm}(object,xnew,gnew,type=c("vector","response","prob","terms"),pow=1,\dots)}
\arguments{
  \item{object}{ an existing \sQuote{s4pm} object.}
  \item{xnew}{ an object of class \sQuote{data.frame}, \sQuote{vector}, or \sQuote{matrix}.
  	This is not always necessary and depends on call context (refer to details below).
  	}
  \item{gnew}{ the \sQuote{matrix} of new graph links between the data to predict and 
  	the data used for training. This is not always necessary and depends on call
  	context (refer to details below).}
  \item{type}{ the type of prediction to return. In the case of terms the fits are returned
  	and the individual functions with f1 ss the graph and f2 as the linear term.}
  \item{pow}{ tuning parameter for the weighted power in the
    interpolation predictions.}
  \item{\dots}{mop up additional arguments.}
}
\details{
	The prediction inputs are dependent upon how one calls the original \code{\link{s4pm}} generic function. 
	The cases are discussed next:
	
	1) y~.: This is the default and most common case. Set xnew to your new hold-out data
	       set and do not initialize gnew.
	       
	2) y~dG(G): The gnew argument will [most likely] be a  non-symmetric 
	\sQuote{matrix} of adjacencies between some new set of observations and the original x data.
	
	3) y~sG(G): The gnew argument will [most likely] be a  non-symmetric \sQuote{matrix}
	of similarity adjacencies [most likely] observed directly.

	4) y~.+dG(G) or y~.+sG(G): The xnew data set must be given, but it need not correspond to the data used
	to generate the graph. 	
    
    5) Non-formula call: xnew will either be provided or NULL depending on context, and 
    gnew will have to be provided in this case.
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

## Typical, call to s4pm and predict

g.s4pm<-s4pm(Class~.,data=Sonar[c(L,U1),])
p.s4pm<-predict(g.s4pm,xnew=Sonar[U,-p])
tab=table(y.true[U],p.s4pm)
1-sum(diag(tab))/sum(tab)
\donttest{
## Predict the graph only case
Dij<-x.scaleL(Sonar[,-p],L)
Dij<-as.matrix(cosineDist(Dij))
Dij1<-Dij[c(L,U1),c(L,U1)]
attr(Dij1,"metric")=attr(Dij,"metric")
attr(Dij1,"distance.graph")=attr(Dij,"distance.graph")
g.s4pm<-s4pm(Class~dG(Dij1),data=Sonar[c(L,U1),])
p.s4pm<-predict(g.s4pm,gnew=Dij[U,c(L,U1)])
tab=table(y.true[U],p.s4pm)
1-sum(diag(tab))/sum(tab)
}
}


\keyword{classes}
\keyword{models}
\keyword{methods}
