\name{predict.agraph}
\alias{predict.agraph}
\alias{predict,agraph-method}
\title{Out-of-Sample Predict Procedure for \code{agraph}}
\description{
  This implements the out-of-sample prediction for an
  \sQuote{agraph} object. 
}
\usage{\S4method{predict}{agraph}(object,xnew,gnew,type=c("vector","response","prob"),\dots)}
\arguments{
  \item{object}{ an existing \sQuote{agraph} object.}
  \item{xnew}{ an object of class \sQuote{data.frame}, \sQuote{vector}, or \sQuote{matrix}.
  	This is not always necessary and depends on call context (refer to details below).
  	}
  \item{gnew}{ the \sQuote{matrix} of new graph links between the data to predict and 
  	the data used for training. This is not always necessary and depends on call
  	context (refer to details below).}
  \item{type}{ the type of prediction to return.}
  \item{\dots}{mop up additional arguments.}
}
\details{
	The prediction inputs are dependent upon how one calls the original \code{\link{agraph}} generic function. 
	The cases are discussed next:
	
	1) y~.: This is the default and most common case. Set xnew to your new hold-out data
	       set and do not initialize gnew.
	       
	2) y~aG(g): Prior to applying the predict, we must first update the \sQuote{anchor} object used
	in the original \code{\link{agraph}} call to the xnew data. To do this, invoke \code{\link{AnchorGraph}} 
	call with fit.g set to the original \sQuote{anchor} object. Fit the updated Anchor Graph as
	the gnew parameter and ignore the xnew parameter.
	
	3)  y~. + aG(g): Prior to applying the predict, we must first update the \sQuote{anchor} object used
	in the original \code{\link{agraph}} call to some new data set that correponds to rows of xnew (it 
	techically doesn't have to be xnew'). Proceed as in case (2) except the xnew data must be 
	provided.
	
    4) Non-formula call: xnew will either be provided or NULL depending on context, and 
    gnew will have to be provided as in case two.
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

## Typical, call to agraph and predict

g.agraph1<-agraph(Class~.,data=Sonar[c(L,U1),])
p.agraph1<-predict(g.agraph1,xnew=Sonar[U,-p])
tab=table(y.true[U],p.agraph1)
1-sum(diag(tab))/sum(tab)
\donttest{
## Predict the graph only case
ctrl=SemiSupervised.control()
scale.x=x.scaleL(Sonar[,-p],L)
g=AnchorGraph(scale.x[c(L,U1),-p],control=ctrl)
g.agraph2<-agraph(Class~aG(g),data=Sonar[c(L,U1),],control=ctrl)
g<-AnchorGraph(scale.x[U,-p],fit.g=g,control=ctrl)
p.agraph2<-predict(g.agraph2,gnew=g)
tab=table(y.true[U],p.agraph2)
1-sum(diag(tab))/sum(tab)
}
}


\keyword{classes}
\keyword{models}
\keyword{methods}
