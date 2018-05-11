library("SemiSupervised")
library("caret")
data("BloodBrain")
n<-nrow(bbbDescr)
set.seed(100)
L<-sample(1:n, ceiling(0.1 * n))
U<-setdiff(1:n, L)
y.fit<-rep(NA, n)
y.fit[L]<-logBBB[L]
msy<-sqrt(mean((logBBB[L] - mean(logBBB[L]))^2))
sRMSEU<-function(f, y, msy = 1.0){ sqrt(mean((f - y)^2)) / msy }

## Safe Combinatorial Laplacian Example

safe.spread<-s4pm(y ~ ., data = data.frame(y = y.fit, bbbDescr))
safe.spread
sRMSEU(fitted(safe.spread)[U], logBBB[U], msy)

## Safe Normalized Laplacian Example

ctrl<-SemiSupervised.control(normalize = FALSE)
safe.lap<-s4pm(y ~ ., data = data.frame(y = y.fit, bbbDescr), control = ctrl)
safe.lap
sRMSEU(fitted(safe.lap)[U], logBBB[U], msy)

##Regularized Laplacians
xscale<-x.scaleL(bbbDescr, L)
rlap<-s4pm(y.fit ~ dG(xscale, metric = "cosine"), gam = 0.0,
           control = ctrl)
rlap
sRMSEU(fitted(rlap)[U], logBBB[U], msy)

##Diffusion Kernel
diffkern<-s4pm(y.fit ~ dG(xscale, metric = "cosine"), gam = 0.0)
sRMSEU(fitted(diffkern)[U], logBBB[U], msy)

##Label Spreading
lspread<-s4pm(y.fit ~ dG(xscale, metric = "cosine"), gam = Inf)
sRMSEU(fitted(lspread)[U], logBBB[U], msy)


## Variable Selection with s4pm (Takes 20 seconds or so)

drp1<-(do.call("c", lapply(names(bbbDescr), function(i){
  form<-as.formula(paste("y ~ . - ", i, sep = ""))
   g1<-s4pm(form,data=data.frame(y = y.fit, bbbDescr),
             hs=hparm(safe.spread), gams=gparm(safe.spread),
             lams=lparm(safe.spread))
  measures(g1)[2]
})) / measures(safe.spread)[2] - 1) * 100
attr(drp1,"names")<-names(bbbDescr)
round(sort(drp1, decreasing=TRUE)[1:5], 4)

####
## Anchor Graph Based
####
ctrl.agr <- SemiSupervised.control(cn = 50)
anchors <- rbind(c(1, 2), c(2, 5), c(3, 4), c(7, 2))
(z <- AnchorGraph(c(2, 3), anchor = anchors, control = ctrl.agr)$Z)
z %*% anchors

z <- AnchorGraph(c(2.5, 1), anchor = anchors, control = ctrl.agr)$Z
as.vector(z)

safe.agraph <- agraph(y ~ ., data = data.frame(y = y.fit, bbbDescr))
safe.agraph
sRMSEU(fitted(safe.agraph)[U], logBBB[U], msy)

agr <- agraph(y.fit ~ aG(AnchorGraph(xscale)))
sRMSEU(fitted(agr)[U], logBBB[U], msy)

#####
## Out-of-Sample Predictions
#####
set.seed(10)
xnew=as.data.frame(jitter(as.matrix(bbbDescr)[1:5, ]))
round(predict(safe.lap, xnew = xnew),3)
round(predict(safe.spread, xnew = xnew),3)
round(predict(safe.agraph, xnew = xnew),3)


sxnew=scale(xnew,center = attr(xscale, "scaled:center"),
                 scale = attr(xscale, "scaled:scale"))
gnew=kgraph.predict(cosineDist(sxnew, xscale))
round(predict(rlap, gnew = gnew), 3)
round(predict(agr, gnew = AnchorGraph(sxnew, fit.g = gmatrix(agr))), 3)

#####
### Observe Graph
#####
onePos<-list(1:2, 1:4, c(2:3, 6), c(2, 4, 5), c(4, 5), c(3, 6))
(W<-do.call("rbind", lapply(onePos, function(i){v = rep(0, 6);v[i] = 1;v})))
y<-as.factor(c("a", "b", NA, NA, "a", "b"))
graph.ex<-s4pm(y ~ sG(W))
graph.ex
fitted(graph.ex)
round(fitted(graph.ex, type = "prob"), 2)


########
### Joint Harmonic Functions
#######
safe.jtharm<-jtharm(y ~ sG(W))
safe.jtharm

########
## The spa
########
library("spa")
safe.spa<-spa(y ~ ., data=data.frame(y = y.fit, bbbDescr))
sRMSEU(fitted(safe.spa)[U], logBBB[U], msy)


