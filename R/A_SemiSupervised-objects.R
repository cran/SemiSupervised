################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: S4 Generic Classes -- The structure of this file
##               was based off of the kernlab library.
##
################################################################################


## The next 4 functions union basic R objects for convience in object
## definitions.

setClassUnion("input", c("matrix","list"))
setClassUnion("listI", c("list","numeric","vector","integer","matrix"))
setClassUnion("output", c("matrix","factor","vector","logical","numeric","list","integer","NULL"))
setClassUnion("parm", c("vector","numeric","integer","NULL"))

## An object of type SemiSupervised. This is a virtual class which is used as 
## the base functionality for the semi-supervised objects.

setClass("SemiSupervised", representation(
  .respinfo="list",
  .cv_str = "output",
  .fitinfo = "output",
  .terminfo = "output",
  .call="call",
  .control = "listI",
  type = "character",
  gmatrix = "ANY",
  ymatrix = "output",
  xmatrix = "output",
  measures = "output",
  fit = "output"),
  contains= "VIRTUAL")

setClass("jtharm", representation(
lparm = "parm",
hparm = "parm",
gparm = "parm"),contains ="SemiSupervised")

setClass("s4pm", representation(
lparm = "parm",
hparm = "parm",
gparm = "parm"),contains ="SemiSupervised")

setClass("agraph", representation(
lparm = "parm",
gparm = "parm"),contains ="SemiSupervised")



## Graph classes are defined below

setClass("lapGraph",representation(graph="listI"),contains="list")
setClass("anchor",representation(graph="listI"),contains="list")

## The routine generic functions are defined next. These functions
## are routine accessor functions to the internal class attributes.

## Generic for virtual class SemiSupervised.

if(!isGeneric("type")){
  if (is.function("type"))
  fun <- type
  else fun <- function(object) standardGeneric("type")
  setGeneric("type", fun)
}
setMethod("type", "SemiSupervised", function(object) object@type)
setGeneric("type<-", function(x, value) standardGeneric("type<-"))
setReplaceMethod("type", "SemiSupervised", function(x, value) {
  x@type <- value
  x
})

setMethod("dim","SemiSupervised",function(x){
  if(!is.null(x@.fitinfo)){
    slot(x,".fitinfo")$dims[-3]
  }
})


setMethod("fitted", "SemiSupervised", function(object, type=c("vector","response","prob"),...){
  if(!is.null(object@.fitinfo)){
    f=object@fit
    ctype=slot(object,"type")
    if(ctype=="r"){
      return(f)
    }
    if(missing(type)){type="response"}
    if(type=="vector"){
      return(f)
    }
    if(type=="prob"){
      return(exp(f)/(1+exp(f)))
    }
    lev=slot(object,".respinfo")$lev
    return(factor(lev[as.numeric(object@fit>0.0)+1],levels=lev))
  }
})

if(!isGeneric("fit")){
  if (is.function("fit"))
  fun <- xmatrix
  else fun <- function(object) standardGeneric("fit")
  setGeneric("fit", fun)
}
setGeneric("fit<-", function(x, value) standardGeneric("fit<-"))
setReplaceMethod("fit", "SemiSupervised", function(x, value) {
  x@fit <- value
  x
})

if(!isGeneric("xmatrix")){
  if (is.function("xmatrix"))
  fun <- xmatrix
  else fun <- function(object) standardGeneric("xmatrix")
  setGeneric("xmatrix", fun)
}
setMethod("xmatrix", "SemiSupervised", function(object) object@xmatrix)
setGeneric("xmatrix<-", function(x, value) standardGeneric("xmatrix<-"))
setReplaceMethod("xmatrix", "SemiSupervised", function(x, value) {
  x@xmatrix <- value
  x
})

if(!isGeneric("ymatrix")){
  if (is.function("ymatrix"))
  fun <- ymatrix
  else fun <- function(object) standardGeneric("ymatrix")
  setGeneric("ymatrix", fun)
}
setMethod("ymatrix", "SemiSupervised", function(object) object@ymatrix)
setGeneric("ymatrix<-", function(x, value) standardGeneric("ymatrix<-"))
setReplaceMethod("ymatrix", "SemiSupervised", function(x, value) {
  x@ymatrix <- value
  x
})

if(!isGeneric("gmatrix")){
  if (is.function("gmatrix"))
  fun <- gmatrix
  else fun <- function(object) standardGeneric("gmatrix")
  setGeneric("gmatrix", fun)
}
setMethod("gmatrix", "SemiSupervised", function(object) object@gmatrix)
setGeneric("gmatrix<-", function(x, value) standardGeneric("gmatrix<-"))
setReplaceMethod("gmatrix", "SemiSupervised", function(x, value) {
  x@gmatrix <- value
  x
})

if(!isGeneric("measures")){
  if (is.function("measures"))
  fun <- zmatrix
  else fun <- function(object) standardGeneric("measures")
  setGeneric("measures", fun)
}
setMethod("measures", "SemiSupervised", function(object) object@measures)
setGeneric("measures<-", function(x, value) standardGeneric("measures<-"))
setReplaceMethod("measures", "SemiSupervised", function(x, value) {
  x@measures <- value
  x
})

## Generic functions for class s4pm 

if(!isGeneric("lparm")){
  if (is.function("lparm"))
  fun <- lparm
  else fun <- function(object) standardGeneric("lparm")
  setGeneric("lparm", fun)
}
setMethod("lparm", "s4pm", function(object) object@lparm)
setGeneric("lparm<-", function(x, value) standardGeneric("lparm<-"))
setReplaceMethod("lparm", "s4pm", function(x, value) {
  x@lparm <- value
  x
})

if(!isGeneric("hparm")){
  if (is.function("hparm"))
  fun <- hparm
  else fun <- function(object) standardGeneric("hparm")
  setGeneric("hparm", fun)
}
setMethod("hparm", "s4pm", function(object) object@hparm)
setGeneric("hparm<-", function(x, value) standardGeneric("hparm<-"))
setReplaceMethod("hparm", "s4pm", function(x, value) {
  x@hparm <- value
  x
})
if(!isGeneric("gparm")){
  if (is.function("gparm"))
  fun <- gparm
  else fun <- function(object) standardGeneric("gparm")
  setGeneric("gparm", fun)
}
setMethod("gparm", "s4pm", function(object) object@gparm)
setGeneric("gparm<-", function(x, value) standardGeneric("gparm<-"))
setReplaceMethod("gparm", "s4pm", function(x, value) {
  x@gparm <- value
  x
})

## Generic functions for class jtharm

if(!isGeneric("lparm")){
  if (is.function("lparm"))
  fun <- lparm
  else fun <- function(object) standardGeneric("lparm")
  setGeneric("lparm", fun)
}
setMethod("lparm", "jtharm", function(object) object@lparm)
setGeneric("lparm<-", function(x, value) standardGeneric("lparm<-"))
setReplaceMethod("lparm", "jtharm", function(x, value) {
  x@lparm <- value
  x
})

if(!isGeneric("hparm")){
  if (is.function("hparm"))
  fun <- hparm
  else fun <- function(object) standardGeneric("hparm")
  setGeneric("hparm", fun)
}
setMethod("hparm", "jtharm", function(object) object@hparm)
setGeneric("hparm<-", function(x, value) standardGeneric("hparm<-"))
setReplaceMethod("hparm", "jtharm", function(x, value) {
  x@hparm <- value
  x
})
if(!isGeneric("gparm")){
  if (is.function("gparm"))
  fun <- gparm
  else fun <- function(object) standardGeneric("gparm")
  setGeneric("gparm", fun)
}
setMethod("gparm", "jtharm", function(object) object@gparm)
setGeneric("gparm<-", function(x, value) standardGeneric("gparm<-"))
setReplaceMethod("gparm", "jtharm", function(x, value) {
  x@gparm <- value
  x
})

## Generic functions for class agraph

if(!isGeneric("gparm")){
  if (is.function("gparm"))
  fun <- gparm
  else fun <- function(object) standardGeneric("gparm")
  setGeneric("gparm", fun)
}
setMethod("gparm", "agraph", function(object) object@gparm)
setGeneric("gparm<-", function(x, value) standardGeneric("gparm<-"))
setReplaceMethod("gparm", "agraph", function(x, value) {
  x@gparm <- value
  x
})

if(!isGeneric("lparm")){
  if (is.function("lparm"))
  fun <- lparm
  else fun <- function(object) standardGeneric("lparm")
  setGeneric("lparm", fun)
}
setMethod("lparm", "agraph", function(object) object@lparm)
setGeneric("lparm<-", function(x, value) standardGeneric("lparm<-"))
setReplaceMethod("lparm", "agraph", function(x, value) {
  x@lparm <- value
  x
})

if(!isGeneric("lparm")){
  if (is.function("lparm"))
  fun <- lparm
  else fun <- function(object) standardGeneric("lparm")
  setGeneric("lparm", fun)
}
setMethod("lparm", "agraph", function(object) object@lparm)
setGeneric("lparm<-", function(x, value) standardGeneric("lparm<-"))
setReplaceMethod("lparm", "agraph", function(x, value) {
  x@lparm <- value
  x
})
