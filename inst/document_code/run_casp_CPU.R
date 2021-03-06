####
## Load Libraries
####

libloader<-function(txt){
  g<-try(library(package=txt,character.only=TRUE,quietly=TRUE),silent=TRUE)
  if(class(g)=="try-error"){
    stop(paste("Must install package: ",txt," in order to proceed.",sep=""))
  }
}
g<-lapply(c("SemiSupervised","glmnet"),libloader)
if(!file.exists("source/data_reader.R")){stop("Must be in the main code directory (refer to README.txt)")}
source("source/data_reader.R")

###
## Run/System Information
###
read.data()
sessionInfo()

####
## Setup Run
####
nruns<-10

iter<-1
nm<-c("data","iter","n","m","p","s4pm","agraph","glmnet","s4pm_time","agraph_time")
final<-matrix(0,nruns,length(nm))
final<-as.data.frame(final)
names(final)<-nm
ctrl.agraph<-agraph.control()
ctrl.s4pm1<-s4pm.control(normalize=FALSE)
lalp<-7
alpha<-seq(0.0,1.0,length=lalp)
filenm<-"csvs/casp_CPU.csv"
####
## Setup Run
####
dat<-read.data(8)
x<-dat$x
y<-dat$y
n<-dat$n
p<-dat$p
m=200
type<-"r"
nms<-dat$nm
t1<-proc.time()
source("source/base.R")
q(save="no")
