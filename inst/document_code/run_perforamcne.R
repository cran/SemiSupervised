####
## Load Libraries
####

libloader<-function(txt){
  g<-try(library(package=txt,character.only=TRUE,quietly=TRUE),silent=TRUE)
  if(class(g)=="try-error"){
    stop(paste("Must install package: ",txt," in order to proceed.",sep=""))
  }
}
g<-lapply(c("SemiSupervised","glmnet","mlbench","ada","bigRR"),libloader)
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
nruns<-100
iter<-1
nm<-c("data","iter","n","m","p","s4pm","agraph","glmnet","s4pm_time","agraph_time")
final<-matrix(0,nruns*10,length(nm))
final<-as.data.frame(final)
names(final)<-nm
filenm<-"csvs/perf_results.csv"
lalp<-7
alpha<-seq(0.0,1.0,length=lalp)

####
## Execute Performance
####
t1<-proc.time()
for(i12 in 1:10){
  dat<- read.data(i12)
  x<-dat$x
  y<-dat$y
  n<-dat$n
  p<-dat$p
  nms<-dat$nm
  type<-strsplit(dat$type,"")[[1]][1]
  if(n<1001){
    m<-floor(0.15*n)
   }else{
    m<-200
  }
   source("source/base.R")
}

q(save="no")
