################################################################################
##      Programmer: Mark Vere Culp
##      Date:       1-13-2017
##      Purpose:    This file prints all Tables in the paper. Except for
##                  Table 3. Only column 4 is printed. The same code can be run
##                  with the alternative R libraries on different R
##                  configurations to get the different BLAS results. These
##                  Tables are created from the csv files, which can be in-turn
##                  regenerated from the run R-scripts.
##
################################################################################

#####
## Create Performance Means (Table 2)/Time Means (Table version of Figure 2)
####
if(!file.exists("source/data_reader.R")){stop("Must be in the main code directory (refer to README.txt)")}

x<-read.csv("csvs/perf_results.csv",T)
nm<-names(x)
pfmeans<-round(sapply(3:10,function(i)tapply(x[,i],x[,1],mean,na.rm=TRUE)),3)
pfmeans<-as.data.frame(pfmeans)
names(pfmeans)<-nm[-c(1:2)]
ord<-order(pfmeans$n)
pfmeans[ord,]

#####
## Create CPU Time (Linux Ubuntu 8 Xeon E5-2600 v4 processors at 3.5 GHz
##                  using Microsoft R Open 3.2.5 with 8 threads) for
##         Table 3 column 4.
####
x<-read.csv("csvs/casp_CPU.csv",T)
round(apply(x[,-1],2,median),5)

#####
## Create Labeled Means (Table 4)
####

x<-read.csv("csvs/casp_change_m.csv",T)
round(sapply(3:10,function(i)tapply(x[,i],x[,4],median,na.rm=TRUE)),3)

q(save="no")
