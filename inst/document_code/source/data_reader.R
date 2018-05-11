library(SemiSupervised)
read.data<-function(dat_num=NA){
  if(is.na(dat_num)){
    cat("Call function with dat_num set as:\nEthanol=1\nBreast Cancer=2\nImage=3\n")
    cat("Solubility=4\nPower=5\nNavy=6\nHouse=7\nCASP=8\nNYC=9\nSong=10\n\n")
    cat("Data sources may change over time which could result in different n,p's than ")
    cat("reported.\n\n")
    cat("Note: Some data sets require that R packages ada,mlbench, and bigRR are installed.\n")
    return(NULL)
  }
  if(dat_num==2){
    ##Breast Cancer
    library("mlbench")
    data(BreastCancer)
    x<-as.matrix(BreastCancer[,-c(1,11)])
    x<-apply(x,2,as.numeric)
    x<-impute.median(x)
    y<-BreastCancer[,11]
    y<-c(-1,1)[as.numeric(y)]
    nm<-"Breast Cancer"
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"class"
  }
  if(dat_num==4){
    ##Solubility
    library("ada")
    data(soldat)
    y<-soldat$y
    x<-as.matrix(soldat[,-73])
    x<-impute.median(x)
    nm<-"Solubility"
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"class"
  }
  if(dat_num==1){
    ##Ethonal Data
    library("bigRR")
    data(Ethanol)
    x<-as.matrix(Z.FTIR)
    x<-impute.median(x)
    y<-ethanol
    nm<-"Ethanol"
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"regress"
  }
  if(dat_num==10){
    ##Song
    temp <- tempfile()
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00203/YearPredictionMSD.txt.zip",temp)
    dat <- read.csv(unz(temp, "YearPredictionMSD.txt"))
    unlink(temp)
    y<-dat[,1]
    x<-model.matrix(~.-1,dat[,-1])
    x<-impute.median(x)
    nm<-"Song"
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"regress"
  }
  if(dat_num==8){
    ##CASP
    temp <- tempfile()
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00265/CASP.csv",temp)
    dat<-read.csv(temp,T)
    unlink(temp)
    y<-dat[,1]
    x<-as.matrix(dat[,-1])
    x<-impute.median(x)
    nm<-"CASP"
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"regress"
  }
  if(dat_num==7){
    ##California Housing
    temp <- tempfile()
    download.file("https://raw.githubusercontent.com/ageron/handson-ml/master/datasets/housing/housing.csv",temp)
    dat<-read.csv(temp,T)
    unlink(temp)
    v<-model.matrix(~ocean_proximity-1,data=dat)
    y<-dat[,9]
    x<-as.matrix(impute.median(dat[,-c(9,10)]))
    x<-cbind(x,v)
    nm<-"California Housing"
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"regress"
  }
  if(dat_num==6){
    ## Navy database
    temp <- tempfile()
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00316/UCI%20CBM%20Dataset.zip",temp)
    dat <- read.table(unz(temp, "UCI CBM Dataset/data.txt"))
    unlink(temp)
    y<-dat[,17]
    x<-dat[,-c(17:18)]
    x<-impute.median(x)
    nm<-"Navy"
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"regress"
  }
  if(dat_num==3){
    ##Image Data (7 classes are reduced to 2 similar to Tighter PAC-Bayes Bounds, NIPS 2007 version of the data).
    temp <- tempfile()
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/image/segmentation.data",temp)
    txt2<-readLines(temp)
    unlink(temp)
    temp <- tempfile()
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/image/segmentation.test",temp)
    txt1<-readLines(temp)
    unlink(temp)
    nm<-txt1[4]
    txt<-c(txt1[-c(1:5)],txt2[-c(1:5)])
    v<-strsplit(txt,",")
    vec<-sapply(v,function(i)i[1])
    x<-t(sapply(v,function(i)as.numeric(i[-1])))[,-3]
    x<-impute.median(x)
    n<-dim(x)[1]
    p<-dim(x)[2]
    y<-rep(1,n)
    y[!is.na(match(vec,c("BRICKFACE","CEMENT","WINDOW")))]=-1
    nm<-"Image"
    type<-"class"
  }
  if(dat_num==5){
    ##Power database
    library(SemiSupervised)
    data("powerplant")
    y<-powerplant[,5]
    x<-as.matrix(impute.median(powerplant[,-5]))
    n<-dim(x)[1]
    p<-dim(x)[2]
    type<-"regress"
    nm<-"Power"
  }
  if(dat_num==9){
    ##NYC database
    temp <- tempfile()
    download.file("http://taxbills.nyc/tax_bills_june15_bbls.csv",temp)
    dat<-read.csv(temp,T)
    unlink(temp)
    
    dat<-dat[,-c(3,9,11)]
    dat<-na.omit(dat)
    n<-dim(dat)[1]
    x1<-paste(dat[,2])
    tab<-table(dat$ownername)
    nm1<-names(sort(tab,dec=TRUE)[1:25])
    
    check<-function(x,nm){
      for(i in 1:length(nm)){
        if(x==nm[i])return(x)
      }
      return("Other")
    }
    x2<-sapply(1:n,function(i)check(x1[i],nm1))
    dat[,2]<-factor(x2,levels=unique(x2))
    y<-dat[,8]
    x<-dat[,-8]
    x<-model.matrix(~.-1,x)
    x<-impute.median(x)
    n<-dim(x)[1]
    p<-dim(x)[2]
    nm<-"NYC Tax"
    type<-"regress"
  }
  return(list(x=x,y=y,n=n,p=p,type=type,nm=nm))
}




