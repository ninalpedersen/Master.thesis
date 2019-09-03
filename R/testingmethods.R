###PACKAGES###
library(fields)
library(Master.thesis)
library(tictoc)
####

setwd("~/Dropbox/Dokumenter/SDU/Newtec/Speciale")

rows = 200
cols = 400
bands = 216
comp <- 3
data.mat <- read.dataCube("simdata.csv",rows,cols,bands)
image.plot(data.mat[])

##TIME TEST##
infomaxT <- c(0,1:19)
for(i in 1:20){
  start_time <- Sys.time()
  ica <- ica.analysis(data.mat,comp,"infomax",50)
  end_time <- Sys.time()
  infomaxT[i] <- end_time - start_time
  rm(ica)
}
deflationT <- c(0,1:19)
for(i in 1:20){
  start_time <- Sys.time()
  ica <- ica.analysis(data.mat,comp,"deflation",50)
  end_time <- Sys.time()
  deflationT[i] <- end_time - start_time
  rm(ica)
}
parallelT <- c(0,1:19)
for(i in 1:20){
  start_time <- Sys.time()
  ica <- ica.analysis(data.mat,comp,"parallel",50)
  end_time <- Sys.time()
  parallelT[i] <- end_time - start_time
  rm(ica)
}

plot(deflationT, type='n', main="", xlab="x", ylab="y",ylim = c(0,max(deflationT,parallelT,infomaxT)))
grid()
lines(infomaxT,col="darkblue" , lwd=3 , pch=18 , type="b" )
lines(deflationT,col="darkred" , lwd=3 , pch=18 , type="b" )
lines(parallelT,col="darkgreen" , lwd=3 , pch=18 , type="b" )

#####

#SUM OF SQUARES


