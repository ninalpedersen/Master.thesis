#' Function for creating data matrix from raw data file.
#'
#' Implementation of the fixed point iteration method where a deflation scheme is used to find ICs.
#'
#' @param HSIdata csv file
#' @param rows number of rows
#'
#' @return data.mat data in matrix format
#' @examples
#' data.mat <- read.dataCube("syntetisk_aeble.csv",rows)
#'
#' @export

read.dataCube <- function(HSIdata, rows){
  dataCube <- read.csv(HSIdata, header = F, sep="")
  dataCube <- dataCube[1:dim(dataCube)[1],]
  data <- as.matrix(dataCube)

  # matrix to not summed array
  data.arr <- array(0, dim = c(rows,400,216))
  i <- 1
  r <- 0
  while(i < dim(dataCube)[1]){
    r <- r + 1
    for(c in 1:400){
      data.arr[r,c,] <- data[i,]
      i <- i + 1
    }
  }

  # array to matrix
  data.mat <- matrix(0,dim(dataCube)[1],216)
  i = 1
  while(i < dim(dataCube)[1]){
    for(r in 1:rows){
      for(c in 1:400){
        data.mat[i,] <- data.arr[r,c,]
        i <- i + 1
      }
    }
  }
  return(data.mat)
}
