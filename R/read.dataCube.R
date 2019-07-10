#' Function for creating data matrix from raw data file.
#'
#' Creates a data matrix from a .csv file containing HSI data. The input data must be (rowsxcols)xbands
#'
#' @param HSIdata csv file
#' @param rows number of rows
#' @param cols number of columns
#' @param bands number of bands
#'
#' @return data.mat data in matrix format
#' @examples
#' data.mat <- read.dataCube("syntetisk_aeble.csv",rows,cols,bands)
#'
#' @export

read.dataCube <- function(HSIdata, rows = 200, cols = 400, bands = 216){

  # =====================================================
  if(typeof(rows) != "double"){
    stop("rows must be a double")
  }
  if(typeof(cols) != "double"){
    stop("cols must be a double")
  }
  if(typeof(bands) != "double"){
    stop("bands must be a double")
  }
  if(typeof(HSIdata) != "character"){
    message("'HSIdata' must be a character string or connection")
  }
  # =====================================================




  dataCube <- read.csv(HSIdata, header = F, sep="")
  dataCube <- dataCube[1:dim(dataCube)[1],]
  data <- as.matrix(dataCube)

  # matrix to not summed array
  data.arr <- array(0, dim = c(rows,cols,bands))
  i <- 1
  r <- 0
  while(i < dim(dataCube)[1]){
    r <- r + 1
    for(c in 1:cols){
      data.arr[r,c,] <- data[i,]
      i <- i + 1
    }
  }

  # array to matrix
  data.mat <- matrix(0,dim(dataCube)[1],bands)
  i = 1
  while(i < dim(dataCube)[1]){
    for(r in 1:rows){
      for(c in 1:cols){
        data.mat[i,] <- data.arr[r,c,]
        i <- i + 1
      }
    }
  }
  return(data.mat)
}
