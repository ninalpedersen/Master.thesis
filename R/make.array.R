#' Function converting 80000xbands matrix to 200x400xbands array.
#'
#' @param mat matrix with 80000xbands dimension
#' @param bands number of bands
#'
#' @return arr array with 200x400xbands dimension
#' @examples
#' arr <- make.array(mat,bands)
#'
#' @export

make.array <- function(mat, rows, cols, bands){

  # =====================================================
  if(missing(rows)){
    stop("Number of rows missing!")
  }
  if(missing(cols)){
    stop("Number of cols missing!")
  }
  if(missing(bands)){
    stop("Number of bands missing!")
  }
  if(typeof(bands) != "double"){
    stop("bands must be a double!")
  }
  # =====================================================
  arr <- array(0, dim = c(rows,cols,bands))
  i <- 1
  r <- 0
  if(bands == 1){
    while(i < rows*cols){
      r <- r + 1
      for(c in 1:cols){
        arr[r,c,] <- mat[i,]
        i <- i + 1
      }
    }
  }else{
    while(i < rows*cols){
      r <- r + 1
      for(c in 1:cols){
        arr[r,c,] <- mat[i,]
        i <- i + 1
      }
    }
  }
  return(arr)
}
