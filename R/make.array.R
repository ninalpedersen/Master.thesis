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

make.array <- function(mat, bands){

  # =====================================================
  if(missing(bands)){
    stop("Number of bands missing!")
  }
  if(typeof(bands) != "double"){
    stop("bands must be a double!")
  }
  if(typeof(mat) != "matrix"){
    stop("'matrix' must be a matrix!")
  }
  # =====================================================
  arr <- array(0, dim = c(200,400,bands))
  i <- 1
  r <- 0
  if(bands == 1){
    while(i < 80000){
      r <- r + 1
      for(c in 1:400){
        arr[r,c,] <- mat[i]
        i <- i + 1
      }
    }
  }else{
    while(i < 80000){
      r <- r + 1
      for(c in 1:400){
        arr[r,c,] <- mat[i,]
        i <- i + 1
      }
    }
  }
  return(mat)
}
