#' fastICA implementation using a parallel scheme
#'
#' Implementation of the fixed point iteration method where a parallel scheme is used to find ICs.
#'
#' @param data.mat data matrix 80.000x216
#' @param p number of components
#' @param max maximum number of iterations
#'
#' @return ica.arr independent component array
#' @examples
#' ica <- fastICAparallel(data,3,200)
#'
#' @export

fastICAparallel <- function(data.mat, p, max = 200){
  r <- nrow(data.mat)#hehe
  c <- ncol(data.mat)
  w.init <- matrix(rnorm(p^2), p, p)
  data.mat <- scale(data.mat, scale = FALSE) # centering
  data.mat <- t(data.mat)

  #whitening
  V <- data.mat %*% t(data.mat)/r # covariance matrix, E(x*x^T) = E*D*E^T
  s <- La.svd(V) #whitening
  D <- diag(c(1/sqrt(s$d))) # eigen values, D^(1/2)
  K <- D %*% t(s$u) # Matrix to make data.mat white with
  K <- matrix(K[1:p, ], p, c) # only he selcted number of components

  data.mat.tilde <- K %*% data.mat # data.mat.tilde is white data.mat - K is the pre-whitening matrix, ie. it whitens data.mat

  #parallel

  #ica
  ##
  ##

  ica.arr <- array(0, dim = c(rows,400,p))
  i <- 1
  r <- 0
  while(i < dim(data.mat)[1]){
    r <- r + 1
    for(c in 1:400){
      ica.arr[r,c,] <- S[i,]
      i <- i + 1
    }
  }
  #return(ica <- list("arr" = ica.arr, "K" = K, "W" = W))
}



