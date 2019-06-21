#' fastICA implementation using a deflation scheme
#'
#' Implementation of the fixed point iteration method where a deflation scheme is used to find ICs.
#'
#' @param data.mat data matrix 80.000x216
#' @param p number of components
#' @param max maximum number of iterations
#'
#' @return ica.arr independent component array
#' @examples
#' ica <- fastICAdeflation(data,3,200)
#'
#' @export

fastICAdeflation <- function(data.mat, p, max = 200){
  r <- nrow(data.mat) #test
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

  #deflation
  c <- ncol(data.mat.tilde)
  W <- matrix(0, p, p)

  for (i in 1:p) {
    w <- matrix(w.init[i, ], p, 1)
    if (i > 1) {
      temp.w <- w
      temp.w[1:length(w)] <- 0
      for (j in 1:(i - 1)){
        temp.w <- temp.w + sum(w * W[j, ]) * W[j, ]
      }
      w <- w - temp.w
    }
    w <- w/sqrt(sum(w^2))
    conv <- 1
    iter <- 1
    while (iter < max) {
      xg <- data.mat.tilde * matrix(tanh(t(w) %*% data.mat.tilde), p, c, byrow = TRUE)
      w1 <- apply(xg, 1, FUN = mean) - mean(1 - (tanh(t(w) %*% data.mat.tilde))^2) * w
      w1 <- matrix(w1, p, 1)
      if (i > 1) {
        temp.w <- w1
        temp.w[1:length(w1)] <- 0
        for (j in 1:(i - 1)){
          temp.w <- temp.w + sum(w1 * W[j, ]) * W[j, ]
        }
        w1 <- w1 - temp.w
      }
      w1 <- w1/sqrt(sum(w1^2))
      w <- matrix(w1, p, 1)
      iter <- iter + 1
    }
    W[i, ] <- w
  }

  #ica
  w <- W %*% K
  S <- w %*% data.mat # source matrix
  A <- t(w) %*% solve(w %*% t(w)) #estimated mixing matrix, x = As
  data.mat <- t(data.mat) #pre-processed data matrix
  K <- t(K) #pre-whitening matrix
  W <- t(W) #estimated un-mixing matrix, y = WAs = Wx
  A <- t(A) #
  S <- t(S) #estimated

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
  return(ica <- list("arr" = ica.arr, "K" = K, "W" = W))
}



