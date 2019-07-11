#' fastICA implementation using a parallel scheme
#'
#' Implementation of the fixed point iteration method where a parallel scheme is used to find ICs.
#'
#' @param data.mat white data matrix
#' @param p number of components
#' @param iter number of iterations
#'
#' @return ica.arr independent component array
#' @examples
#' ica <- fastICAparallel(data,3,200)
#'
#' @export

fastICAparallel <- function(data.mat, p, iter = 200){
  w.init <- matrix(rnorm(p^2), p, p)
  f <- ncol(data.mat)
  W <- w.init
  svdW1 <- svd(W)
  W <- svdW1$u %*% diag(1/svdW1$d) %*% t(svdW1$u) %*% W
  w <- W
  it <- 1
  while (it < iter) {
    w1 <- tanh(W %*% data.mat) %*% t(data.mat)/f
    w2 <- diag(apply((1 - (tanh(W %*% data.mat))^2), 1, FUN = mean)) %*% W
    w <- w1 - w2
    svdW2 <- svd(w)
    w <- svdW2$u %*% diag(1/svdW2$d) %*% t(svdW2$u) %*% w
    W <- w
    it <- it + 1
  }
  return(W)
}



