#' fastICA implementation using a deflation scheme
#'
#' Implementation of the fixed point iteration method where a deflation scheme is used to find ICs.
#'
#' @param data.mat white data matrix
#' @param p number of components
#' @param iter number of iterations
#'
#' @return ica.arr independent component array
#' @examples
#' ica <- fastICAdeflation(data,3,200)
#'
#' @export

fastICAdeflation <- function(data.mat, p, iter = 200){
  w.init <- matrix(rnorm(p^2), p, p)
  c <- ncol(data.mat)
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
    while (iter < iter) {
      xg <- data.mat * matrix(tanh(t(w) %*% data.mat), p, c, byrow = TRUE)
      w1 <- apply(xg, 1, FUN = mean) - mean(1 - (tanh(t(w) %*% data.mat))^2) * w
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
  return(W)
}



