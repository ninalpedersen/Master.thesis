#' Infomax implementation
#'
#' Implementation of the infomax algorithm to find ICs.
#'
#' @param data.mat data matrix where the columns are bands
#' @param p number of components
#' @param max maximum number of iterations
#'
#' @return ica.arr independent component array
#' @examples
#' ica <- infomax(data,3,200)
#'
#' @export

infomax <- function(data.mat, p, maxit = 200){
  w.init <- diag(p)
  alpha <- 1
  L.init <- ncol(data.mat)*log(abs(det(w.init)))+sum(-w.init%*%data.mat-2*log(1+exp(-w.init%*%data.mat)))
  deltaW <- solve(t(w.init))+((1-2*(1/(1+exp(-w.init%*%data.mat))))%*%t(data.mat))/ncol(data.mat)
  iter <- 1
  W <- w.init
  while (iter < maxit) {
    alpha <- 2*alpha
    w <- W + alpha*deltaW
    L <- ncol(data.mat)*log(abs(det(w)))+sum(-w%*%data.mat-2*log(1+exp(-w%*%data.mat)))
    while (L <= L.init) {
      alpha <- alpha/2
      w <- W+alpha*deltaW
      L <- ncol(data.mat)*log(abs(det(w)))+sum(-w%*%data.mat-2*log(1+exp(-w%*%data.mat)))
    }
    deltaW <- solve(t(w))+((1-2*(1/(1+exp(-w%*%data.mat))))%*%t(data.mat))/ncol(data.mat)
    L.init <- L
    W <- w
    iter <- iter + 1
  }
  return(W)
}



