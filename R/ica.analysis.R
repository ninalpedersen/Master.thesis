#' Function to perform the fastICA.
#'
#' Function that estimates the independent components.
#'
#' @param data data matrix to be whiten
#' @param p number of components
#' @param scheme choose the method for the fastICA
#' @param iter number of iterations
#'
#' @return white data for ICA analysis
#' @examples
#' data.white <- ica.analysis(data,3)
#'
#' @export


ica.analysis <- function(data, p, scheme = c("deflation","parallel"),iter = 200){

  # =====================================================
  if(scheme != "deflation" && scheme != "parallel"){
    stop("Method must be either 'defaltion' or 'parallel'")
  }
  if(missing(p)){
    stop("number of components must be chosen")
  }
  if(typeof(p) != "double"){
    stop("p must be a double")
  }
  if(missing(iter)){
    message("Number of iterations missing - default used.")
  }
  # =====================================================


  r <- nrow(data)
  c <- ncol(data)
  data.c <- scale(data, scale = FALSE) # centering
  data.c <- t(data.c)

  #whitening
  V <- data.c %*% t(data.c)/r # covariance matrix, E(x*x^T) = E*D*E^T
  s <- La.svd(V) #whitening
  D <- diag(c(1/sqrt(s$d))) # eigen values, D^(1/2)
  K <- D %*% t(s$u) # Matrix to make data.mat white with
  K <- matrix(K[1:p, ], p, c) # only he selcted number of components

  data.c.tilde <- K %*% data.c # data.mat.tilde is white data.mat - K is the pre-whitening matrix, ie. it whitens data.mat
  W <- if (scheme == "deflation"){
    fastICAdeflation(data.c.tilde,p,iter)
  } else if(scheme == "parallel"){
    fastICAparallel(data.c.tilde,p,iter)
  }

  #ica
  w <- W %*% K
  S <- w %*% data.c # source matrix
  A <- t(w) %*% solve(w %*% t(w)) #estimated mixing matrix, x = As
  data.c <- t(data.c) #pre-processed data matrix
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



