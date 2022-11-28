#' Perform Dynamic Mode Decomposition (DMD) on a matrix of snapshots of a system.
#'
#' @description Dynamic mode composition (DMD) is a method for approximating the eigenvalue and eigenvectors of
#'  the Koopman operator of a system. That is, DMD yields the eigendecomposition of the linear map of a system
#'  from time \eqn{t} to time \eqn{t+1}.
#' @param x A matrix of snapshots of a system.
#' @param y A matrix of snapshot paired to \code{x}.
#' @param r An integer; the specific number of singular vectors to include.
#' @param dt Numeric; the time-lag between two subsequent time series measurements.
#' @return An object of class 'dmd' with the following components: \itemize{
#' \item{\code{dmdResult} - }{TBD}
#' \item{\code{timeDynamics} - }{TBD}
#' \item{\code{dmdModes} - }{TBD}
#' \item{\code{dmdAmplitudes} - }{TBD}
#' \item{\code{dtEigen} - }{TBD}
#' \item{\code{ctEigen} - }{TBD}}
#' @references Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L. (2016).
#' Dynamic mode decomposition: data-driven modeling of complex systems.
#' Society for Industrial and Applied Mathematics.
#' @examples
#'library(pracma)
#'
#'# Generate data
#'x <- seq(-5, 5, length.out = 128)
#'t <- seq(0, 4*pi, length.out = 256)
#'
#' grids <- meshgrid(x, t)
#'
#' # First periodic function
#' f1xt <- t(sech(grids$X + 3)*exp(2.3i*grids$Y))
#'
#' # Second periodic function
#' f2xt <- t(2*sech(grids$X)*tanh(grids$X)*exp(2.8i*grids$Y))
#'
#' # Observed values
#' fxt <- f1xt + f2xt
#'
#' dmd(fxt, r = 2, dt = 1)
###################################

## S3 method for class "dmd"
#' @export

dmd <- function(x, y = NULL, r, dt) {

  if (!is.null(y)) {
    X <- x
    Y <- y
  } else {
    X <- x[ , 1:(ncol(x) - 1)]
    Y <- x[ , 2:ncol(x)]
  }

  Xsvd <- svd(X)
  U <- Xsvd$u
  V <- Xsvd$v
  sig <- diag(Xsvd$d)

  Ur <- U[, 1:r]
  Vr <- V[, 1:r]
  sigr <- sig[1:r, 1:r]
  sigrInv <- solve(sigr)

  Atilde <- t(Ur) %*% Y %*% Vr %*% sigrInv
  eigenAtilde <- eigen(Atilde)

  lambda <- eigenAtilde$values
  lambdaMat <- matrix(0, nrow = r, ncol = r)
  diag(lambdaMat) <- lambda

  omega <- log(lambda) / dt
  omegaMat <- matrix(0, nrow = r, ncol = r)
  diag(omegaMat) <- omega

  Phi <- Y %*% Vr %*% sigrInv %*% eigenAtilde$vectors

  svdPhi <- svd(t(Phi) %*% Phi)

  isTol <- svdPhi$d > max(.Machine$double.eps^(2/3) * svdPhi$d[1], 0)

  if (all(isTol)) {

    invRes <- svdPhi$v %*% (1/svdPhi$d * t(svdPhi$u))

  } else if (any(isTol)) {

    invRes  <- svdPhi$v[, isTol, drop=FALSE] %*% (1/svdPhi$d[isTol] * t(svdPhi$u[, isTol, drop=FALSE]))

  } else {

    invRes  <- matrix(0, nrow=ncol(t(Phi) %*% Phi), ncol=nrow(t(Phi) %*% Phi))
  }

  b <- invRes %*% t(Phi) %*% as.matrix(X[ ,1])

  xCols <- dim(X)[2]
  timeDynamics <- matrix(rep(1:xCols, r), nrow = r, ncol = xCols, byrow = TRUE) * dt
  timeDynamics <- c(b) * (exp(omegaMat %*% timeDynamics))

  Xdmd <- Phi %*% timeDynamics

  res <- list("dmdResult" = Xdmd,
              "timeDynamics" = timeDynamics,
              "dmdModes" = Phi,
              "dmdAmplitudes" = b,
              "dtEigen" = lambda,
              "ctEigen" = omega)

  class(res) <- "dmd"

  return(res)

}
