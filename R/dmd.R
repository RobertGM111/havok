#' Dynamic Mode Decomposition (DMD)
#'
#' @description DMD Description
#' @param x A matrix of snapshots
#' @param y A matrix of snapshot paired to \code{x}
#' @param r An integer; specific number of singular vectors to include.
#' @param dt A numeric value indicating the time-lag between two subsequent time series measurements.
#' @return An object of class 'dmd' with the following components: \itemize{
#' \item{\code{dmdResult} - }{A matrix representing the DMD recreation of \code{x} or \code{xy} if paired.}
#' \item{\code{timeDynamics} - }{A matrix of parameter values used for this function.}
#' \item{\code{dmdModes} - }{A matrix of first order derivatives of the first r columns of the V matrix with respect to time.}
#' \item{\code{dmdAmplitudes} - }{Estimated optimal number singular vectors to include into analysis up to \code{rmax}.}
#' \item{\code{dtEigen} - }{The first r columns of the V matrix of the SVD of the Hankel matrix of \code{xdat}.}
#' \item{\code{ctEigen} - }{HAVOK model represented in state-space form.}
#' @references Brunton, S. L., Proctor, J. L., & Kutz, J. N. (2016). Discovering
#' governing equations from data by sparse identification of nonlinear dynamical
#' systems. Proceedings of the National Academy of Sciences, 113(15), 3932-3937.
#' @examples
#' \dontrun{
#' sparsify_dynamics(Theta, dXdt, lambda, n)
#' sparsify_dynamics(Theta, dXdt, 0.1, 10)
#' sparsify_dynamics(pool_data(yIn, 15, 5, TRUE), dXdt, 0, 15)
#' }
###################################

## S3 method for class "dmd"
#' @export
mldivide <- function (A, B, pinv = TRUE) {
  stopifnot(is.numeric(A) || is.complex(A), is.numeric(B) ||
              is.complex(B))
  if (is.vector(A))
    A <- as.matrix(A)
  if (is.vector(B))
    B <- as.matrix(B)
  if (nrow(A) != nrow(B))
    stop("Matrices 'A' and 'B' must have the same number of rows.")
  if (pinv) {
    pinv(t(A) %*% A) %*% t(A) %*% B
  }
  else {
    qr.solve(A, B)
  }
}

pinv <- function (A, tol = .Machine$double.eps^(2/3)) {
  stopifnot(is.numeric(A) || is.complex(A), length(dim(A)) == 2, is.matrix(A))

  s <- svd(A)

  p <- ( s$d > max(tol * s$d[1], 0) )
  if (all(p)) {
    mp <- s$v %*% (1/s$d * t(s$u))
  } else if (any(p)) {
    mp <- s$v[, p, drop=FALSE] %*% (1/s$d[p] * t(s$u[, p, drop=FALSE]))
  } else {
    mp <- matrix(0, nrow=ncol(A), ncol=nrow(A))
  }

  return(mp)
}



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

  b <- mldivide(Phi, X[ ,1], pinv = TRUE)

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
