#' Hankel Alternative View of Koopman (HAVOK) Analysis
#'
#' @description Coefficient determining optimal location of Hard Threshold for Matrix
#'  Denoising by Singular Values Hard Thresholding when noise level is known or
#'  unknown.  Recreation of matlab code by Matan Gavish and David Donoho.
#' @param xdat A number.
#' @param dt A number.
#' @param stackman A number.
#' @param lambda A number.
#' @param rmax A number.
#' @param polyorder A number.
#' @param useSine A number.
#' @param n A number.
#' @return  A matrix of sparse coefficients
#' @examples
#'#Generate Data
#'library(deSolve)
#'##Set Lorenz Parameters
#'parameters <- c(s = 10, r = 28, b = 8/3)
#'n <- 3
#'state <- c(X = -8, Y = 8, Z =2 7) ##Inital Values
#'
#'#Intergrate
#'dt <- 0.001
#'tspan <- seq(dt, 200, dt)
#'N <- length(tspan)
#'
#'Lorenz <- function(t, state, parameters) {
#'  with(as.list(c(state, parameters)), {
#'    dX <- s * (Y - X)
#'    dY <- X * (r - Z) - Y
#'    dZ <- X * Y - b * Z
#'    list(c(dX, dY, dZ))
#'  })
#'}
#'
#'out <- ode(y = state, times = tspan, func = Lorenz, parms = parameters, rtol = 1e-12, atol = 1e-12)
#'xdat <- out[, "X"]
#'t <- out[, "time"]
#'hav <- <- havok(xdat = xdat, dt = dt, stackmax = 100, lambda = 0, rmax = 15, polyorder = 1, useSine = FALSE, n = 1)
#'#'
###################################

havok <- function(xdat, dt = 1, stackmax = 100, lambda = 0,
                  rmax = 15, polyorder = 1, useSine = FALSE, n = 1) {

  H <- matrix(0, nrow = stackmax, ncol = length(xdat) - stackmax)

  for (k in 1:stackmax) {
    H[k,] <- xdat[k:(length(xdat) - stackmax - 1 + k)]
  }

  USV <- svd(H)
  U <- USV$u
  sigs <- USV$d
  V <- USV$v
  beta <- nrow(H) / ncol(H)
  thresh <- optimal_SVHT_coef(beta, FALSE) * median(sigs)
  r <- length(sigs[which(sigs > thresh)])
  r <- min(rmax, r)

  # COMPUTE DERIVATIVES
  dV <- matrix(0, nrow = max(dim(V)) - 5, ncol = r)
  for (i in 3:(nrow(V) - 3)) {
    for (k in 1:r) {
      dV[(i-2),k] <- (1 / (12 * dt)) *
        (-V[i + 2, k] + 8 * V[i + 1, k] -
           8 * V[i - 1, k] + V[i - 2, k])
    }
  }


  # concatenate
  x <- V[3:(nrow(V) - 3), 1:r]
  dx <- dV

  polyorder <- 1
  Theta <- pool_data(x, r, 1, useSine)


  # normalize columns of Theta (required in new time-delay coords)
  normTheta <- rep(NA, dim(Theta)[2])

  for (k in 1:dim(Theta)[2]) {
    normTheta[k] <- sqrt(sum(Theta[ , k]^2))
    Theta[ , k] <- Theta[ , k]/normTheta[k]
  }

  m <- dim(Theta)[2]

  #compute Sparse regression: sequential least squares
  #requires different lambda parameters for each column
  Xi <- matrix(NA,
             nrow = nrow(sparsify_dynamics(Theta,dx[ , 1], lambda * 1, 1)),
             ncol = r - 1)

  for (k in 1:(r - 1)) {
    Xi[ , k] <- sparsify_dynamics(Theta, dx[ , k], lambda * k, 1)  # lambda = 0 gives better results
  }

  Theta <- pool_data(x, r , 1, useSine)

  for (k in 1:max(dim(Xi))) {
    Xi[k, ] <- Xi[k, ] / normTheta[k]
  }

  A <- t(Xi[2:(r + 1), 1:(r - 1)])
  B <- A[, r]
  A <- A[ , 1:(r - 1)]
  L <- 1:nrow(x)
  #Need State space models!
  sys <- control::ss(A, B, pracma::eye(r - 1), 0 * B)
  HAVOK <- control::lsim(sys, x[L, r], dt * (L - 1), x[1, 1:(r - 1)])

  res <- list(HAVOK, r, x, sys, Theta, Xi, U, sigs, V)
  names(res) <- c("havok", "r", "x", "sys", "theta", "Xi", "U", "sigs", "V")
  return(res)
}

