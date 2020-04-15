#' Hankel Alternative View of Koopman (HAVOK) Analysis
#'
#' @description Data-driven decomposition of chaotic time series into an intermittently
#' forced linear system. HAVOK combines delay embedding and Koopman theory to decompose
#' chaotic dynamics into a linear model in the leading delay coordinates with forcing by
#' low-energy delay coordinates. Forcing activity demarcates coherent phase space regions
#' where the dynamics are approximately linear from those that are strongly nonlinear.
#' @param xdat A vector of equally spaced measurements over time.
#' @param dt A numeric value indicating the time-lag between two subsequent time series measurements.
#' @param stackmax An integer; number of shift-stacked rows.
#' @param lambda A numeric value; sparsification threshold.
#' @param center Logical; Should \code{xdat} be centered around 0?
#' @param rmax An integer; maximum number of singular vectors to include.
#' @param rset An integer; specific number of singular vectors to include.
#' @param polyOrder An integer from 0 to 5 indicating the highest degree of polynomials
#' included in the matrix of candidate functions.
#' @param useSine A logical value indicating whether sine and cosine functions
#' of variables should be added to the library of potential candidate functions.
#' If TRUE, candidate function matrix is augmented with sine and cosine functions
#' of integer multiples 1 through 10 of all the variables in \code{yIn}.
#' @param discrete Logical; Is the underlying system discrete?
#' @return An object of class 'havok' with the following components: \itemize{
#' \item{\code{havokSS} - }{A HAVOK analysis generated state space model with its time history.}
#' \item{\code{dVrdt} - }{A matrix of first order derivatives of the first r columns of the V matrix with respect to time.}
#' \item{\code{r} - }{Estimated optimal number singular vectors to include into analysis up to \code{rmax}.}
#' \item{\code{Vr} - }{The first r columns of the V matrix of the SVD of the Hankel matrix of \code{xdat}.}
#' \item{\code{sys} - }{HAVOK model represented in state-space form.}
#' \item{\code{normTheta} - }{Normalized matrix of candidate functions obtained from \code{\link{pool_data}}.}
#' \item{\code{Xi} - }{A matrix of sparse coefficients obtained from \code{\link{sparsify_dynamics}}.}
#' \item{\code{U} - }{The U matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{sigs} - }{Values of the diagonal of the \eqn{\Sigma} matrix of the SVD of the Hankel matrix of the time series.}}
#' \item{\code{V} - }{The V matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{dt} - }{ime-lag between two subsequent time series measurements.}
#' @references S. L. Brunton, B. W. Brunton, J. L. Proctor, E. Kaiser, and J. N. Kutz,
#' "Chaos as an intermittently forced linear system," Nature Communications, 8(19):1-9, 2017.
#' @examples
#' \donttest{
#'#Lorenz Attractor
#'#Generate Data
#'##Set Lorenz Parameters
#'parameters <- c(s = 10, r = 28, b = 8/3)
#'n <- 3
#'state <- c(X = -8, Y = 8, Z =27) ##Inital Values
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
#'
#'hav <- havok(xdat = xdat, dt = dt, stackmax = 100, lambda = 0,
#'             rmax = 15, polyOrder = 1, useSine = FALSE)
#'
#'# ECG Example
#'
#'data(ECG_measurements)
#'
#'xdat <- ECG_measurements[,"channel1"]
#'dt <- ECG_measurements[2,"time"] - ECG_measurements[1,"time"]
#'stackmax <- 25
#'rmax <- 5
#'lambda <- .001
#'hav <- havok(xdat = xdat, dt = dt, stackmax = stackmax, lambda = lambda,
#'             rmax = 5, polyOrder = 1, useSine = FALSE)
#'plot(hav)
#'}
###################################

#' @export
havok <- function(xdat, dt = 1, stackmax = 100, lambda = 0, center = TRUE,
                  rmax = 15, rset = NA, polyOrder = 1, useSine = FALSE, discrete = FALSE) {

  if (center == TRUE){
    xdat <- xdat - mean(xdat)
  }

  H <- build_hankel(x = xdat, stackmax = stackmax)

  USV <- svd(H)
  U <- USV$u
  sigs <- USV$d
  V <- USV$v

  if (is.na(rmax) & is.na(rset)) {
    stop("Either 'rmax' or 'rset' must be a positive integer")
  }

  if (!is.na(rmax) & !is.na(rset)) {
    stop("Please only give values for 'rmax' or 'rset', not both")
  }

  if (!is.na(rmax)) {
    beta <- nrow(H) / ncol(H)
    thresh <- optimal_SVHT_coef(beta, FALSE) * stats::median(sigs)
    r <- length(sigs[which(sigs > thresh)])
    r <- min(rmax, r)
  }

  if (!is.na(rset)) {
    r <- rset
  }

  if (discrete == FALSE){
    dV <- compute_derivative(x = V, dt = dt, r = r)

    x <- V[3:(nrow(V) - 3), 1:r]
    dx <- dV

    Theta <- pool_data(x, r, polyOrder = polyOrder, useSine)

    normTheta <- rep(NA, dim(Theta)[2])

    for (k in 1:dim(Theta)[2]) {
      normTheta[k] <- sqrt(sum(Theta[ , k]^2))
      Theta[ , k] <- Theta[ , k]/normTheta[k]
    }

    m <- dim(Theta)[2]

    Xi <- matrix(NA,
                 nrow = nrow(sparsify_dynamics(Theta,dx[ , 1], lambda * 1)),
                 ncol = r - 1)

    for (k in 1:(r - 1)) {
      Xi[ , k] <- sparsify_dynamics(Theta, dx[ , k], lambda * k)
    }

    for (k in 1:max(dim(Xi))) {
      Xi[k, ] <- Xi[k, ] / normTheta[k]

      A <- t(Xi[2:(r + 1), 1:(r - 1)])
      B <- A[, r]
      A <- A[ , 1:(r - 1)]
      L <- 1:nrow(x)

      sys <- control::ss(A, B, pracma::eye(r - 1), 0 * B)
      HAVOK <- control::lsim(sys, x[L, r], dt * (L - 1), x[1, 1:(r - 1)])

      res <- list(HAVOK, dx, r, x, sys, Theta, Xi, U, sigs)
      names(res) <- c("havokSS", "dVrdt", "r", "Vr", "sys", "normTheta", "Xi", "U", "sigs")
      class(res) <- "havok"
      return(res)

    }

  } else if (discrete == TRUE) {
    # concatenate
    x <- V[1:(nrow(V) - 1), 1:r]
    dx <- V[2:nrow(V), 1:r]

    Xi <- pracma::mldivide(dx, x)
    B <- Xi[1:(r-1), r]
    A <- Xi[1:(r-1), 1:(r-1)]
    L <- 1:nrow(x)

    sys <- control::ss(A, B, pracma::eye(r-1), 0*B, dt)
    HAVOK <- control::lsim(sys, x[L,r], dt*(L-1), x[1, 1:r-1])

    res <- list(HAVOK, dx, r, x, sys, Xi, U, sigs, V, dt)
    names(res) <- c("havokSS", "dVrdt", "r", "Vr", "sys", "Xi", "U", "sigs", "V", "dt")
    class(res) <- "havok"
    return(res)
  }

}

# Copyright 2020 Robert Glenn Moulder Jr. & Elena Martynova
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
