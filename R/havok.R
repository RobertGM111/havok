#' Hankel Alternative View of Koopman (HAVOK) Analysis
#'
#' @description Data-driven decomposition of chaotic time series into an intermittently
#' forced linear system. HAVOK combines delay embedding and Koopman theory to decompose
#' chaotic dynamics into a linear model in the leading delay coordinates with forcing by
#' low-energy delay coordinates. Forcing activity demarcates coherent phase space regions
#' where the dynamics are approximately linear from those that are strongly nonlinear.
#' @param xdat A vector of equally spaced measurements over time.
#' @param dt A numeric value indicating the time-lag between two subsequent time series measures.
#' @param stackmax An integer; number of shift-stacked rows.
#' @param lambda A numeric value; sparsification threshold.
#' @param center Logical; Should \code{xdat} be centered around 0?
#' @param rmax An integer; maximum number of singular vectors to include.
#' @param polyOrder An integer from 0 to 5 indicating the highest degree of polynomials
#' included in the matrix of candidate functions.
#' @param useSine A logical value indicating whether sine and cosine functions
#' of variables should be added to the library of potential candidate functions.
#' If TRUE, candidate function matrix is augmented with sine and cosine functions
#' of integer multiples 1 through 10 of all the variables in \code{yIn}.
#' @return  A matrix of sparse coefficients
#' @references S. L. Brunton, B. W. Brunton, J. L. Proctor, E. Kaiser, and J. N. Kutz,
#' "Chaos as an intermittently forced linear system," Nature Communications, 8(19):1-9, 2017.
#' @examples
#' \donttest{
#'#Lorenz Attractor
#'#Generate Data
#'library(deSolve)
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
#'t <- out[, "time"]
#'hav <- havok(xdat = xdat, dt = dt, stackmax = 100, lambda = 0,
#'             rmax = 15, polyOrder = 1, useSine = FALSE)
#'
#'# EEG Example
#'
#'
#'
#'
#'
#'
#'
#'
#'}
###################################

#' @export
havok <- function(xdat, dt = 1, stackmax = 100, lambda = 0, center = TRUE,
                  rmax = 15, polyOrder = 1, useSine = FALSE, discrete = FALSE) {

  if (center == TRUE){
    xdat <- xdat - mean(xdat)
  }

  H <- build_hankel(x = xdat, stackmax = stackmax)

  USV <- svd(H)
  U <- USV$u
  sigs <- USV$d
  V <- USV$v
  beta <- nrow(H) / ncol(H)
  thresh <- optimal_SVHT_coef(beta, FALSE) * stats::median(sigs)
  r <- length(sigs[which(sigs > thresh)])
  r <- min(rmax, r)

  # COMPUTE DERIVATIVES
  dV <- compute_derivative(x = V, dt = dt, r = r)

  # concatenate
  x <- V[3:(nrow(V) - 3), 1:r]
  dx <- dV

  Theta <- pool_data(x, r, polyOrder = polyOrder, useSine)

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
             nrow = nrow(sparsify_dynamics(Theta,dx[ , 1], lambda * 1)),
             ncol = r - 1)

  for (k in 1:(r - 1)) {
    Xi[ , k] <- sparsify_dynamics(Theta, dx[ , k], lambda * k)
  }

  for (k in 1:max(dim(Xi))) {
    Xi[k, ] <- Xi[k, ] / normTheta[k]
  }

  A <- t(Xi[2:(r + 1), 1:(r - 1)])
  B <- A[, r]
  A <- A[ , 1:(r - 1)]
  L <- 1:nrow(x)

  sys <- control::ss(A, B, pracma::eye(r - 1), 0 * B)
  HAVOK <- control::lsim(sys, x[L, r], dt * (L - 1), x[1, 1:(r - 1)])

  res <- list(HAVOK, dx, r, x, sys, Theta, Xi, U, sigs)
  names(res) <- c("havok", "dxdt", "r", "Vr", "sys", "normTheta", "Xi", "U", "sigs")
  class(res) <- "havok"
  return(res)
}

