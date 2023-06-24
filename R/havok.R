#' Hankel Alternative View of Koopman (HAVOK) Analysis
#'
#' @description Data-driven decomposition of chaotic time series into an intermittently
#' forced linear system. HAVOK combines delay embedding and Koopman theory to decompose
#' chaotic dynamics into a linear model in the leading delay coordinates with forcing by
#' low-energy delay coordinates. Forcing activity demarcates coherent phase space regions
#' where the dynamics are approximately linear from those that are strongly nonlinear.
#' @usage havok(xdat, dt = 1, stackmax = 100, lambda = 0, loops = 1, 
#'              lambda_increase = FALSE, center = TRUE,
#'              rmax = 15, rset = NA, rout = NA,
#'              polyOrder = 1, useSine = FALSE,
#'              discrete = FALSE, devMethod = c("FOCD", "GLLA"),
#'              gllaEmbed = NA, alignSVD = TRUE)
#' @param xdat A vector of equally spaced measurements over time.
#' @param dt A numeric value indicating the time-lag between two subsequent time series measurements.
#' @param stackmax An integer; number of shift-stacked rows.
#' @param lambda A numeric value; sparsification threshold.
#' @param center Logical; should \code{xdat} be centered around 0?
#' @param rmax An integer; maximum number of singular vectors to include.
#' @param rset An integer; specific number of singular vectors to include.
#' @param rout An integer or vector of integers; excludes columns of singular values from analysis.
#' @param lambda_increase Logical; should the sparsification threshold be multiplied by the sparsified column index?
#' @param polyOrder An integer from 0 to 5; if useSINDy = TRUE, the highest degree of polynomials
#' included in the matrix of candidate functions.
#' @param useSine Logical; should sine and cosine functions of variables be added to the library of potential candidate functions?
#' If TRUE, candidate function matrix is augmented with sine and cosine functions of integer multiples 1 through 10.
#' @param discrete Logical; is the underlying system discrete?
#' @param devMethod A character string; One of either \code{"FOCD"} for fourth order central difference or \code{"GLLA"} for generalized local linear approximation.
#' @param gllaEmbed An integer; the embedding dimension used for \code{devMethod = "GLLA"}.
#' @param alignSVD Logical; Whether the singular vectors should be aligned with the data.
#' @return An object of class 'havok' with the following components: \itemize{
#' \item{\code{havokSS} - }{A HAVOK analysis generated state space model with its time history.}
#' \item{\code{params} - }{A matrix of parameter values used for this function.}
#' \item{\code{dVrdt} - }{A matrix of first order derivatives of the reduced rank V matrix with respect to time.}
#' \item{\code{r} - }{Estimated optimal number singular vectors to include in analysis.}
#' \item{\code{sys} - }{HAVOK model represented in state-space form.}
#' \item{\code{normTheta} - }{Normalized matrix of candidate functions obtained from \code{\link{pool_data}}.}
#' \item{\code{Xi} - }{A matrix of sparse coefficients obtained from \code{\link{sparsify_dynamics}}.}
#' \item{\code{Vr} - }{The reduced rank V matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{Ur} - }{The reduced rank U matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{sigsr} - }{Values of the diagonal of the reduced rank \eqn{\Sigma} matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{Vr_aligned} - }{Vr truncated based upon \code{devMethod}}}
#' \item{\code{R2} - }{Squared correlation between the model predicted v_1 and v_1 exracted from SVD. A model fit estimate.}
#' @references S. L. Brunton, B. W. Brunton, J. L. Proctor, E. Kaiser, and J. N. Kutz,
#' "Chaos as an intermittently forced linear system," Nature Communications, 8(19):1-9, 2017.
#' @examples
#' \dontrun{
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
#'             lambda_increase = TRUE, loops = 10, rmax = 5, polyOrder = 1, useSine = FALSE)
#'plot(hav)
#'}
###################################

#' @export
havok <- function(xdat, dt = 1, stackmax = 100, lambda = 0,
                  lambda_increase = FALSE, loops = 1,
                  center = TRUE, rmax = 15, rset = NA, rout = NA,
                  polyOrder = 1, useSine = FALSE,
                  discrete = FALSE, devMethod = c("FOCD", "GLLA"),
                  gllaEmbed = NA, alignSVD = TRUE) {
  
  # Error catch
  if (is.na(rmax) & is.na(rset)) {
    stop("Either 'rmax' or 'rset' must be a positive integer")
  }
  
  if (!is.na(rmax) & !is.na(rset)) {
    stop("Please only give values for 'rmax' or 'rset', not both")
  }
  
  # Center time series
  if (center == TRUE){
    xdat <- xdat - mean(xdat)
  }
  
  # Create Hankel matrix
  H <- build_hankel(x = xdat, nrows = stackmax)
  
  # Perform SVD of Hankel matrix
  USV <- svd(H)
  U <- USV$u
  sigs <- USV$d
  V <- USV$v
  
  # Determine number of retained singular vectors
  if (!is.na(rmax)) {
    beta <- nrow(H) / ncol(H)
    thresh <- optimal_SVHT_coef(beta, FALSE) * stats::median(sigs)
    r <- length(sigs[which(sigs > thresh)])
    r <- min(rmax, r)
  }
  
  if (!is.na(rset)) {
    r <- rset
  }
  
  # Align and reduce rank of SVD
  if (alignSVD == TRUE){
    USV <- svd_align(H, r = r)
    U <- USV$u
    sigs <- USV$d
    V <- USV$v
  } else {
    U <- U[,1:r]
    sigs <- sigs [1:r]
    V <- V[,1:r]
  }
  
  if (!is.na(rout)) {
    V <- V[,-rout]
    r <- r - length(rout)
  }
  
  # Numerical calculation of derivatives
  devMethod <- toupper(devMethod)
  
  if (all(devMethod == c("FOCD", "GLLA"))){
    warning('Agument "devMethod" not selected. Defaulting to devMethod = "FOCD"')
    devMethod <- "FOCD"
  }
  
  if (!devMethod %in% c("FOCD", "GLLA") | length(devMethod) > 1){
    stop('devMethod must be one of either "FOCD" or "GLLA"')
  }
  
  if (discrete == FALSE){
    
    if (devMethod == "GLLA"){
      
      dV <- matrix(NA, nrow = nrow(V) - (gllaEmbed - 1), ncol = r)
      devList <- compute_derivative(V, dt, devMethod = "GLLA",
                                    gllaEmbed = gllaEmbed,
                                    gllaTau = 1,
                                    gllaOrder = 1)
      for (i in 1:r){
        dV[,i] <- devList[[i]][,2]
      }
      
      x <- V[ceiling(gllaEmbed/2):(nrow(V) - floor(gllaEmbed/2)),]
      dx <- dV
    }
    
    if (devMethod == "FOCD"){
      dV <- compute_derivative(x = V, dt = dt, r = r, devMethod = "FOCD")
      x <- V[3:(nrow(V) - 3),]
      dx <- dV
    }
    
    if (devMethod == "FINITE"){
      dV <- compute_derivative(x = V, dt = dt, r = r, devMethod = "FINITE")
      x <- V[2:nrow(V),]
      dx <- dV
    }
    
    # SINDy application and vector normalization
    Theta <- pool_data(x, nVars = r, polyOrder = polyOrder, useSine)
    
    normTheta <- rep(NA, dim(Theta)[2])
    
    for (k in 1:dim(Theta)[2]) {
      normTheta[k] <- sqrt(sum(Theta[ , k]^2))
      Theta[ , k] <- Theta[ , k]/normTheta[k]
    }
    
    m <- dim(Theta)[2]
    
    Xi <- matrix(NA,
                 nrow = nrow(sparsify_dynamics(Theta,dx[ , 1], lambda * 1, loops)),
                 ncol = r - 1)
    
    if(lambda > 0 & lambda_increase == TRUE){
      for (k in 1:(r - 1)) {
        Xi[ , k] <- sparsify_dynamics(Theta, dx[ , k], lambda * k, loops)
      }
    }
    
    if(lambda > 0 & lambda_increase == FALSE){
      for (k in 1:(r - 1)) {
        Xi[ , k] <- sparsify_dynamics(Theta, dx[ , k], lambda, loops)
      }
    }
    else {
      Xi <- pracma::mldivide(Theta, dx)
    }
    
    # State-space model reconstruction
    for (k in 1:max(dim(Xi))) {
      Xi[k, ] <- Xi[k, ] / normTheta[k]
    }
    
    A <- t(Xi[2:(r + 1), 1:(r - 1)])
    B <- A[, r]
    A <- A[ , 1:(r - 1)]
    L <- 1:nrow(x)
    
    sys <- control::ss(A, B, pracma::eye(r - 1), 0 * B)
    HAVOK <- control::lsim(sys, x[L, r], dt * (L - 1), x[1, 1:(r - 1)])
    
    R2 <- cor(x[,1],HAVOK$y[1,])^2
    
    params <- matrix(c(dt,
                       stackmax,
                       lambda,
                       lambda_increase,
                       loops,
                       center,
                       rmax,
                       rset,
                       rout,
                       polyOrder,
                       useSine,
                       discrete),
                     nrow = 12,
                     ncol = 1)
    
    colnames(params) <- "Values"
    
    rownames(params) <- c("dt",
                          "stackmax",
                          "lambda",
                          "lambda_increase",
                          "loops",
                          "center",
                          "rmax",
                          "rset",
                          "rout",
                          "polyOrder",
                          "useSine",
                          "discrete")
    
    
    res <- list(HAVOK, params, dx, r, sys, Theta, Xi, V, U, sigs, x, R2)
    names(res) <- c("havokSS", "params", "dVrdt", "r", "sys", "normTheta", "Xi", "Vr", "Ur", "sigsr", "Vr_aligned", "R2")
    class(res) <- "havok"
    return(res)
    
    
    
  } else if (discrete == TRUE) {
    # concatenate
    x <- V[1:(nrow(V) - 1),]
    dx <- V[2:nrow(V),]
    
    Xi <- pracma::mldivide(dx, x)
    B <- Xi[1:(r-1), r]
    A <- Xi[1:(r-1), 1:(r-1)]
    L <- 1:nrow(x)
    
    sys <- control::ss(A, B, pracma::eye(r-1), 0*B, dt)
    HAVOK <- control::lsim(sys, x[L,r], dt*(L-1), x[1, 1:(r-1)])
    
    R2 <- cor(x[,1],HAVOK$y[1,])^2
    
    params <- matrix(c(dt,
                       stackmax,
                       lambda,
                       lambda_increase,
                       loops,
                       center,
                       rmax,
                       rset,
                       rout,
                       polyOrder,
                       useSine,
                       discrete),
                     nrow = 12,
                     ncol = 1)
    colnames(params) <- "Values"
    rownames(params) <- c("dt",
                          "stackmax",
                          "lambda",
                          "lambda_increase",
                          "loops",
                          "center",
                          "rmax",
                          "rset",
                          "rout",
                          "polyOrder",
                          "useSine",
                          "discrete")
    
    res <- list(HAVOK, params, dx, r, sys, Xi, U, sigs, V, x, R2)
    names(res) <- c("havokSS", "params", "dVrdt", "r", "sys", "Xi", "Ur", "sigsr", "Vr","Vr_aligned", "R2")
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