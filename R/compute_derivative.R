#' Compute derivative using fourth order central difference
#'
#' @description Estimates the derivative of a time series using fourth order central difference.
#' @param x A vector of measurements over time.
#' @param dt A number representing the change in time between successive measurements.
#' @param r A number representing the number of time series in x used to calculate dXdt.
#' @return \code{dXdt - } A vector or matrix of first order derivatives of the variables of interest with respect to time.
#' @examples
#' \donttest{
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
#'dXdt <- compute_derivative(x = xdat, dt = dt)
#'}
###################################

#' @export
compute_derivative <- function(x, dt, r = min(dim(as.matrix(x)))){

  if (is.vector(x)){
    x <- as.matrix(x)
  }

  dXdt <- matrix(NA, nrow = max(dim(x)) - 5, ncol = r)

  for (i in 3:(nrow(x) - 3)) {
    for (k in 1:r) {
      dXdt[(i-2),k] <- (1 / (12 * dt)) *
        (-x[i + 2, k] + 8 * x[i + 1, k] -
           8 * x[i - 1, k] + x[i - 2, k])
    }
  }

  if (min(dim(dXdt)) == 1){
    return(as.vector(dXdt))
  } else {
    return(dXdt)
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
