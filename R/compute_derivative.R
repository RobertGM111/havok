#' Compute derivative using fourth order central difference
#'
#' @description Estimates the derivative of a time series using fourth order central difference.
#' @param x A vector of measurements over time.
#' @param dt A number representing the change in time between successive measurements.
#' @param r A number representing the number of time series in x used to calculate \code{dXdt}.
#' @param gllaEmbed NEED DESCRIPTION
#' @param gllaTau NEED DESCRIPTION
#' @param gllaOrder NEED DESCRIPTION
#' @param devMethod NEED DESCRIPTION - either "FOCD" or "GLLA"
#' @return \code{dXdt - } A vector or matrix of first order derivatives of the variables of interest with respect to time.
#' @examples
#' \donttest{
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
#'dXdt <- compute_derivative(x = xdat, dt = dt)
#'}
###################################

#' @export
compute_derivative <- function(x, dt, r = min(dim(as.matrix(x))),
                               gllaEmbed = 3, gllaTau = 1, gllaOrder = 1,
                               devMethod = "FOCD"){

  if (!devMethod %in% c("FOCD", "GLLA")){
    stop('devMethod must be either "FOCD" or "GLLA"')
  }

  if (devMethod == "FOCD"){

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

  if (devMethod == "GLLA"){

    L <- rep(1, gllaEmbed)

    for(i in 1:gllaOrder) {
      L <- cbind(L, (((c(1:gllaEmbed) - mean(1:gllaEmbed)) *
                        gllaTau * dt)^i) / factorial(i))
    }

    W <- L %*% solve(t(L) %*% L)

    minLen <- (gllaTau + 1 + ((gllaEmbed - 2) * gllaTau))

    if (gllaEmbed < 2 | is.na(gllaEmbed) | gllaTau < 1 | is.na(gllaTau)){
      return(NA)
    }

    embedMat <- function(y, gllaEmbed, gllaTau){

      embeddedMatrix <- matrix(NA, length(y) + (gllaEmbed*gllaTau), gllaEmbed + 1)

      tLen <- length(y) - minLen

      embeddedMatrix[1:(1 + tLen), 1] <- 1

      for (j in 1:gllaEmbed) {
        k <- 1 + ((j - 1) * gllaTau)
        embeddedMatrix[1:(1 + tLen), j + 1] <- y[k:(k + tLen)]
      }

      tRow <- tLen + 2

      return(embeddedMatrix[1:(tRow - 1), 2:(gllaEmbed + 1)])
    }


    if (r == 1){
      dXdt <- embedMat(x, gllaEmbed = gllaEmbed, gllaTau = gllaTau) %*% W
      return(dXdt)
    } else {
      dXdt <- rep(list(NA), r)
      for (i in 1:r) {
        dXdt[[i]] <-  embedMat(x[,i], gllaEmbed = gllaEmbed, gllaTau = gllaTau) %*% W
      }
      return(dXdt)
    }

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
