#' Compute numeric derivatives of a time series.
#'
#' @description Estimates the derivative of a time series using numeric methods.
#' @usage compute_derivative(x, dt, r = min(dim(as.matrix(x))),
#'                           devMethod = c("FOCD", "GLLA", "FINITE"),
#'                           gllaEmbed = NA, gllaTau = NA, gllaOrder = NA)
#' @param x A vector or matrix of measurements over time.
#' @param dt Numeric; the change in time between successive measurements.
#' @param r An integer; the number of time series in \code{x} used to calculate \code{dXdt}.
#' @param devMethod A character string. One of either \code{"FOCD"} for fourth order central difference, \code{"GLLA"} for generalized local linear approximation,
#' or \code{"FINITE"} for simple finite difference.
#' @param gllaEmbed An integer; t embedding dimension used for \code{devMethod = "GLLA"}.
#' @param gllaTau An integer; the time delay used for \code{devMethod = "GLLA"}.
#' @param gllaOrder An integer; the embedding dimension used for \code{devMethod = "GLLA"}.
#' @return
#' \itemize{
#' \item If \code{devMethod = "FOCD"} - returns a vector or matrix of first order derivatives the first \code{r} columns of \code{x} with respect to time.
#' \item If \code{devMethod = "GLLA"} - returns a matrix or list of matrices of derivatives of the first \code{r} columns of \code{x} with respect to time.
#' Derivatives returned by \code{devMethod = "GLLA"} are up to order \code{"gllaOrder"} and include the "\eqn{0^{th}}" derivative.
#' \item If \code{devMethod = "FINITE"} - returns a vector or matrix of first order derivatives the first \code{r} columns of \code{x} with respect to time.
#' }
#' @references Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010). Generalized local linear approximation of derivatives
#'  from time series.
#'  In Chow S, Ferrer E, and Hsieh F, editors, Statistical methods for modeling human dynamics: An interdisciplinary dialogue.
#' @examples
#'data(ECG_measurements)
#'
#'xdat <- ECG_measurements[,"channel1"]
#'dt <- ECG_measurements[2,"time"] - ECG_measurements[1,"time"]
#'
#'compute_derivative(x = xdat, dt = dt, devMethod = "FOCD")
#'
#'\dontrun{
#'compute_derivative(x = xdat, dt = dt, devMethod = "GLLA", gllaEmbed = 5, gllaTau = 1, gllaOrder = 1)
#'}
###################################

#' @export
compute_derivative <- function(x, dt, r = min(dim(as.matrix(x))),
                               devMethod = c("FOCD", "GLLA", "FINITE"),
                               gllaEmbed = NA, gllaTau = NA, gllaOrder = NA){

  devMethod <- toupper(devMethod)

  if (all(devMethod == c("FOCD", "GLLA", "FINITE"))){
    warning('Agument "devMethod" not selected. Defaulting to devMethod = "FOCD"')
    devMethod <- "FOCD"
  }

  if (!devMethod %in% c("FOCD", "GLLA", "FINITE") | length(devMethod) > 1){
    stop('devMethod must be one of either "FOCD", "GLLA", or "FINITE')
  }

  if (is.vector(x)){
    x <- as.matrix(x)
  }

  if (devMethod == "FOCD"){

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

  if (devMethod == "FINITE"){

    dXdt <- matrix(NA, nrow = max(dim(x)) - 1, ncol = r)

    for (i in 2:nrow(x)) {
      for (k in 1:r) {
        dXdt[(i-1),k] <- (x[i, k] - x[(i-1), k]) / dt
      }
    }

    if (min(dim(dXdt)) == 1){
      return(as.vector(dXdt))
    } else {
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
