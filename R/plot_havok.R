#' Plotting functions for the havok package of class “havok"
#'
#' @description Generic plotting function for object of class ("havok")
#' @param x A "havok" object.
#' @param what \itemize {
#' \item{interactive - }{An interactive plotting function.}
#' \item{reconstruction - }{Reconstruction of the major component of a time-series.}
#' \item{forcing - }{Forcing vector derived from HAVOK.}
#' \item{both - }{A combination of 'reconstruction' and 'forcing'.}
#' \item{U-modes - }{U modes of the reconstructed time series.}
#' \item{embedded - }{A 2D reconstruction of the attractor.}
#' \item{nonlinear - }{A 2D reconstruction of the attractor with nonlinear regions colored red.}
#' @param ... Other calls to plot
#' @examples
#'
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
#'t <- out[, "time"]
#'hav <- havok(xdat = xdat, dt = dt, stackmax = 100, lambda = 0,
#'             rmax = 15, polyOrder = 1, useSine = FALSE)
#'
#'plot(hav)
#'}
###################################

## S3 method for class "havok"
#' @export
plot.havok <- function(x, what = "interactive", ...) {

  if (what == "interactive"){

    cat("--- Please select a plot type by number ---\n
        Plot Types:
        1 - Reconstruction     4 - U-modes
        2 - Forcing            5 - Embedded Attractor
        3 - Both               6 - Nonlinear Region")

    repeat {
      whatPlot <- readline("Please select a number (press esc to exit): ")

      if (!whatPlot %in% 1:6){
        stop("Please pick a number between 1 and 5")
      }

      if (whatPlot == 1){
        graphics::par(mai = c(0.9, 0.8, 0.1, 0.1))
        graphics::plot(x$havokSS$t, x$Vr[,1], type = "l", xlab = "Time", ylab = "Value", ...)
        graphics::par(mai = c(1.02, 0.82, 0.82, 0.42))
      }

      if (whatPlot == 2){
        graphics::par(mai = c(0.9, 0.8, 0.1, 0.1))
        graphics::plot(x$havokSS$t, x$Vr[,x$r], type = "l", xlab = "Time", ylab = "Forcing", ...)
        graphics::par(mai = c(1.02, 0.82, 0.82, 0.42))
      }

      if (whatPlot == 3){
        graphics::par(mfrow=c(2,1), mai = c(0.5, 1.0, 0.1, 0.1))

        graphics::plot(x$havokSS$t, x$Vr[,1], type = "l", xlab = NA, ylab = "Value", ...)

        graphics::par(mai = c(1, 1.0, 0.1, 0.1))

        graphics::plot(x$havokSS$t, x$Vr[,x$r], type = "l", xlab = "Time", ylab = "Forcing", ...)

        graphics::par(mfrow=c(1,1), mai = c(1.02, 0.82, 0.82, 0.42))
      }

      if (whatPlot == 4){
        graphics::par(mai = c(0.9, 0.8, 0.1, 0.1))
        plotBuff <- (.1 * (max(x$U[, 1:x$r]) - min(x$U[, 1:x$r])))
        graphics::plot(x$U[,1], ylab = "Ur", xlab = "Time",
                       type = "l",
                       xaxt = "n",
                       col = grDevices::rainbow(x$r)[1],
                       ylim = c(min(x$U[,1:x$r]) - (.1*plotBuff), max(x$U[,1:x$r]) + (.1*plotBuff)))
        graphics::axis(1, at = seq(1, ncol(x$U), length.out = 10),
             labels = round(seq(1, max(x$havokSS$t), length.out = 10)))

        for (i in 1:x$r){
          graphics::lines(x$U[,i], col = grDevices::rainbow(x$r, alpha = 1/sqrt(i))[i])
        }

        if (x$r < 6){
          graphics::legend("topleft",
                 fill = grDevices::rainbow(x$r, alpha = 1/sqrt(1:x$r)),
                 legend = paste("r = ", 1:x$r, sep = ""))
        }

        if (x$r >= 6){
          graphics::legend("topleft",
                 fill = grDevices::rainbow(x$r, alpha = 1/sqrt(1:x$r))[c(1:6,x$r)],
                 legend = c(paste("r = ", 1:5, sep = ""),
                            "...",
                            paste("r = ", x$r, sep = ""))
                 )
        }
        graphics::par(mai = c(1.02, 0.82, 0.82, 0.42))

      }

      if (whatPlot == 5){
        graphics::par(mai = c(0.9, 0.8, 0.1, 0.1))
        graphics::plot(x$Vr[,1], x$Vr[,2], type = "l", xlab = "V1", ylab = "V2", ...)
        graphics::par(mai = c(1.02, 0.82, 0.82, 0.42))
      }

      if (whatPlot == 6){
        havForce <- active_forcing(x)
        graphics::par(mai = c(0.9, 0.8, 0.1, 0.1))
        graphics::plot(x$Vr[,1], x$Vr[,2], type = "l", xlab = "V1", ylab = "V2", ...)
        segments(head(x$Vr[,1], -1), head(x$Vr[,2], -1), x$Vr[,1][-1], x$Vr[,2][-1], ifelse(havForce$active==1,"red","black"))
        graphics::par(mai = c(1.02, 0.82, 0.82, 0.42))
      }

    }
  }


  if (what == "reconstruction"){
    graphics::plot(x$havokSS$t, x$Vr[,1], type = "l", xlab = NA, ylab = "Value", ...)
  }

  if (what == "forcing"){
    graphics::plot(x$havokSS$t, x$Vr[,x$r], type = "l", xlab = "Time", ylab = "Forcing", ...)
  }

  if (what == "both") {
     graphics::par(mfrow=c(2,1), mai = c(0.5, 1.0, 0.1, 0.1))

     graphics::plot(x$havokSS$t, x$Vr[,1], type = "l", xlab = NA, ylab = "Value", ...)

     graphics::par(mai = c(1, 1.0, 0.1, 0.1))

     graphics::plot(x$havokSS$t, x$Vr[,x$r], type = "l", xlab = "Time", ylab = "Forcing", ...)

     graphics::par(mfrow=c(1,1), mai = c(1.02, 0.82, 0.82, 0.42))

  }

  if (what == "U-modes") {
    plotBuff <- (.1 * (max(x$U[, 1:x$r]) - min(x$U[, 1:x$r])))
    graphics::plot(x$U[,1], ylab = "Ur", xlab = "Time",
                   type = "l",
                   xaxt = "n",
                   col = grDevices::rainbow(x$r)[1],
                   ylim = c(min(x$U[,1:x$r]) - (.1*plotBuff), max(x$U[,1:x$r]) + (.1*plotBuff)))
    graphics::axis(1, at = seq(1, ncol(x$U), length.out = 10),
                   labels = round(seq(1, max(x$havokSS$t), length.out = 10)))

    for (i in 1:x$r){
      graphics::lines(x$U[,i], col = grDevices::rainbow(x$r, alpha = 1/sqrt(i))[i])
    }

    if (x$r < 6){
      graphics::legend("topleft",
             fill = grDevices::rainbow(x$r, alpha = 1/sqrt(1:x$r)),
             legend = paste("r = ", 1:x$r, sep = ""))
    }

    if (x$r >= 6){
      graphics::legend("topleft",
             fill = grDevices::rainbow(x$r, alpha = 1/sqrt(1:x$r))[c(1:6, x$r)],
             legend = c(paste("r = ", 1:5, sep = ""),
                        "...",
                        paste("r = ", x$r, sep = ""))
      )
    }

  }

  if (what == "embedded"){

    graphics::plot(x$Vr[,1], x$Vr[,2], type = "l", xlab = "V1", ylab = "V2", ...)

  }

  if (what == "nonlinear"){
    havForce <- active_forcing(x)
    graphics::plot(x$Vr[,1], x$Vr[,2], type = "l", xlab = "V1", ylab = "V2", ...)
    segments(head(x$Vr[,1], -1), head(x$Vr[,2], -1), x$Vr[,1][-1], x$Vr[,2][-1], ifelse(havForce$active==1,"red","black"))

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




