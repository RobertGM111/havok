#' Plotting functions for the havok package of class “havok"
#'
#' @description Generic plotting function for object of class ("havok")
#' @param x A "havok" object.
#' @param ... Other calls to plot
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
#'t <- out[, "time"]
#'hav <- havok(xdat = xdat, dt = dt, stackmax = 100, lambda = 0,
#'             rmax = 15, polyOrder = 1, useSine = FALSE)
#'
#'plot(hav)
#'}
###################################

## S3 method for class "havok"
#' @export
plot.havok <- function(x,...){

  graphics::par(mfrow=c(2,1), mai = c(0.5, 1.0, 0.1, 0.1))

  graphics::plot(x$havok$t, x$havok$x[1,], type = "l", xlab = NA, ylab = "Value")

  graphics::par(mai = c(1, 1.0, 0.1, 0.1))

  graphics::plot(x$havok$t, x$x[,x$r], type = "l", xlab = "Time", ylab = "Forcing")

}








