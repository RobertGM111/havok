#' Matrix Sparcification for SINDy Algorithm
#'
#' @description Coefficient determining optimal location of Hard Threshold for Matrix
#'  Denoising by Singular Values Hard Thresholding when noise level is known or
#'  unknown.  Recreation of matlab code by Matan Gavish and David Donoho.
#' @param beta A number.
#' @param sigma_known A number.
#' @return  A matrix of sparse coefficients
#' @examples
#' add(1, 1)
#' add(10, 1)
###################################

optimal_SVHT_coef <- function(beta, sigma_known = FALSE) {
  if (sigma_known == TRUE) {
    coef = optimal_SVHT_coef_sigma_known(beta)
  } else {
    coef = optimal_SVHT_coef_sigma_unknown(beta)
  }
}

optimal_SVHT_coef_sigma_unknown <- function(beta) {
  coef.unknown <- optimal_SVHT_coef_sigma_known(beta)
  MPmedian <- rep(0, length(beta))
  for (i in 1:length(beta)) {
    MPmedian[i] <- median_marcenko_pastur(beta[i])
  }
  omega <- coef.unknown / sqrt(MPmedian)
  return(omega)
}

optimal_SVHT_coef_sigma_known <- function(beta) {
  w <- (8 * beta) / (beta + 1 + sqrt(beta^2 + 14 * beta + 1))
  lambda_star <- sqrt(2 * (beta + 1) + w)
  return(lambda_star)
}

inc_mar_pas <- function(x0, beta, gamma) {
  if (beta > 1) stop("Beta must be <= 1")
  topSpec <- (1 + sqrt(beta))^2
  botSpec <- (1 - sqrt(beta))^2
  MarPas <- function(x) ifelse((topSpec - x) * (x - botSpec) > 0,
                             sqrt((topSpec - x) * (x - botSpec)) /
                               (beta * x) / (2 * pi), 0)
  if (gamma != 0) {
    fun <- function(x) (x^gamma * MarPas(x))
  } else {
    fun <- function(x) MarPas(x)
  }
  Int <- integrate(fun, x0, topSpec)
  return(Int$value)
}

median_marcenko_pastur <- function(beta) {
  MarPas <- function(x) 1 - inc_mar_pas(x, beta, 0)
  lobnd <- (1 - sqrt(beta))^2
  hibnd <- (1 + sqrt(beta))^2
  change <- 1
  while (change & (hibnd - lobnd > .001)) {
    change <- 0
    x <- seq(lobnd, hibnd, length.out = 5)
    y <- rep(NA,length(x))
    for (i in 1:length(x)) {
      y[i] <- MarPas(x[i])
    }
    if (any(y < 0.5)) {
      lobnd = max(x[which(y < 0.5)])
      change = 1
    }
    if (any(y > 0.5)){
      hibnd = min(x[which(y > 0.5)])
      change = 1
    }
  }
  med = (hibnd + lobnd) / 2
  return(med)
}

marcenko_pastur_integral <- function(x, beta) {
  if (beta <= 0 | beta > 1) stop("beta must be > 0 and <=1")
  lobnd <- (1 - sqrt(beta))^2
  hibnd <- (1 + sqrt(beta))^2
  if (x < lobnd | x > hibnd) stop("X is out of bounds")
  dens <- function(t) sqrt((hibnd - t) * (t - lobnd)) / (2 * pi * beta * t)
  Int <- integrate(dens, lobnd, x)
  print(cat("\r", paste("x=", x, " beta=", beta,
                        " I=", round(Int$value, 5), sep = "")))
  return(Int$value)
}









