#' Optimal Hard Threshold for Singular Values
#'
#' @description A function for the calculation of the coefficient determining optimal
#' location of hard threshold for matrix denoising by singular values hard thresholding
#' when noise level is known or unknown. Recreation of MATLAB code by Matan Gavish and
#' David Donoho.
#' @param beta A single value or a vector that represents aspect ratio m/n of the matrix to
#' be denoised. 0<\code{beta}<=1.
#' @param sigma_known A logical value. TRUE if noise level known, FALSE if unknown.
#' @return  Optimal location of hard threshold, up the median data singular value (sigma
#' unknown) or up to \code{sigma*sqrt(n)} (sigma known); a vector of the same dimension
#' as \code{beta}, where \code{coef[i]} is the coefficient corresponding to \code{beta[i]}.
#' @references Gavish, M., & Donoho, D. L. (2014). The optimal hard threshold for singular
#' values is 4/sqrt(3). IEEE Transactions on Information Theory, 60(8), 5040-5053.
#' @examples
#' \donttest{
#' # Usage in known noise level:
#' # Given an m-by-n matrix \code{Y} known to be low rank and observed in white noise
#' # with mean zero and known variance \code{sigma^2}, form a denoised matrix \code{Xhat} by:
#'  USV <- svd(Y)
#'  y <- USV$d
#'  y(y < (optimal_SVHT_coef(m/n,1) * sqrt(n) * sigma) ) <- 0
#'  Xhat <- USV$u * diag(y) * t(USV$v)
#'
#'  # Given an m-by-n matrix \code{Y} known to be low rank and observed in white
#'  # noise with mean zero and unknown variance, form a denoised matrix \code{Xhat} by:
#'  USV <- svd(Y)
#'  y <- USV$d
#'  y(y < (optimal_SVHT_coef(m/n,0) * median(y)) ) <- 0
#'  Xhat <- USV$u * diag(y) * t(USV$v)
#'  }
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
  Int <- stats::integrate(fun, x0, topSpec)
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
  Int <- stats::integrate(dens, lobnd, x)
  print(cat("\r", paste("x=", x, " beta=", beta,
                        " I=", round(Int$value, 5), sep = "")))
  return(Int$value)
}









