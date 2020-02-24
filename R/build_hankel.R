#' Build a Hankel matrix from a time series.
#'
#' @description Create a Hankel matrix from a vector of of measurements over time.
#' @param x A vector or matrix of measurements over time.
#' @param stackmax An integer; the number of shift-stacked rows.
#' @return A time-lagged matrix of stacked measurements over time.
#' @examples
#' \donttest{
#' build_hankel(x = xdat, stackmax = 15)
#' }
###################################
#' @export
build_hankel <- function(x, stackmax){

  H <- matrix(0, nrow = stackmax, ncol = length(x) - stackmax)

  for (k in 1:stackmax) {
    H[k,] <- x[k:(length(x) - stackmax - 1 + k)]
  }

  return(H)

}
