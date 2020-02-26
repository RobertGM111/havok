#' Determine when forcing is active for an object of class "havok"
#'
#' @description This function uses a minimum value to determin when forcing
#' is active in a fitted "havok" object.
#' @param x An object of class "havok"
#' @param minVal A cutoff value for determining when forcing is active. Defualts
#' to one standard deviation of the forcing term.
#' @return  A matrix of sparse coefficients.
#' @examples
#' \donttest{
#' hav <- havok(xdat = xdat, dt = dt)
#' active_forcing(x = hav)
#' }
###################################
#' @export


active_forcing <- function(x, minVal = stats::sd(x$x[,x$r])){
  if (class(x) != "havok"){
    stop("Object x must be of class \"havok\"")
  }

  allForcing <- x$x[,x$r]

  forcingOn <- ifelse(abs(allForcing) >= minVal, 1, 0)

  res <- list("forcing" = allForcing, "active" = forcingOn)

  return(res)
}





