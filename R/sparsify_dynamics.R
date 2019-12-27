#' Matrix Sparcification for SINDy Algorithm
#'
#' @description Sparcification function for sparsifying dynamics as shown in
#' the SINDy algorithm as shown in "Discovering
#' governing equations from data: Sparse identification of nonlinear dynamical
#' systems" (Brunton, Proctor, & Kutz, 2015).
#' @param Theta A number.
#' @param dXdt A number.
#' @param lambda A number.
#' @param n A number.
#' @return  A matrix of sparse coefficients
#' @examples
#' add(1, 1)
#' add(10, 1)
###################################

sparsify_dynamics <- function(Theta, dXdt, lambda, n){

  # Original regression result
  Xi <- pracma::mldivide(Theta, dXdt)

  # lambda is our sparsification knob.
  for (k in 1:10) {
    smallinds <- abs(Xi) < lambda    #find small coefficients
    Xi[smallinds] <- 0                # and threshold
    for (ind in 1:n) {                   # n is state dimension
      biginds <- !smallinds[ , ind]
      if (n == 1) {
        Xi[biginds,ind] <- pracma::mldivide(Theta[ , biginds], dXdt)
      } else {
        Xi[biginds,ind] <- pracma::mldivide(Theta[ , biginds], dXdt[ , ind])
      }
    }
  }
  return(Xi)
}
