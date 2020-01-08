#' Matrix Sparsifcation for SINDy Algorithm
#'
#' @description Sparsification function based on sequential thresholded least-squares
#' as shown in the SINDy algorithm in "Discovering governing equations from data:
#' Sparse identification of nonlinear dynamical systems" (Brunton, Proctor, & Kutz, 2016).
#' @param Theta A matrix of candidate functions.
#' @param dXdt A matrix of first order derivatives of the variables of interest with respect to time.
#' @param lambda A numeric value; sparsification threshold.
#' @param n An integer that indicates the number of variables.
#' @return  A matrix of sparse coefficients.
#' @references Brunton, S. L., Proctor, J. L., & Kutz, J. N. (2016). Discovering
#' governing equations from data by sparse identification of nonlinear dynamical
#' systems. Proceedings of the National Academy of Sciences, 113(15), 3932-3937.
#' @examples
#' sparsify_dynamics(Theta, dXdt, lambda, n)
#' sparsify_dynamics(Theta, dXdt, 0.1, 10)
#' sparsify_dynamics(pool_data(yIn, 15, 5, TRUE), dXdt, 0, 15)
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
