#' Sparse Identification of Nonlinear Dynamics (SINDy)
#'
#' @description Sparsification function based on sequential thresholded least-squares
#' as shown in the SINDy algorithm in "Discovering governing equations from data:
#' Sparse identification of nonlinear dynamical systems" (Brunton, Proctor, & Kutz, 2016).
#' @param x A vector or matrix of measurments over time.
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


sindy <- function(x, lambda, polyOrder = 5, useSine = FALSE, nVars = ncol(xdat)) {

  Theta <- pool_data(x, nVars = nVars, polyOrder = polyOrder, useSine = useSine)

  Xi <- sparsifyDynamics(Theta,dx,lambda,n)



}











