#' Matrix Sparsification for SINDy Algorithm
#'
#' @description Sparsification function based on sequential thresholded least-squares
#' as shown in the SINDy algorithm in "Discovering governing equations from data:
#' Sparse identification of nonlinear dynamical systems" (Brunton, Proctor, & Kutz, 2016).
#' @param Theta A matrix of candidate functions.
#' @param dXdt A matrix of first order derivatives of the variables of interest with respect to time.
#' @param lambda A numeric value; sparsification threshold.
#' @return  A matrix of sparse coefficients.
#' @references Brunton, S. L., Proctor, J. L., & Kutz, J. N. (2016). Discovering
#' governing equations from data by sparse identification of nonlinear dynamical
#' systems. Proceedings of the National Academy of Sciences, 113(15), 3932-3937.
#' @examples
#' \dontrun{
#' sparsify_dynamics(Theta, dXdt, lambda, n)
#' sparsify_dynamics(Theta, dXdt, 0.1, 10)
#' sparsify_dynamics(pool_data(yIn, 15, 5, TRUE), dXdt, 0, 15)
#' }
###################################
#' @export

svd_sign_flip <- function(x, U, S, V){



}


#X <- rmvnorm(500, mean = c(0,1,2,3,4,5,6,7,8,9,10))
#svd_X <- svd(X)





#   for m = 1:order % for each mode
# for f=1:F(m) % for each component
# s=[];
# a = loads{m}(:,f);
# a = a /(a'*a);
#             x = subtract_otherfactors_tucker(X, loads, m, f);
#             for i=1:size(x(:,:),2) % for each column
#                 s(i)=(a'*x(:,i));
# s(i)=sign(s(i))*power(s(i),2);
# end
# S(m,f) =sum(s);
# end
# end
# sgns = sign(S);
#
#
#
#
#
