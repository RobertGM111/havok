#' Sparse Identification of Nonlinear Dynamics (SINDy)
#'
#' @description Sparsification function based on sequential thresholded least-squares
#' as shown in the SINDy algorithm in "Discovering governing equations from data:
#' Sparse identification of nonlinear dynamical systems" (Brunton, Proctor, & Kutz, 2016).
#' @param x A vector or matrix of measurements over time.
#' @param dt A numeric value indicating the time-lag between two subsequent time series measures.
#' @param lambda A numeric value; sparsification threshold.
#' @param polyOrder An integer from 0 to 5 indicating the highest degree of polynomials
#' included in the matrix of candidate functions.
#' @param useSine A logical value indicating whether sine and cosine functions
#' of variables should be added to the library of potential candidate functions.
#' If TRUE, candidate function matrix is augmented with sine and cosine functions
#' of integer multiples 1 through 10 of all the variables in \code{yIn}.
#' @param normalize Logical; Should the columns of the matrix of candidate functions be normalized?
#' @param nVars An integer that indicates the number of variables.
#' @return  A matrix of candidate functions and a matrix of sparse coefficients.
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
sindy <- function(x, dt, lambda, polyOrder = 5, useSine = FALSE,
                  normalize = FALSE, nVars = ncol(as.matrix(x))) {

  x <- as.matrix(x)

  dXdt <- compute_derivative(x = x, dt = dt)

  x <- x[3:(nrow(x) - 3), ]

  Theta <- pool_data(x, nVars = nVars, polyOrder = polyOrder, useSine = useSine)

  if (normalize == TRUE){
    normTheta <- rep(NA, dim(Theta)[2])

    for (k in 1:dim(Theta)[2]) {
      normTheta[k] <- sqrt(sum(Theta[ , k]^2))
      Theta[ , k] <- Theta[ , k]/normTheta[k]
    }
  }

  Xi <- sparsify_dynamics(Theta, dXdt = dXdt, lambda = lambda)

  res <- list("candidateFunctions" = Theta,
              "sparse" = Xi)

  class(res) <- "sindy"
  return(res)
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





