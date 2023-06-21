#' Matrix Sparsification for SINDy Algorithm
#'
#' @description Sparsification function based on sequential thresholded least-squares
#' as shown in the SINDy algorithm in "Discovering governing equations from data:
#' Sparse identification of nonlinear dynamical systems" (Brunton, Proctor, & Kutz, 2016).
#' @param Theta A matrix of candidate functions.
#' @param dXdt A matrix of first order derivatives of the variables of interest with respect to time.
#' @param lambda A numeric value; sparsification threshold.
#' @param loops ADD DESCRIPTION
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
par_sparsify_dynamics <- function(Theta, dXdt, lambda, loops){
  
  dXdt <- as.matrix(dXdt)
  termin <- T
  
  # Original regression result
  Xi <- pracma::mldivide(Theta, dXdt)
  
  if (lambda!=0) {
    dXdt <- as.data.frame(dXdt)
    
    reregress <- function (x,d) {
      biginds <- abs(x) > lambda
      d <- as.matrix(d)
      x[biginds] <- pracma::mldivide(Theta[,biginds],d)
      x
    }
    
    # lambda is our sparsification knob.
    for (k in 1:loops) {
      smallinds <- abs(Xi) < lambda    #find small coefficients
      Xi[smallinds] <- 0                # and threshold
      Xi <- as.data.frame(Xi)
      if (sum(sapply(Xi,function(x) all(x==0))>0)) {
        termin <- F 
        break}
      Xi <- mapply(reregress, Xi, dXdt)
    }
  }
  outp <- list(Xi,termin)
  return(outp)
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
