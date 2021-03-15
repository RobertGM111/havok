#' Align singular vectors with majority of data vectors
#' @description Singular value decomposition is specified up to the sign (+/-) of the singular vectors.
#' This algorithm attempts to align the singular vectors with the majority of vectors of the input matrix.
#' @usage svd_align(x, r = NA)
#' @param x A data matrix.
#' @param r An integer; number of columns to include as a reduced rank SVD solution. Defaults to full rank.
#' @return \code{svd_flip} Aligned svd
#' @references Bro, R., Acar, E., & Kolda, T. G. (2008). Resolving the sign ambiguity in the singular value decomposition. Journal of Chemometrics: A Journal of the Chemometrics Society, 22(2), 135-140.
#' @examples
#' # Generate data
#' x <- cbind(runif(10, -0.5, 2), runif(10, -0.5, 2))
#'
#'
#' svd_flip <- svd_align(x)
#'
#' # Extract components of svd and plot
#' U <- svd_flip$u
#' Sigma <- svd_flip$d
#' VT <- t(svd_flip$v)
#'
#' # Plot data
#' plot(0, 0, xlim=c(-2,2), ylim=c(-2,2), type="n")
#'
#' for (i in 1:dim(x)[1]){
#'   arrows(0, 0, x[i,1], x[i,2], lwd=2)
#' }
#'
#' arrows(0, 0, VT[1, 1], VT[1, 2], col='red')
#'
#' @export
svd_align <- function(x, r = NA) {

  # Perform vanilla svd
  svd_flip <- svd(x)

  # Truncate results
  if (is.na(r)){
    U <- svd_flip$u
    Sigma <- svd_flip$d
    VT <- t(svd_flip$v)
    Sigma_mat <- diag(Sigma)
  } else {
    U <- svd_flip$u[ ,1:r]
    Sigma <- svd_flip$d[1:r]
    VT <- t(svd_flip$v[ ,1:r])
    Sigma_mat <- diag(Sigma)
  }

  Y <- matrix(NA, dim(x))

  # Svd signflip algorithm
  for(k in 1:dim(Sigma_mat)[1]){
    k = 1
    Sigma_mat_k = diag(Sigma)
    Sigma_mat_k[k,k] = 0

    # Correction for correlations
    Y <- x - U %*% Sigma_mat_k %*% VT

    # Calculate left and right signs
    UY <- t(U[,k]) %*% Y
    VTY <- VT[k,] %*% t(Y)
    sleft <- sign(UY)%*%t(UY**2)
    sright <- sign(VTY)%*%t(VTY**2)

    if(sleft*sright < 0){
      if(sleft < sright){
        sleft <- -sleft
      }
      else{
        sleft <- -sright
      }
    }

    # Update the signs for the kth left and right singular vector
    U[,k] <- U[,k] %*% sign(sleft)
    VT[k,] <- VT[k,] %*% sign(sright)
  }

  # Modify signs in svd object and return
  svd_flip$u <- U
  svd_flip$v <- t(VT)
  svd_flip$d <- Sigma

  return(svd_flip)
}


# Copyright 2021 Robert Glenn Moulder Jr., Elena Martynova, & Shashwat Kumar
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
