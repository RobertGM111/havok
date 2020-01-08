#' Data Pooling into the Matrix of Candidate Functions for SINDy Algorithm
#'
#' @description Pooling function for contruction of the library of potential right
#' hand-side candidate functions as shown in the SINDy algorithm in "Discovering
#' governing equations from data: Sparse identification of nonlinear dynamical
#' systems" (Brunton, Proctor, & Kutz, 2016).
#' @param yIn A matrix/data frame of time-dependent variables.
#' @param nVars An integer that indicates the number of variables.
#' @param polyOrder An integer from 0 to 5 indicating the highest degree of polynomials
#' included in the matrix of candidate functions.
#' @param useSine A logical value indicating whether sine and cosine functions
#' of variables should be added to the library of potential candidate functions.
#' If TRUE, candidate function matrix is augmented with sine and cosine functions
#' of integer multiples 1 through 10 of all the variables in \code{yIn}.
#' @return  A matrix of candidate functions.
#' @references Brunton, S. L., Proctor, J. L., & Kutz, J. N. (2016). Discovering
#' governing equations from data by sparse identification of nonlinear dynamical
#' systems. Proceedings of the National Academy of Sciences, 113(15), 3932-3937.
#' @examples
#' pool_data(yIn, nVars, polyOrder, useSine)
#' pool_data(yIn, 15, 5, TRUE)
#' pool_data(yIn, 3, 1, 0)
###################################
pool_data <- function(yIn, nVars, polyOrder, useSine) {
  n <- dim(yIn)[1]

  ind <- 1

  # poly order 0
  yOut <- matrix(NA, nrow = n,
                 ncol = 1 + nVars + (nVars * (nVars + 1) / 2) +
                   (nVars * (nVars + 1) * (nVars + 2) / (2 * 3)) + 11)

  yOut[ , ind] <- 1
  ind <- ind + 1

  # poly order 1
  for (i in 1:nVars) {
    yOut[ , ind] <- yIn[ , i]
    ind <- ind + 1
  }

  if (polyOrder >= 2) {
    # poly order 2
    for (i in 1:nVars) {
      for (j in 1:nVars) {
        yOut[ , ind] <- yIn[ , i] * yIn[ , j]
        ind <- ind + 1
      }
    }
  }

  if (polyOrder >= 3) {
    # poly order 3
    for (i in 1:nVars) {
      for (j in i:nVars) {
        for (k in j:nVars) {
          yOut[ , ind] <- yIn[ , i] * yIn[ , j] * yIn[ , k]
          ind <- ind + 1
        }
      }
    }
  }

  if (polyOrder >= 4) {
    # poly order 4
    for (i in 1:nVars) {
      for (j in i:nVars) {
        for (k in j:nVars) {
          for (l in k:nVars) {
            yOut[ , ind] <- yIn[ , i] * yIn[ , j] * yIn[ , k] * yIn[ , l]
            ind <- ind + 1
          }
        }
      }
    }
  }

  if (polyOrder >= 5) {
    # poly order 5
    for (i in 1:nVars) {
      for (j in i:nVars) {
        for (k in j:nVars) {
          for (l in k:nVars) {
            for (m in l:nVars) {
              yOut[ , ind] <- yIn[ , i] * yIn[ , j] * yIn[, k] * yIn[ , l] * yIn[ , m]
              ind <- ind + 1
            }
          }
        }
      }
    }
  }

  if (useSine == TRUE) {
    for(k in 1:10){
      yOut <- cbind(yOut, sin(k * yIn), cos(k * yIn))
    }
  }
  yOut <- yOut[, apply(yOut, 2, function(x) !any(is.na(x)))]
  return(yOut)
}
