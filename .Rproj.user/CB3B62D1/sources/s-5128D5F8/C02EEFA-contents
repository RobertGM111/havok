#' Data Pooling for SINDy Algorithm
#'
#' @description Pooling function for sparsifyIng dynamics as shown in
#' the SINDy algorithm as shown in "Discovering
#' governing equations from data: Sparse identification of nonlinear dynamical
#' systems" (Brunton, Proctor, & Kutz, 2015).
#' @param yIn A number.
#' @param nVars A number.
#' @param polyOrder A number.
#' @param useSine A number.
#' @return  A matrix of candidate functions
#' @examples
#' add(1, 1)
#' add(10, 1)
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
