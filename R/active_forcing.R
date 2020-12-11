#' Determine when forcing is active for an object of class \code{havok}.
#'
#' @description This function uses a threshold to determine when forcing
#' is active in a fitted \code{havok} object.
#' @usage active_forcing(x, thresh = stats::sd(x$Vr[, x$r]))
#' @param x An object of class \code{havok}.
#' @param thresh A cutoff value for determining when forcing is active. Defaults
#' to one standard deviation of the forcing term.
#' @return  A list of forcing values with their corresponding activity status.
#' @examples
#'data(ECG_measurements)
#'
#'xdat <- ECG_measurements[,"channel1"]
#'dt <- ECG_measurements[2,"time"] - ECG_measurements[1,"time"]
#'
#'stackmax <- 25
#'rmax <- 5
#'lambda <- .001
#'
#'hav <- havok(xdat = xdat, dt = dt, stackmax = stackmax, lambda = lambda, rmax = rmax)
#'active_forcing(x = hav)
###################################
#' @export

active_forcing <- function(x, thresh = stats::sd(x$Vr[,x$r])){
  if (class(x) != "havok"){
    stop("Object x must be of class \"havok\"")
  }

  allForcing <- x$Vr[,x$r]

  forcingOn <- ifelse(abs(allForcing) >= thresh, 1, 0)

  res <- list("forcing" = allForcing, "active" = forcingOn)

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
