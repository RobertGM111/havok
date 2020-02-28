#' Determine when forcing is active for an object of class "havok"
#'
#' @description This function uses a minimum value to determin when forcing
#' is active in a fitted "havok" object.
#' @param x An object of class "havok"
#' @param minVal A cutoff value for determining when forcing is active. Defualts
#' to one standard deviation of the forcing term.
#' @return  A matrix of sparse coefficients.
#' @examples
#' \donttest{
#' hav <- havok(xdat = xdat, dt = dt)
#' active_forcing(hav)
#' }
###################################
#' @export


active_forcing <- function(x, minVal = stats::sd(x$Vr[,x$r])){
  if (class(x) != "havok"){
    stop("Object x must be of class \"havok\"")
  }

  allForcing <- x$Vr[,x$r]

  forcingOn <- ifelse(abs(allForcing) >= minVal, 1, 0)

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

