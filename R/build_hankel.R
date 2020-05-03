#' Build a Hankel matrix from a time series.
#'
#' @description Create a Hankel matrix from a vector of measurements over time.
#' @param x A vector of measurements over time.
#' @param stackmax An integer; the number of shift-stacked rows.
#' @return A Hankel matrix of \code{x}.
#' @examples
#' \donttest{
#' build_hankel(x = xdat, stackmax = 15)
#' }
###################################
#' @export
build_hankel <- function(x, stackmax){

  H <- matrix(NA, nrow = stackmax, ncol = length(x) - stackmax)

  for (k in 1:stackmax) {
    H[k,] <- x[k:(length(x) - stackmax + k)]
  }

  return(H)

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
