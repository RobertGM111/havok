#' Parallelized Hankel Alternative View of Koopman (pHAVOK) Analysis 
#' 
#' @description Parallel HAVOK (phavok) is a parallelized and optimized version of the HAVOK procedure. 
#' It estimates multiple HAVOK models simultaneously across multiple selected or 
#' randomized hyperparameter sets. phavok() is intended for model 
#' selection and inspection of the model fit surfaces. Once the model of interest is selected, 
#' it should be refit with the havok().
#' 
#' Note: If the model selected with phavok() has a sparsification parameter ‘lambda’ larger 
#' than 0, and the set of its hyperparameters does not result in an approximately equivalent 
#' model fit with havok(), a slight adjustment in the ‘lambda’ parameter when fitting the model 
#' with havok() will be required to achieve approximately equivalent model fit. To select appropriate 
#' adjustment magnitude we recommend inspecting nearby ‘lambda’ values in the ‘stackmax’ and ‘r’ 
#' combination that corresponds to the selected model.
#' 
#' @param xdat A vector of measurements over time.
#' @param dt A numeric value indicating the time-lag between two subsequent time series measurements.
#' @param stackmaxes A vector of 'stackmax' hyperparameter values, where 'stackmax' stands for the number of shift-stacked 
#' rows in the Hankel matrix.
#' @param rs A vector of 'r' hyperparameter values or NA, where 'r' stands for the number of singular vectors to include (also 
#' known as model degree or truncation parameter). If NA is selected, HAVOK models with all possible 'r' hyperparameters within
#' selected 'stackmaxes' will be fit.
#' @param random A numeric value from 0 to 1; what proportion of 'stackmaxes' and 'rs' combinations should be selected randomly? 
#' Both 0 and 1 result in fitting all possible models.
#' @param sparsify Logical; should models be sparsified?
#' @param sparserandom A numeric value from 0 to 1; what proportion of sparsification parameters should be selected randomly? 
#' @param loops An integer; number of times sequential thresholded least-squares procedure is repeated.
#' @param devMethod A character string; One of either \code{"FOCD"} for fourth order central difference or \code{"GLLA"} for generalized local linear approximation.
#' @param gllaEmbed An integer; the embedding dimension used for \code{devMethod = "GLLA"}.
#' @param alignSVD Logical; Whether the singular vectors should be aligned with the data.
#' @param numCores  An integer; number of cores to be used by phavok(). If not specified, defaults to the number of 
#' cores detected - 2. 
#' @return An dataframe with the following columns: \itemize{
#' \item{\code{stackmax} - }{'stackmax' hyperparameter.}
#' \item{\code{r} - }{'r' hyperparameter}
#' \item{\code{R2} - }{Squared correlation between the model predicted v_1 and v_1 exracted from SVD. A model fit estimate.}
#' \item{\code{Lambda} - }{Sparsification threshold.}
#' \item{\code{LambdaDeletions} - }{Number of model coefficient matrix elements that were expected to be truncated by the sparsification threshold.}
#' \item{\code{TrueLambdaDeletions} - }{Number of model coefficient matrix elements that were truncated by the sparsification threshold.}
#' \item{\code{kurtosis} - }{Pearson's measure of kurtosis of the forcing term value distribution.}
#' \item{\code{prop2sd} - }{Proportion of the forcing values exceeding the threshold of +-2 standard deviations.}
#' \item{\code{prop1.5sd} - }{Proportion of the forcing values exceeding the threshold of +-1.5 standard deviations.}
#' \item{\code{prop1sd} - }{Proportion of the forcing values exceeding the threshold of +-1 standard deviation.}}
#' @references S. L. Brunton, B. W. Brunton, J. L. Proctor, E. Kaiser, and J. N. Kutz,
#' "Chaos as an intermittently forced linear system," Nature Communications, 8(19):1-9, 2017.
#' @examples
#' \dontrun{
#' # Russian Twitter Troll Activity Example
#' 
#' library(plotly)
#' 
#' load("~/Downloads/results.all.trolls.RData")  # Contains time series of Russian Twitter 
#' # troll activity extracted 4 times per day during the US presidential  election year 2016 on 11 different topics 
#' right <- results.all.truncated[results.all.truncated$Type=="Right",] # only right-wing trolls
#' xdat <- right$Topic3  # Russian Twitter troll posting activity on the topic of Racial Justice/Black Lives Matter
#' dt <- 0.25   # 4 measurements per day
#' 
#' # All possible rs within specified stackmax range, no sparsification dimension
#' 
#' results <- phavok(xdat = xdat, dt = dt, stackmaxes = 28:58)
#' 
#' 
#' 
#' plot_ly(data = results, type = "scatter", x = ~stackmax,
#'         y = ~r, colors= "PiYG" , color = ~ R2, size = ~ R2,
#'         mode = "markers", text = ~R2,
#'         hovertemplate = paste('stackmax = %{x}',
#'                               '<br>r = %{y}<br>',
#'                               'R2 = %{text:.2f}',
#'                               '<extra></extra>')) %>%
#'   layout(title = 'R2') 
#' 
#' 
#' 
#' plot_ly(data = results, type = "scatter", x = ~stackmax,
#' y = ~r, colors= "BrBG" , color = ~ kurtosis, size = ~ kurtosis*(-1),
#' mode = "markers", text = ~kurtosis,
#' hovertemplate = paste('stackmax = %{x}',
#'                       '<br>r = %{y}<br>',
#'                       'kurtosis = %{text:.2f}',
#'                       '<extra></extra>')) %>%
#'   layout(title = 'kurtosis')
#' 
#' plot_ly(data = results, x = ~stackmax, y = ~r,
#'colors = "viridis" , color = ~prop2sd, size = ~prop2sd * (-1),
#'mode = "markers", text = ~prop2sd,
#'hovertemplate = paste('stackmax = %{x}',
#'                      '<br>r = %{y}<br>',
#'                      'prop2sd = %{text:.2f}',
#'                      '<extra></extra>')) %>%
#'  layout(title = 'prop2sd')
#' 
#' }
#' @export

phavok <- function(xdat, dt = 1, stackmaxes = NA, rs = NA, random = 0, 
                   sparsify = FALSE, sparserandom = 0, loops = 1,
                   devMethod = "FOCD", gllaEmbed = NA, alignSVD = TRUE,
                   numCores = parallel::detectCores(all.tests = FALSE, logical = TRUE)-2){
  
  if (random ==1) {
    random <- 0
  }
  
  if (random > 0) {
    
    if (is.na(rs[1])) {    
      rs <- 2:max(stackmaxes)
    } 
    
    srcombos <- matrix(0, nrow = max(rs), ncol = max(stackmaxes))
    
    for (sss in stackmaxes){
      for (rrr in rs) {
        if (sss > rrr-1) {
          srcombos[rrr,sss] <- 1 
        }
      }
    }
    
    # randomly select given the specified proportion
    
    howmany <- max((sum(srcombos==1)*random),1) 
    
    # randomly select howmany slots from srcombos = 1, randomly sample indices (rows of which array.ind)
    
    eligibleIndices <- which(srcombos==1, arr.ind=TRUE)
    
    selected <- sample(nrow(eligibleIndices), howmany, replace=F)
    
    stackmaxes <- unique(eligibleIndices[selected,][,2])
    
    rstackmax <- eligibleIndices[selected,]
    
  } 
  
  if (!is.na(rs[1]) & random ==0){
    rstackmax <- rs
  }  
  
  if (is.na(rs[1]) & random ==0) {
    rstackmax <- NA
  }
  
  stackmaxes <- sort(stackmaxes)
  samples_distributed <- vector(length=ceiling(length(stackmaxes)/numCores)*numCores)
  reverse <- FALSE
  batch_size <- ceiling(length(stackmaxes)/numCores)
  
  for(i in 1:batch_size){
    for (j in 1:numCores){
      original_index <- j + (i-1)*numCores
      if(reverse){
        target_index <- (numCores-j)*batch_size + i
        samples_distributed[target_index] <- stackmaxes[original_index]
      }
      else{
        target_index <- i + (j-1)*batch_size
        samples_distributed[target_index] <- stackmaxes[original_index]
      }
    }
    reverse <- !reverse
  }
  
  cl <- parallel::makeCluster(numCores)
  res <- parallel::parSapply(cl, samples_distributed, loop_havok, xdat,
                             dt, rstackmax, sparsify, sparserandom, loops,
                             devMethod, gllaEmbed, alignSVD)
  parallel::stopCluster(cl)
  
  res <- res[lapply(res,function(x) all(is.na(x))) == FALSE]
  
  res <- as.data.frame(do.call(rbind, res))
  
  class(res) <- c("data.frame","phavok")
  
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
