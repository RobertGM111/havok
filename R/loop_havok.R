#' Hankel Alternative View of Koopman (HAVOK) Analysis
#'
#' @description Data-driven decomposition of chaotic time series into an intermittently
#' forced linear system. HAVOK combines delay embedding and Koopman theory to decompose
#' chaotic dynamics into a linear model in the leading delay coordinates with forcing by
#' low-energy delay coordinates. Forcing activity demarcates coherent phase space regions
#' where the dynamics are approximately linear from those that are strongly nonlinear.
#' @param xdat A vector of equally spaced measurements over time.
#' @param dt A numeric value indicating the time-lag between two subsequent time series measurements.
#' @param stackmax An integer; number of shift-stacked rows.
#' @param sparsify ADD
#' @param sparserandom ADD
#' @param loops ADD
#' @param center Logical; should \code{xdat} be centered around 0?
#' @param rstackmax ADD
#' @param polyOrder NEED CLEANED
#' @param useSine NEED CLEANED
#' @param discrete Logical; is the underlying system discrete?
#' @param devMethod A character string; One of either \code{"FOCD"} for fourth order central difference or \code{"GLLA"} for generalized local linear approximation.
#' @param gllaEmbed An integer; the embedding dimension used for \code{devMethod = "GLLA"}.
#' @param alignSVD Logical; Whether the singular vectors should be aligned with the data.
#' @return An object of class 'havok' with the following components: \itemize{
#' \item{\code{havokSS} - }{A HAVOK analysis generated state space model with its time history.}
#' \item{\code{params} - }{A matrix of parameter values used for this function.}
#' \item{\code{dVrdt} - }{A matrix of first order derivatives of the reduced rank V matrix with respect to time.}
#' \item{\code{r} - }{Estimated optimal number singular vectors to include in analysis.}
#' \item{\code{sys} - }{HAVOK model represented in state-space form.}
#' \item{\code{normTheta} - }{Normalized matrix of candidate functions obtained from \code{\link{pool_data}}.}
#' \item{\code{Xi} - }{A matrix of sparse coefficients obtained from \code{\link{sparsify_dynamics}}.}
#' \item{\code{Vr} - }{The reduced rank V matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{Ur} - }{The reduced rank U matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{sigsr} - }{Values of the diagonal of the reduced rank \eqn{\Sigma} matrix of the SVD of the Hankel matrix of the time series.}
#' \item{\code{Vr_aligned} - }{Vr truncated based upon \code{devMethod}}}
#' @references S. L. Brunton, B. W. Brunton, J. L. Proctor, E. Kaiser, and J. N. Kutz,
#' "Chaos as an intermittently forced linear system," Nature Communications, 8(19):1-9, 2017.
#' @examples
#' \dontrun{
#'#Lorenz Attractor
#'#Generate Data
#'##Set Lorenz Parameters
#'parameters <- c(s = 10, r = 28, b = 8/3)
#'n <- 3
#'state <- c(X = -8, Y = 8, Z =27) ##Inital Values
#'
#'#Intergrate
#'dt <- 0.001
#'tspan <- seq(dt, 200, dt)
#'N <- length(tspan)
#'
#'Lorenz <- function(t, state, parameters) {
#'  with(as.list(c(state, parameters)), {
#'    dX <- s * (Y - X)
#'    dY <- X * (r - Z) - Y
#'    dZ <- X * Y - b * Z
#'    list(c(dX, dY, dZ))
#'  })
#'}
#'
#'out <- ode(y = state, times = tspan, func = Lorenz, parms = parameters, rtol = 1e-12, atol = 1e-12)
#'xdat <- out[, "X"]
#'
#'hav <- havok(xdat = xdat, dt = dt, stackmax = 100, lambda = 0,
#'             rmax = 15, polyOrder = 1, useSine = FALSE)
#'
#'# ECG Example
#'
#'data(ECG_measurements)
#'
#'xdat <- ECG_measurements[,"channel1"]
#'dt <- ECG_measurements[2,"time"] - ECG_measurements[1,"time"]
#'stackmax <- 25
#'rmax <- 5
#'lambda <- .001
#'hav <- havok(xdat = xdat, dt = dt, stackmax = stackmax, lambda = lambda,
#'             rmax = 5, polyOrder = 1, useSine = FALSE)
#'plot(hav)
#'}
###################################

#' @export
loop_havok <- function(stackmax, xdat = xdat, dt = dt,
                       rstackmax = rstackmax, sparsify = sparsify, 
                       sparserandom = sparserandom, loops = loops, 
                       devMethod = devMethod, gllaEmbed = gllaEmbed, 
                       alignSVD = alignSVD, center = TRUE, polyOrder = 1,
                       useSine = FALSE, discrete = FALSE) {
  
  ################################################### 
  library(havok)
  library(deSolve)
  library(pracma)
  library(control)
  library(rsvd)
  library(moments)
  
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
  
  rsvd_align <- function(x, r = NA) {
    
    # Perform rsvd
    
    svd_flip <- rsvd::rsvd(x,k=r,nu=r, nv=r, p=10,q=2,sdist="normal")
    U <- svd_flip$u
    Sigma <- svd_flip$d
    VT <- t(svd_flip$v)
    Sigma_mat <- diag(Sigma)
    
    
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
  
  ##############################################################################
  
  if (is.na(stackmax)) {
    if(exists("res")){
      res2 <- cbind(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      colnames(res2) <- c("stackmax", "rr", "R2", "lambda", "lambdaDeletions", "CascadeLambdaDeletions", "kurtosis", "prop2sd", "prop1.5sd", "prop1sd")
      res <- rbind(res, res2)
    }else{
      return(NA)
    }
  } else {
    
    if (center == TRUE){
      xdat <- xdat - mean(xdat)
    }
    
    H <- build_hankel(x = xdat, nrows = stackmax)
    
    if (class(rstackmax)[1]=="logical") {rr_samples <- stackmax:2 }
    if (class(rstackmax)[1]=="numeric" | class(rstackmax)[1]=="integer") {rr_samples <- min(max(rstackmax),stackmax):max(min(rstackmax),2)}
    if (class(rstackmax)[1]=="matrix") {rr_samples <- rstackmax[,1][rstackmax[,2]==stackmax]}
    
    
    rr_samples <- sort(rr_samples,decreasing=TRUE)
    rset <- max(rr_samples)
    
    # Rsvd, align and reduce rank 
    if (alignSVD == TRUE){
      USV <- rsvd_align(H, r = stackmax)
      U <- USV$u
      sigs <- USV$d
      V <- USV$v
    } else {
      USV <- rsvd(H,k=stackmax,nu=stackmax, nv=stackmax, p=10,q=2,sdist="normal")
      U <- USV$u
      sigs <- USV$d
      V <- USV$v
    }
    
    
    if (!is.na(rset)) {
      r <- rset
    }
    
    
    if (discrete == FALSE){
      
      if (devMethod == "GLLA"){
        
        dV <- matrix(NA, nrow = nrow(V) - (gllaEmbed - 1), ncol = r)
        devList <- compute_derivative(V, dt, devMethod = "GLLA",
                                      gllaEmbed = gllaEmbed,
                                      gllaTau = 1,
                                      gllaOrder = 1)
        for (i in 1:r){
          dV[,i] <- devList[[i]][,2]
        }
        
        x <- V[ceiling(gllaEmbed/2):(nrow(V) - floor(gllaEmbed/2)), 1:r]
        dx <- dV
        
      }
      
      if (devMethod == "FOCD"){
        dV <- compute_derivative(x = V, dt = dt, r = r, devMethod = "FOCD")
        x <- V[3:(nrow(V) - 3), 1:r]
        dx <- dV
      }
      
      for (rr in rr_samples) { 
        
        lambda <- 0
        
        U <- USV$u[,1:rr]
        sigs <- USV$d[1:rr]
        V <- USV$v[,1:rr]
        x <- x[,1:rr]
        dx <- dx[,1:rr]
        r <- rr
        
        Theta <- pool_data(x, r, polyOrder = polyOrder, useSine)
        
        normTheta <- rep(NA, dim(Theta)[2])
        
        for (k in 1:dim(Theta)[2]) {
          normTheta[k] <- sqrt(sum(Theta[ , k]^2))
          Theta[ , k] <- Theta[ , k]/normTheta[k]
        }
        
        m <- dim(Theta)[2]
        
        sparse <- par_sparsify_dynamics(Theta,dx,lambda, loops)
        
        Xi <- sparse[[1]]
        
        for (k in 1:max(dim(Xi))) {
          Xi[k, ] <- Xi[k, ] / normTheta[k]
        }
        A <- t(Xi[2:(r + 1), 1:(r - 1)])
        B <- A[, r]
        A <- A[ , 1:(r - 1)]
        L <- 1:nrow(x)
        
        sys <- control::ss(A, B, pracma::eye(r - 1), 0 * B)
        HAVOK <- control::lsim(sys, x[L, r], dt * (L - 1), x[1, 1:(r - 1)])
        
        
        R2 <- cor(x[,1],HAVOK$y[1,])^2
        kurtosis <- moments::kurtosis(x[,r])
        prop2sd <- mean((x[,r]>(mean(x[,r])+2*sd(x[,r])))|(x[,r]<(mean(x[,r])-2*sd(x[,r]))))
        prop1.5sd <- mean((x[,r]>(mean(x[,r])+1.5*sd(x[,r])))|(x[,r]<(mean(x[,r])-1.5*sd(x[,r]))))
        prop1sd <- mean((x[,r]>(mean(x[,r])+1*sd(x[,r])))|(x[,r]<(mean(x[,r])-1*sd(x[,r]))))
        
        lambdaDeletions <- 0
        CascadeLambdaDeletions <- 0
        
        if(rr == rr_samples[1]){
          res <- cbind(stackmax, rr, R2, lambda, lambdaDeletions, CascadeLambdaDeletions,kurtosis, prop2sd, prop1.5sd, prop1sd)
          colnames(res) <- c("stackmax", "r", "R2", "Lambda", "LambdaDeletions", "TrueLambdaDeletions", "kurtosis", "prop2sd", "prop1.5sd", "prop1sd")
        } else {
          res2 <- cbind(stackmax, rr, R2, lambda, lambdaDeletions, CascadeLambdaDeletions,kurtosis, prop2sd,prop1.5sd,prop1sd)
          colnames(res2) <- c("stackmax", "r", "R2", "Lambda", "LambdaDeletions", "TrueLambdaDeletions", "kurtosis", "prop2sd", "prop1.5sd", "prop1sd")
          res <- rbind(res, res2)
        }
        
        if (sparsify==T & r > 2){
          
          orderedCoefs <- sort(abs(c(as.vector(A), as.vector(B))))
          
          newLambdas <- (orderedCoefs[2:((r-1)*r - (r-2)*2)] + orderedCoefs[1:((r-1)*r - (r-2)*2 - 1)])/2
          
          # How many of the newLambdas to select?
          
          if (sparserandom > 0) {
            
            howmany <- max((length(newLambdas)*sparserandom),1)
            
            sampledLambdas <- sample(newLambdas,howmany, replace = F)
            sampledLambdas <- sort(sampledLambdas, decreasing = F)
            
          }
          
          if (sparserandom == 0) {
            sampledLambdas <- newLambdas
          }
          
          for (lambda in sampledLambdas){ 
            
            lambdaDeletions <- which(newLambdas==lambda)
            
            sparse <- par_sparsify_dynamics(Theta,dx,lambda,loops)
            
            if (sparse[[2]]==F) {break}
            Xi <- sparse[[1]]
            
            for (k in 1:max(dim(Xi))) {
              Xi[k, ] <- Xi[k, ] / normTheta[k]
            }
            A <- t(Xi[2:(r + 1), 1:(r - 1)])
            B <- A[, r]
            A <- A[ , 1:(r - 1)]
            L <- 1:nrow(x)
            
            AB <- cbind(A, B)
            CascadeLambdaDeletions <- sum(AB==0)
            
            sys <- control::ss(A, B, pracma::eye(r - 1), 0 * B)
            HAVOK <- control::lsim(sys, x[L, r], dt * (L - 1), x[1, 1:(r - 1)])
            
            
            R2 <- cor(x[,1],HAVOK$y[1,])^2
            kurtosis <- kurtosis(x[,r])
            prop2sd <- mean((x[,r]>(mean(x[,r])+2*sd(x[,r])))|(x[,r]<(mean(x[,r])-2*sd(x[,r]))))
            prop1.5sd <- mean((x[,r]>(mean(x[,r])+1.5*sd(x[,r])))|(x[,r]<(mean(x[,r])-1.5*sd(x[,r]))))
            prop1sd <- mean((x[,r]>(mean(x[,r])+1*sd(x[,r])))|(x[,r]<(mean(x[,r])-1*sd(x[,r]))))
            
            
            res2 <- cbind(stackmax, rr, R2, lambda, lambdaDeletions, CascadeLambdaDeletions,kurtosis, prop2sd,prop1.5sd,prop1sd)
            colnames(res2) <- c("stackmax", "r", "R2", "Lambda", "LambdaDeletions", "TrueLambdaDeletions", "kurtosis", "prop2sd", "prop1.5sd", "prop1sd")
            res <- rbind(res, res2)
            
            
          }
        }
      }
    }
  }
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
