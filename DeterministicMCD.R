# --------------------------------------------------------------------
# Author:

# Group 24
# Thijs de Bruin  -
# Colin Huliselan -
# Sanne Lin - 426416
# Mate Varadi - 495556
#
# --------------------------------------------------------------------

## Use this code skeleton to implement the deterministic MCD and the plug-in
## robust regression estimator.  Please submit your implementation together
## with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each group's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.

## Preparations
library(MASS)
library(robustbase)
#library(psych)

## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators


# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
  Y=apply(z,2,tanh)
  S1=cor(Y)
  return(S1)
}

# spearman correlation matrix
corSpearman <- function(z) {
  ranks<- apply(z,2,rank)
  S2=cor(ranks,method="spearman")
  return(S2)
}

# correlation matrix based on normal scores of the ranks
corNSR <- function(z) {
  ranks <- apply(z,2,rank)
  n=nrow(z)
  normalscores <- qnorm((ranks - 1/3)/(n + 1/3))
  S3=cor(normalscores)
  return(S3)
}

# modified spatial sign covariance matrix
covMSS <- function(z) {
  norms=apply(z,1, function(x) sqrt(sum(x^2)))
  k=z/norms
  
  #Avoid Nan when summing
  k[which(norms==0),]<-0
  
  n=nrow(z)
  #n=length(z)
  S4=(t(k)%*%k)/n
  return(S4)
}

# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  norms <- apply(z,1,function(x) sqrt(sum(x^2)))
  #Subset of the rank of the norms always smaller than or equal to the (rounded up) half of the obs
  subset <- z[rank(norms)<=ceiling(0.5*nrow(z)),]
  S5=cov(subset)
  return(S5)
}

# raw OGK estimator of the covariance matrix with median and Qn
rawCovOGK <- function(z) {
  
  output <- covOGK(z, sigmamu = s_Qn)
  S6 <- output$cov
  return(S6)
  # Hint: have a look at function covOGK() in package robustbase
}

## ___________________________________________________________________

# compare & run
# example data matix
set.seed(0)
y = rnorm(5000, 100, 40)
x = mvrnorm(n=5000,mu=c(0,2),Sigma=matrix(c(1,0.6,0.6,1),2,2))
x2= mvrnorm(n=250,mu=c(10,2),Sigma=matrix(c(1,0.6,0.6,1),2,2))
# sigma <- matrix(c(1,0.6,-0.4, 0.6, 1, 0.3, -0.4, 0.3, 1),3,3)
# x = mvrnorm(n=5000,mu=c(0,2,3),Sigma=sigma)
# x2= mvrnorm(n=250,mu=c(10,2,3),Sigma=sigma)
x <- rbind(x, x2)

# resultX <- covDetMCD(x, alpha=alpha)$raw.cov
# resultX_package <- DetMCD::DetMCD(x, alpha=alpha)$raw.cov
# resultX
# resultX_package
# 
# resultX_rw <- covDetMCD(x, alpha=alpha)$cov
# resultX_package_rw <- DetMCD::DetMCD(x, alpha=alpha)$cov
# resultX_rw
# resultX_package_rw

# standardize with median and Qn
z=apply(x,2,function(y) (y-median(y))/Q_n(y))


#Tools---------
#' Mahanalobis distance calculation
#'
#' @param X input matrix
#' @param t location parameter
#' @param s scale parameter
Mdistance <- function(X, t, s) {
  X <- t(apply(X, 1, function(x) x-t(t)))
  return(apply(X, 1, function(x) sqrt(t(x)%*%solve(s)%*%x)))
}

Q_n <- function(x){
  constant = 2.21914
  h <- round(length(x)/2)+1
  k <- choose(h,2)
  distance <- dist(x)
  return(constant*distance[order(distance)][k])
}

## ______________________________________________________________
## Main function for deterministic MCD algorithm

# Input:
# x ....... data matrix
# alpha ... proportion of observations to be used for subset size
# anything else you need

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# any other output you want to return

covDetMCD <- function(x, alpha, ...) {
  X_og <- x
  
  #standardize input
  x <- apply(x,2,function(y) (y-median(y))/Q_n(y))
  
  # Compute subset size
  h <- h.alpha.n(alpha, nrow(x), ncol(x))
  
  #Calculate initial estimates of the scaling matrixes
  initialcovs <- list()
  initialcovs[["S1"]] <-corHT(x)
  initialcovs[["S2"]] <-corNSR(x)
  initialcovs[["S3"]] <-corSpearman(x)
  initialcovs[["S4"]] <-covMSS(x)
  initialcovs[["S5"]] <-covBACON1(x)
  initialcovs[["S6"]] <-rawCovOGK(x)
  
  #prepare a temporary list for the calculation of MDCraw
  MDCRawtemp <- list()
  
  #For each of the 6 estimates do the following
  for (S in initialcovs){
    #3.1.(1) Compute the matrix E of eigenvectors of Sk
    E <- eigen(S)$vectors
    B <- x%*%E
    
    #3.1.(2) Estimate the covariance of Z
    L <- diag(apply(B,2,Q_n)^2)
    sigmahat <- E%*%L%*%t(E)
    
    #3.1.(3) Estimate the location paramater
    sigmahalf <- E%*%sqrt(L)%*%t(E)
    # sigmahalf <- E%*%sqrt(L)
    muhat <- sigmahalf%*%apply(x%*%solve(sigmahalf),2, median)
    
    #Select initial subset based on Mahanalobis distance
    h0 <- x[rank(Mdistance(x, muhat, sigmahat))<=ceiling(0.5*nrow(x)),]
    sigmahat <- cov(h0)
    muhat <- colMeans(h0)
    
    #Derive the new parameters based on a subset of the non-standardized data
    subset <- X_og[rank(Mdistance(x, muhat, sigmahat))<=h,]
    sigmahat <- cov(subset)
    muhat <- colMeans(subset)
    
    #Initialize sigmahat_old to always satisfy convergence criteria
    sigmahat_old <- NA
    indices <- rep(NA,h)
    
    #C-step until convergence, convergence defined as smaller than 1% decrease in the determinant of sigmahat
    while(is.na(sigmahat_old) || abs((det(sigmahat)-det(sigmahat_old))/det(sigmahat_old))>0.000001){
      # if (length(sigmahat_old)>1) {
      #   print(paste("Percentage change:", (det(sigmahat)-det(sigmahat_old))/det(sigmahat_old)))
      # }
      
      # Create subset of observations with smallest distance
      indices <- which(rank(Mdistance(X_og, muhat, sigmahat)) <= h)
      subset <- X_og[indices,]
      
      # Update estimates
      sigmahat_old <- sigmahat
      sigmahat <- cov(subset)
      muhat <- colMeans(subset)
    }
    
    output <- list()
    output[["sigmahat"]] <- sigmahat
    output[["muhat"]] <- muhat
    output[["indices"]] <- indices
    
    MDCRawtemp[[paste(c("S",S),collapse = "")]] <- output
  }
  
  #Now we find the smallest determinant of the covariance of all 6 methods
  determinants <- rep(NA,6)
  for(i in 1:6){
    determinants[i] <- det(MDCRawtemp[[i]][["sigmahat"]])
  }
  
  #paste the center, cov and indices in the final output list
  output[["raw.center"]] <- MDCRawtemp[[which(determinants == min(determinants))[1]]][["muhat"]]
  output[["raw.cov"]] <- MDCRawtemp[[which(determinants == min(determinants))[1]]][["sigmahat"]]
  output[["best"]] <- MDCRawtemp[[which(determinants == min(determinants))[1]]][["indices"]]
  # cat("Fraction of data used in unweighted: ",length(output[["best"]])/nrow(X_og),'\n')
  
  #Reweighting step
  chisquare <- qchisq(0.975, df = ncol(X_og))
  subset2 <- X_og[Mdistance(X_og, output[["raw.center"]], output[["raw.cov"]])^2<=chisquare,]
  
  #Adding the reweighted parameters & weights to the output list
  output[["center"]] <- colMeans(subset2)
  output[["cov"]] <- cov(subset2)
  output[["weights"]] <- as.integer(as.logical(Mdistance(X_og, output[["raw.center"]], output[["raw.cov"]])^2<=chisquare))
  # cat("Fraction of data used in reweighted: ", sum(output[["weights"]])/nrow(X_og))
  return(output)
  
  # Please note that the subset sizes for the MCD are not simply fractions of 
  # the number of observations in the data set, as discussed in the lectures.
  # You can use function h.alpha.n() from package robustbase to compute the 
  # subset size.
}

## Function for regression based on the deterministic MCD

# Input:
# x ........ matrix of explanatory variables
# y ........ response variable
# alpha .... proportion of observations to be used for the subset size in the 
#            MCD estimator
# anything else you need

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMCD <- function(x, y, alpha, ...) {
  #Combine the two data sets
  X <- cbind(y, x)
  
  #Prepare output list
  output <- list()
  
  #Perform MCD on the full dataset
  output[["MCD"]] <- covDetMCD(X, alpha = alpha)
  
  #Extract the cov and the mu
  cov <- output[["MCD"]]$cov
  mu <- output[["MCD"]]$center
  
  covxx <- cov[-1,-1]
  covxy <- cov[-1,1]
  muy <- mu[1]
  mux <- mu[-1]
  
  #calculate coefficients
  coefficients <- list()
  coefficients[["slope"]] <- solve(covxx)%*%covxy
  coefficients[["intercept"]] <- muy - t(mux)%*%coefficients[["slope"]]
  
  #append coefficients to output
  output[["coefficients"]] <- coefficients
  
  #calculate fitted values
  output[["fitted.values"]] <- x%*%coefficients[["slope"]]+rep(coefficients[["intercept"]],nrow(X))
  
  #calculate residuals
  output[["residuals"]] <- y-output[["fitted.values"]]
  
  return(output)
}