# Preparation
library(MASS)
set.seed(17)
source("DeterministicMCD.R")
norm_vec <- function(x) sqrt(sum(x^2))
MSE<- function(y_hat,y) mean((y_hat-y)^2)

# Setting sample size and true parameter values 
n=1000
trueParams=c(1,2)


simulation <- function(epsilon=0.05,contamination=0,R=100,k=NaN) {
    # Runs a simulation with the following paramaters:
    # @epsilon: contamination level
    # @contamination: contamination type (0,1,2 or 3 as defined in the report)
    # @R: number of replications
    
  # Define matrices for storing results
  biases=matrix(nrow = R, ncol = 3)
  MSEs=matrix(nrow = R, ncol = 3)
  
  # Simulation loop
  for (r in 1:R) {
    # Data generation
    x=rnorm(n,0,1)
    y=trueParams[1]+trueParams[2]*x+rnorm(n,0,1)
    
    # Generate test set for prediction
    x_test = rnorm(20,0,1)
    y_test=trueParams[1]+trueParams[2]*x_test+rnorm(20,0,1)
    
    # Data contamination
    if (is.nan(k)) k=runif(1,-6,-3)+rbinom(1,1,0.5)*runif(1,9,12) # contamination parameter
    contaminated=sample(1:n,round(n*epsilon)) #subsample to be contaminated

    if (contamination==1) {
    # 1) Good levarage points
      x[contaminated]=x[contaminated]+k
      y[contaminated]=y[contaminated]+trueParams[2]*k
      
    } else if (contamination==2) {
    # 2) Vertical outliers
      #x[contaminated]=rnorm(round(n*epsilon),0,0.5)
      y[contaminated]=y[contaminated]+k
    }  else if (contamination==3) {
    # 3) Bad leverage points
      x[contaminated]=x[contaminated]-k
      y[contaminated]=-y[contaminated]+trueParams[2]*k
      
    }
    
    # Estimation
    alpha=0.7
    ols <- lm(y ~ x)
    lts <- ltsReg(x, y, alpha)
    #mcd <- lmDetMCD(x, y, alpha)
    
    # Evaluation
    
    # Bias
    biases[r,1]=norm_vec(abs(trueParams-ols$coefficients))
    biases[r,2]=norm_vec(abs(trueParams-lts$coefficients))
    #bias[r,3]=norm_vec(abs(trueParams-mcd$coefficients))
    
    # MSE
    y_hat_ols=ols$coefficients[1]+ols$coefficients[2]*x_test
    y_hat_lts=lts$coefficients[1]+lts$coefficients[2]*x_test
    #y_hat_mcd=mcd$coefficients[1]+mcd$coefficients[2]*x_test
    MSEs[r,1] = MSE(y_hat_ols,y_test)
    MSEs[r,2] = MSE(y_hat_lts,y_test)
    #MSE[r,3] = MSE(y_hat_mcd,y_test)
  } 
  return(list("bias"=colMeans(biases),
              "MSE"=colMeans(MSEs),
              "y"=y,"x"=x,"contaminated"=contaminated,"k"=k))
}


# Plotting the different contamination processes
par(mfrow=c(2,2),
    mai = c(0.7, 0.7, 0.2, 0.2))
for (c in 0:3) {
  sim=simulation(R=1,contamination=c,k=4)
  if (c==1) {
    title="Good levarage points"
  } else if (c==2) {
    title="Vertical outliers"
  } else if (c==3) {
    title="Bad leverage points"
  } else {
    title="No contamination"
  }
  plot(sim$x[-sim$contaminated],sim$y[-sim$contaminated],pch = 16, cex = 0.4, col = "blue",
       xlab="x",ylab="y",ylim=c(min(sim$y)-1,max(sim$y)+1), xlim=c(min(sim$x)-1,max(sim$x)+1),
       main=title)
  points(sim$x[sim$contaminated],sim$y[sim$contaminated],pch = 16, cex = 0.4, col = "red")
  abline(lm(sim$y[-sim$contaminated] ~ sim$x[-sim$contaminated]))
}
par(mfrow=c(1,1))



