## Preparation
library(MASS)
library(ggplot2)
library(latex2exp)
set.seed(17)
source("DeterministicMCD.R")
norm_vec <- function(x) sqrt(sum(x^2))
MSE<- function(y_hat,y) mean((y_hat-y)^2)

# Setting sample size and true parameter values 
n=1000
trueParams=c(1,2)

#' Runs a simulation with the following paramaters:
#' @param epsilon: contamination level
#' @param contamination: contamination type (0,1,2 or 3 as defined in the report)
#' @param R: number of replications
simulation <- function(epsilon=0.05,contamination=0,R=100,k=NaN) {

  # Define matrices for storing results
  intercept_biases=matrix(nrow = R, ncol = 5)
  slope_biases=matrix(nrow=R,ncol=5)
  MSEs=matrix(nrow = R, ncol = 5)
  
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
    ols <- lm(y ~ x)
    lts1 <- ltsReg(x, y, alpha=0.5)
    mcd1 <- lmDetMCD(x, y, alpha=0.5)
    lts2 <- ltsReg(x, y, alpha=0.75)
    mcd2 <- lmDetMCD(x, y, alpha=0.75)
    
    # Evaluation
    
    # Bias
    intercept_biases[r,1]=abs(trueParams[1]-ols$coefficients[1])
    intercept_biases[r,2]=abs(trueParams[1]-lts1$coefficients[1])
    intercept_biases[r,3]=abs(trueParams[1]-mcd1$coefficients[["intercept"]])
    intercept_biases[r,4]=abs(trueParams[1]-lts2$coefficients[1])
    intercept_biases[r,5]=abs(trueParams[1]-mcd2$coefficients[["intercept"]])
    slope_biases[r,1]=abs(trueParams[2]-ols$coefficients[2])
    slope_biases[r,2]=abs(trueParams[2]-lts1$coefficients[2])
    slope_biases[r,3]=abs(trueParams[2]-mcd1$coefficients[["slope"]])
    slope_biases[r,4]=abs(trueParams[2]-lts2$coefficients[2])
    slope_biases[r,5]=abs(trueParams[2]-mcd2$coefficients[["slope"]])
    
    # MSE
    y_hat_ols=ols$coefficients[1]+ols$coefficients[2]*x_test
    y_hat_lts1=lts1$coefficients[1]+lts1$coefficients[2]*x_test
    y_hat_mcd1=as.vector(mcd1$coefficients[["intercept"]])+as.vector(mcd1$coefficients[["slope"]])*x_test
    y_hat_lts2=lts2$coefficients[1]+lts2$coefficients[2]*x_test
    y_hat_mcd2=as.vector(mcd2$coefficients[["intercept"]])+as.vector(mcd2$coefficients[["slope"]])*x_test
    MSEs[r,1] = MSE(y_hat_ols,y_test)
    MSEs[r,2] = MSE(y_hat_lts1,y_test)
    MSEs[r,3] = MSE(y_hat_mcd1,y_test)
    MSEs[r,4] = MSE(y_hat_lts2,y_test)
    MSEs[r,5] = MSE(y_hat_mcd2,y_test)
  } 
  return(list("intercept_bias"=colMeans(intercept_biases),
              "slope_bias"=colMeans(slope_biases),
              "MSE"=colMeans(MSEs),
              "y"=y,"x"=x,"contaminated"=contaminated,"k"=k))
}


## Plotting the different contamination processes
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
       main=title,cex.main=0.7)
  points(sim$x[sim$contaminated],sim$y[sim$contaminated],pch = 16, cex = 0.4, col = "red")
  abline(lm(sim$y[-sim$contaminated] ~ sim$x[-sim$contaminated]))
}
par(mfrow=c(1,1))



## Running simulations

# 1) contamination level: 10%
results1=matrix(nrow = 15, ncol = 4)
rownames(results1)=c("OLS intercept bias","OLS slope bias","OLS MSE",
                    "LTS(alpha=0.5) intercept bias","LTS(alpha=0.5) slope bias","LTS(alpha=0.5) MSE",
                    "MCD(alpha=0.5) intercept bias","MCD(alpha=0.5) slope bias","MCD(alpha=0.5) MSE",
                    "LTS(alpha=0.75) intercept bias","LTS(alpha=0.75) slope bias","LTS(alpha=0.75) MSE",
                    "MCD(alpha=0.75) intercept bias","MCD(alpha=0.75) slope bias","MCD(alpha=0.75) MSE")
colnames(results1)=c(0,1,2,3)
for (c in 0:3) {
  sim=simulation(R=100,epsilon=0.10,contamination=c)
  results1[1,(c+1)]=round(sim$intercept_bias[1],3)
  results1[2,(c+1)]=round(sim$slope_bias[1],3)
  results1[3,(c+1)]=round(sim$MSE[1],3)
  results1[4,(c+1)]=round(sim$intercept_bias[2],3)
  results1[5,(c+1)]=round(sim$slope_bias[2],3)
  results1[6,(c+1)]=round(sim$MSE[2],3)
  results1[7,(c+1)]=round(sim$intercept_bias[3],3)
  results1[8,(c+1)]=round(sim$slope_bias[3],3)
  results1[9,(c+1)]=round(sim$MSE[3],3)
  results1[10,(c+1)]=round(sim$intercept_bias[4],3)
  results1[11,(c+1)]=round(sim$slope_bias[4],3)
  results1[12,(c+1)]=round(sim$MSE[4],3)
  results1[13,(c+1)]=round(sim$intercept_bias[5],3)
  results1[14,(c+1)]=round(sim$slope_bias[5],3)
  results1[15,(c+1)]=round(sim$MSE[5],3)
}

#2) contamination level: 30%

results2=matrix(nrow = 15, ncol = 4)
rownames(results2)=c("OLS intercept bias","OLS slope bias","OLS MSE",
                     "LTS(alpha=0.5) intercept bias","LTS(alpha=0.5) slope bias","LTS(alpha=0.7) MSE",
                     "MCD(alpha=0.5) intercept bias","MCD(alpha=0.5) slope bias","MCD(alpha=0.7) MSE",
                     "LTS(alpha=0.75) intrecept bias","LTS(alpha=0.75) slope bias","LTS(alpha=0.85) MSE",
                     "MCD(alpha=0.75) intrecept bias","MCD(alpha=0.75) slope bias","MCD(alpha=0.85) MSE")
colnames(results2)=c(0,1,2,3)
for (c in 0:3) {
  sim=simulation(R=100,epsilon=0.30,contamination=c)
  results2[1,(c+1)]=round(sim$intercept_bias[1],3)
  results2[2,(c+1)]=round(sim$slope_bias[1],3)
  results2[3,(c+1)]=round(sim$MSE[1],3)
  results2[4,(c+1)]=round(sim$intercept_bias[2],3)
  results2[5,(c+1)]=round(sim$slope_bias[2],3)
  results2[6,(c+1)]=round(sim$MSE[2],3)
  results2[7,(c+1)]=round(sim$intercept_bias[3],3)
  results2[8,(c+1)]=round(sim$slope_bias[3],3)
  results2[9,(c+1)]=round(sim$MSE[3],3)
  results2[10,(c+1)]=round(sim$intercept_bias[4],3)
  results2[11,(c+1)]=round(sim$slope_bias[4],3)
  results2[12,(c+1)]=round(sim$MSE[4],3)
  results2[13,(c+1)]=round(sim$intercept_bias[5],3)
  results2[14,(c+1)]=round(sim$slope_bias[5],3)
  results2[15,(c+1)]=round(sim$MSE[5],3)
}




# Look at results
as.table(results1)
as.table(results2)
# Save results
xresults=cbind.data.frame(results1,results2[,2:4])
xtab<-xtable(as.table(as.matrix(xresults)))
write.csv(xresults,"simulationResults.csv")


# Plotting MCD results
par(mfrow=c(1,3),
    mai = c(0.7, 0.7, 0.2, 0.2))
for (c in 1:3) {
  if (c==1) {
    title="Type 1"
  } else if (c==2) {
    title="Type 2"
  } else if (c==3) {
    title="Type 3"
  }
  sim=simulation(R=1,epsilon=0.10,contamination=c,k=6)
  data=cbind(sim$x,sim$y)
  mcd <- covDetMCD(data, 0.75)
  data2 <- as.data.frame(cbind(data, mcd[["weights"]]))
  trendline=lm(data2$V2[as.logical(data2$V3)] ~ data2$V1[as.logical(data2$V3)])
  data2$V3[data2$V3==1]="blue"
  data2$V3[data2$V3==0]="black"
  plot(data2$V1,data2$V2,col=data2$V3,pch=20,xlab="x",ylab="y",
       main=title)
  abline(trendline,col="blue")
  abline(lm(sim$y ~ sim$x))
  legend("topright", legend=c("OLS",expression(paste("Plug-in, ",alpha, "=75%"))),
         col=c("black", "blue"),cex=0.8, text.font=1,lty=1,lwd=2)
}
par(mfrow=c(1,1))

length(data2$V3[data2$V3=="blue"])/length(data2$V3)


