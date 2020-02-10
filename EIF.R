
library(robustbase)
library(tictoc)

# Function for plug in robust estimator is loaded from
source("DeterministicMCD.R")

#Eredivisie28 <- load(file.choose())

# Checking outliers
n=length(Eredivisie28$Age)

LS <- lm(Eredivisie28$MarketValue ~ Eredivisie28$Age)
studres(LS)
plot(studres(LS),ylim=c(-2,5.2),ylab="Studentized Residuals")
lines(studres(LS),ylim=c(-2,5.2))
abline(h=qt(0.975,df=n-3),col="blue")
abline(h=-qt(0.975,df=n-3),col="blue")

# Plotting original data and regression lines
plugin <- lmDetMCD(Eredivisie28$Age, Eredivisie28$MarketValue, alpha=0.8)
plot(Eredivisie28$Age, Eredivisie28$MarketValue, xlab="Age", ylab="Market Value", type="p")
points(Eredivisie28$Age[abs(studres(LS))>qt(0.975,df=n-3)],
       Eredivisie28$MarketValue[abs(studres(LS))>qt(0.975,df=n-3)],pch=16,col="orange")
abline(lm(MarketValue ~ Age, data=Eredivisie28), col="red")
abline(ltsReg(MarketValue ~ Age, data=Eredivisie28, alpha=0.8), col="blue")
abline(a = plugin$coefficients[["intercept"]], b=plugin$coefficients[["slope"]], col="darkgreen")
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=FALSE, bty="n", ncol=3)

# # Function for estimating OLS coefficients
# ols <- function(x, y, add_intercept=TRUE) {
#   x <- as.matrix(x)
#   y <- as.matrix(y)
#   if (add_intercept) {
#     x <- cbind(rep(1,nrow(x)),x)
#   }
#   beta <- solve(t(x)%*%x)%*%(t(x)%*%y)
#   return(beta)
# }

# Function for creating replacements
repl <- function(y, length=200) {
  range <- max(y)-min(y)
  seq( (median(y)-range), (median(y)+range), length.out = length)
}

# Function for estimating EIF by changing the y-value of one observation
# y: dependent variable
# x: predictor variable
# dependent: boolean indicating whether we change x or y.
# alpha: percentage of observations to use for ltsReg (do we need this though?)
eif <- function(x, y, dependent=TRUE, alpha) {
  
  # Determining the test range, that is the range of values to plug in
  if (dependent) {
    replacements= repl(y)
  } else  {
    replacements= repl(x)
  }
  
  # Randomly selecting the value to change
  n=length(x)
  i=sample(1:n,1)
  
  # Baseline estimates (i.e. no change in observed values)
  baseline_ols <- lm(y ~ x)$coefficients
  baseline_lts <- ltsReg(x, y, alpha = alpha)$coefficients
  baseline_mcd <- rep(NA,2)
  lmMCD <- lmDetMCD(x, y, alpha = alpha)
  baseline_mcd[1] <- lmMCD$coefficients[["intercept"]]
  baseline_mcd[2] <- lmMCD$coefficients[["slope"]]
  
  # Define matrices for storing coefficients
  beta_ols <- matrix(nrow = length(replacements), ncol = 2)
  beta_lts <- matrix(nrow = length(replacements), ncol = 2)
  beta_mcd <- matrix(nrow = length(replacements), ncol = 2)
  
  # Looping through the range of replacements and plugging in different values
  # to the ith datapoint and then running the three estimators
  for (j in 1:length(replacements)) {
    
    # Replace y-value if dependent is True
    if (dependent) {
      adj_y=y
      adj_y[i]=replacements[j]
      adj_x=x
    } else {
    # Replace x-value if dependent is False
      adj_x=x  
      adj_x[i]=replacements[j]
      adj_y=y
    } 
    
    # Save estimates in matrices
    beta_ols[j,] <- lm(adj_y ~ adj_x)$coefficients
    beta_lts[j,] <- ltsReg(adj_x, adj_y, alpha = alpha)$coefficients
    lmMCD <- lmDetMCD(adj_x, adj_y, alpha = alpha)
    beta_mcd[j,1] <- lmMCD$coefficients[["intercept"]]
    beta_mcd[j,2] <- lmMCD$coefficients[["slope"]]
  }
  
  # Compute EIFs
  eif_ols <- n*apply(beta_ols, 1, function(x) x-baseline_ols)
  eif_lts <- n*apply(beta_lts, 1, function(x) x-baseline_lts)
  eif_mcd <- n*apply(beta_mcd, 1, function(x) x-baseline_mcd)
  
  return(list(beta_ols=beta_ols, beta_lts=beta_lts,beta_mcd=beta_mcd,
              eif_ols=eif_ols, eif_lts=eif_lts,eif_mcd=eif_mcd,
              observation=i, replacements = replacements))
  
}


set.seed(0)
# Main code (for changing y)
tic()
result <- eif(Eredivisie28$Age, Eredivisie28$MarketValue, TRUE, alpha=0.75)
toc()

replace <- result$replacements

# Plotting for y
options(scipen=999)
par(mfrow=c(1,2))
plot(replace/1000000, result$eif_ols[1,],col="red",ylab="Change in intercept (times n)",xlab="Market Value (in millions)",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$MarketValue)),x=Eredivisie28$MarketValue/1000000,col="gray")
points(y=0,x=Eredivisie28$MarketValue[result$observation]/1000000,pch=16)
lines(replace/1000000, result$eif_lts[1,],col="blue",type = "l",lwd=2)
lines(replace/1000000, result$eif_mcd[1,],col="darkgreen",type = "l",lwd=2)
legend("topleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2)
       #,inset=c(0,1), xpd=TRUE, horiz=FALSE, bty="n", ncol=3)

plot(replace/1000000, result$eif_ols[2,],col="red",ylab="Change in slope (times n)",xlab="Market Value (in millions)",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$MarketValue)),x=Eredivisie28$MarketValue/1000000,col="gray")
points(y=0,x=Eredivisie28$MarketValue[result$observation]/1000000,pch=16)
lines(replace/1000000, result$eif_lts[2,],col="blue",type = "l",lwd=2)
lines(replace/1000000, result$eif_mcd[2,],col="darkgreen",type = "l",lwd=2)
legend("topright", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2)
       #,inset=c(0,1), xpd=TRUE, horiz=FALSE, bty="n", ncol=3)
par(mfrow=c(1,1))


# Main code (for changing x)
set.seed(0)
tic()
result <- eif(Eredivisie28$Age, Eredivisie28$MarketValue, FALSE, alpha=0.75)
toc()

replace <- result$replacements

# Plotting for x
par(mfrow=c(1,2))
plot(replace, result$eif_ols[1,],col="red",ylab="Change in intercept (times n)",xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[result$observation],pch=16)
lines(replace, result$eif_lts[1,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[1,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", ncol=3)

plot(replace, result$eif_ols[2,],col="red",ylab="Change in slope (times n)",xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[result$observation],pch=16)
lines(replace, result$eif_lts[2,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[2,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", ncol3)
par(mfrow=c(1,1))

# ######################################

set.seed(0)
# Main code (for changing y)
tic()
result <- eif(Eredivisie28$Age, log(Eredivisie28$MarketValue), TRUE, alpha=0.75)
toc()

replace <- result$replacements

# Plotting for y
par(mfrow=c(1,2))
plot(replace, result$eif_ols[1,],col="red",ylab="Change in intercept (times n)",xlab="(Logarithm of) Market Value",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$MarketValue)),x=log(Eredivisie28$MarketValue),col="gray")
points(y=0,x=log(Eredivisie28$MarketValue[result$observation]),pch=16)
lines(replace, result$eif_lts[1,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[1,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=FALSE, bty="n", ncol=3)

plot(replace, result$eif_ols[2,],col="red",ylab="Change in slope (times n)",xlab="(Logarithm of) Market Value",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$MarketValue)),x=log(Eredivisie28$MarketValue),col="gray")
points(y=0,x=log(Eredivisie28$MarketValue[result$observation]),pch=16)
lines(replace, result$eif_lts[2,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[2,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=FALSE, bty="n", ncol=3)
par(mfrow=c(1,1))


# Main code (for changing x)
set.seed(0)
tic()
result <- eif(Eredivisie28$Age, log(Eredivisie28$MarketValue), FALSE, alpha=0.75)
toc()

replace <- result$replacements

# Plotting for x
par(mfrow=c(1,2))
plot(replace, result$eif_ols[1,],col="red",ylab="Change in intercept (times n)",xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[result$observation],pch=16)
lines(replace, result$eif_lts[1,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[1,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", ncol=3)

plot(replace, result$eif_ols[2,],col="red",ylab="Change in slope (times n)",xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[result$observation],pch=16)
lines(replace, result$eif_lts[2,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[2,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", ncol3)
par(mfrow=c(1,1))

