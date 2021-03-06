library(robustbase)

# Function for plug in robust estimator is loaded from
source("DeterministicMCD.R")

par(mfrow=c(1,2))
# Checking outliers
n=length(Eredivisie28$Age)
LS <- lm(log(Eredivisie28$MarketValue) ~ Eredivisie28$Age)
studres(LS)
plot(studres(LS),ylim=c(-3,3),ylab="Studentized Residuals")
lines(studres(LS),ylim=c(-2,3))
abline(h=qt(0.975,df=n-3),col="blue")
abline(h=-qt(0.975,df=n-3),col="blue")

# Plotting original data and regression lines
plugin <- lmDetMCD(Eredivisie28$Age, log(Eredivisie28$MarketValue), alpha=0.75)
plot(Eredivisie28$Age, log(Eredivisie28$MarketValue), xlab="Age", ylab="(Log) Market Value", type="p")
points(Eredivisie28$Age[abs(studres(LS))>qt(0.975,df=n-3)],
       log(Eredivisie28$MarketValue)[abs(studres(LS))>qt(0.975,df=n-3)],pch=16,col="orange")
abline(lm(log(Eredivisie28$MarketValue) ~ Eredivisie28$Age), col="red")
abline(ltsReg(log(Eredivisie28$MarketValue) ~ Eredivisie28$Age, alpha=0.75), col="blue")
abline(a = plugin$coefficients[["intercept"]], b=plugin$coefficients[["slope"]], col="darkgreen")
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.5, text.font=1,lty=1,lwd=2,pt.cex=1,
       bty="n",
       inset=c(0,1),xpd=TRUE,horiz=TRUE)


# Function for estimating EIF by changing the y-value of one observation
# y: dependent variable
# x: predictor variable
# dependent: boolean indicating whether we change x or y.
# alpha: percentage of observations to use for ltsReg (do we need this though?)
eif <- function(x, y, dependent=TRUE, alpha) {
  
  # Function for creating replacements
  repl <- function(y, length=200,dependent) {
    if (dependent) {
      return(seq(0,30,length.out = length))
    } else {
      return(seq(0,60,length.out = length))  
    }
  }
  
  # Determining the test range, that is the range of values to plug in
  if (dependent) {
    replacements= repl(y,length=200,TRUE)
  } else  {
    replacements= repl(x,length=200,FALSE)
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

#___________________________________________________
## Creating EIF plots

set.seed(0)
# Main code (for changing y)
result <- eif(Eredivisie28$Age, Eredivisie28$MarketValue, TRUE, alpha=0.75)

replace <- result$replacements

# Plotting for y
options(scipen=999)
par(mfrow=c(1,2))
plot(replace/1000000, result$eif_ols[1,],col="red",ylab="Change in intercept (times n)",xlab="Market Value (in millions)",
     type = "l",lty=1,lwd=2)
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
result <- eif(Eredivisie28$Age, Eredivisie28$MarketValue, FALSE, alpha=0.5)

replace <- result$replacements

# Plotting for x
par(mfrow=c(1,2))
plot(replace, result$eif_ols[1,],col="red",ylab="Change in intercept (times n)",
     ylim=c(-20000000,50000000),xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[result$observation],pch=16)
lines(replace, result$eif_lts[1,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[1,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", ncol=3)

plot(replace, result$eif_ols[2,],col="red",ylab="Change in slope (times n)"
     ,ylim=c(-3000000,3000000),xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[result$observation],pch=16)
lines(replace, result$eif_lts[2,],col="blue",type = "l",lwd=2)
lines(replace, result$eif_mcd[2,],col="darkgreen",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", ncol3)
par(mfrow=c(1,1))

#___________________________________________________
# Using Logged Marketvalue

set.seed(0)
# Main code (for changing y)
eifylog <- eif(Eredivisie28$Age, log(Eredivisie28$MarketValue), TRUE, alpha=0.50)
replace_y <- eifylog$replacements

# Plotting for y
par(mfrow=c(1,2))
plot(replace_y, eifylog$eif_ols[1,],col="red",ylab="Change in intercept (times n)",
     xlab="(Log of) Market Value",type = "l",lty=1,lwd=2,ylim=c(-100,100))
points(y=rep(0,length(Eredivisie28$MarketValue)),x=log(Eredivisie28$MarketValue),col="gray")
points(y=0,x=log(Eredivisie28$MarketValue[eifylog$observation]),pch=16)
lines(replace_y, eifylog$eif_lts[1,],col="blue",type = "l",lwd=2)
lines(replace_y, eifylog$eif_mcd[1,],col="darkgreen",type = "l",lwd=2)
legend("topleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2)
plot(replace_y, eifylog$eif_ols[2,],col="red",ylab="Change in slope (times n)",
     xlab="(Log of) Market Value",type = "l",lty=1,lwd=2,ylim=c(-1,2))
points(y=rep(0,length(Eredivisie28$MarketValue)),x=log(Eredivisie28$MarketValue),col="gray")
points(y=0,x=log(Eredivisie28$MarketValue[eifylog$observation]),pch=16)
lines(replace_y, eifylog$eif_lts[2,],col="blue",type = "l",lwd=2)
lines(replace_y, eifylog$eif_mcd[2,],col="darkgreen",type = "l",lwd=2)
legend("topright", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2)
par(mfrow=c(1,1))


# Main code (for changing x)
set.seed(0)
eifxlog <- eif(Eredivisie28$Age, log(Eredivisie28$MarketValue), FALSE, alpha=0.5)
replace_x <- eifxlog$replacements

# Plotting for x
par(mfrow=c(1,2))
plot(replace_x, eifxlog$eif_ols[1,],col="red",ylab="Change in intercept (times n)",
     ylim=c(-250,240),xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[eifxlog$observation],pch=16)
lines(replace_x, eifxlog$eif_lts[1,],col="blue",type = "l",lwd=2)
lines(replace_x, eifxlog$eif_mcd[1,],col="darkgreen",type = "l",lwd=2)
legend("topleft", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2)
plot(replace_x, eifxlog$eif_ols[2,],col="red",ylab="Change in slope (times n)",
     ylim=c(-20,20),xlab="Age",type = "l",lty=1,lwd=2)
points(y=rep(0,length(Eredivisie28$Age)),x=Eredivisie28$Age,col="gray")
points(y=0,x=Eredivisie28$Age[eifxlog$observation],pch=16)
lines(replace_x, eifxlog$eif_lts[2,],col="blue",type = "l",lwd=2)
lines(replace_x, eifxlog$eif_mcd[2,],col="darkgreen",type = "l",lwd=2)
legend("topright", legend=c("OLS", "LTS", "Plug-in"),
       col=c("red", "blue", "darkgreen"),cex=0.8, text.font=1,lty=1,lwd=2)
par(mfrow=c(1,1))

