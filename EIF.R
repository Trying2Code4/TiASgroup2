# Eredivisie28 <- load("~/Master QM/Block 3/Topics in Advanced Statistics/Group Assignment/Eredivisie28.RData")
Eredivisie28 <- load(file.choose())
plot(Eredivisie28$Age, Eredivisie28$MarketValue, type="p")

library(robustbase)
library(tictoc)

# Function for estimating OLS coefficients
ols <- function(x, y, add_intercept=TRUE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (add_intercept) {
    x <- cbind(rep(1,nrow(x)),x)
  }
  beta <- solve(t(x)%*%x)%*%(t(x)%*%y)
}

# Function for estimating EIF by changing the y-value of one observation
eif <- function(x, y, shift, obs, alpha) {
  
  # Baseline estimates (i.e. no change in observed values)
  baseline_ols <- ols(x, y)
  baseline_lts <- ltsReg(x, y, alpha)$coefficients
  # baseline_mcd <- lmDetMCD(x, y, alpha)$coefficients
  
  # Define matrices for storing coefficients
  beta_ols <- matrix(nrow = length(shift), ncol = 2)
  beta_lts <- matrix(nrow = length(shift), ncol = 2)
  # beta_mcd <- matrix(nrow = length(shift), ncol = 2)
  
  for (i in 1:length(shift)) {
    
    # Shift y-value of selected observation
    adj_y <- y
    adj_y[obs] <- adj_y[obs] + shift[i]
    
    # Save estimates in matrices
    beta_ols[i,] <- ols(as.matrix(x), as.matrix(adj_y))
    beta_lts[i,] <- ltsReg(x, adj_y, alpha)$coefficients
    # beta_mcd[i,] <- lmDetMCD(x, adj_y, alpha)$coefficients
  }
  
  # Compute EIFs
  eif_ols <- length(x)*apply(beta_ols, 1, function(x) x-baseline_ols)
  eif_lts <- length(x)*apply(beta_lts, 1, function(x) x-baseline_lts)
  # eif_mcd <- nrow(Eredivisie28)*apply(beta_mcd, 1, function(x) x-baseline_mcd)

  return(list(beta_ols=beta_ols, beta_lts=beta_lts, eif_ols=eif_ols, eif_lts=eif_lts))
}

# Main code
alpha = 0.75
shift <- seq(-5000000, 5000000, length.out = 200)
random_int <- sample(1:nrow(Eredivisie28),1)
tic()
result <- eif(Eredivisie28$Age, Eredivisie28$MarketValue, shift, obs=random_int, alpha=alpha)
toc()
# Plotting
par(mfrow=c(1,2))
plot(shift/1000000, result$eif_ols[1,],col="red",ylab="Intercept",xlab="Change in observation (x10^6)",type = "l",lty=1,lwd=2)
lines(shift/1000000, result$eif_lts[1,],col="blue",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS"),
       col=c("red", "blue"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")

plot(shift/1000000, result$eif_ols[2,],col="red",ylab="Slope",xlab="Change in observation (x10^6)",type = "l",lty=1,lwd=2)
lines(shift/1000000, result$eif_lts[2,],col="blue",type = "l",lwd=2)
legend("bottomleft", legend=c("OLS", "LTS"),
       col=c("red", "blue"),cex=0.8, text.font=1,lty=1,lwd=2,
       inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
par(mfrow=c(1,1))
