library("mvtnorm")
data <- read.csv("Disease.csv")

Covs <- c(1:6)
Nobs <- nrow(data)
Npar <- ncol(data) -1
lambda <- 1# scaling factor for the prior of beta 
# Setting up the prior
mu <- as.matrix(rep(0,Npar)) # Prior mean vector
Sigma <- (1/lambda)*diag(Npar) # Prior covariance matrix

y <- data$class_of_diagnosis
X <- data[,1:6]
Xnames <- colnames(X)
X <- as.matrix(X)
# Functions that returns the log posterior for the logistic and probit regression.
# First input argument of this function must be the parameters we optimize on, 
# i.e. the regression coefficients beta.

LogPostLogistic <- function(betas,y,X,mu,Sigma){
  print(dim(betas))
  print(dim(X))
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

initVal <- matrix(0,Npar,1)

OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

approxPostMode <- matrix(OptimRes$par,1, length(Covs))
colnames(approxPostMode) <- Xnames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(solve(-OptimRes$hessian))) # Computing approximate standard deviations.
names(approxPostStd) <- Xnames # Naming the coefficient by covariates
Cred_int <- matrix(0,2,length(Covs)) # Create 95 % approximate credibility intervals for each coefficient
Cred_int[1,] <- approxPostMode - 1.96*approxPostStd
Cred_int[2,] <- approxPostMode + 1.96*approxPostStd
colnames(Cred_int) <- Xnames
rownames(Cred_int) <- c("LCI","UCI") # Lower and upper limits of the approximate credibility intervals
print('The posterior mode is:')
print(approxPostMode)
print('The approximate posterior standard deviation is:')
print(approxPostStd)
print('95 % approximate credibility intervals for each coefficient:')
print(Cred_int)#
print("Negative Hessian Matrix evaluated at Posterior Median")
print(-OptimRes$hessian)
print("J(θ)")
print(solve(-OptimRes$hessian))