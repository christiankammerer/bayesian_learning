library(asbio)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(matlib)


# 1

## 1.1 

data <- read.csv("temp_linkoping.csv")

beta_prior_f <- function(n, mu0, sigma_sq, omega_0_inv){
  return(rmvnorm(n, mean = mu0, sigma = sqrt(sigma_sq * omega_0_inv)))
}

sigma_prior_f <- function(n, v0, sigma_sq_0){
  rinvchisq(n, df = v0, scale = sigma_sq_0)
}

joint_prior <- function(n=1, mu0 = c(20, -100, 100), 
                        omega_0_inv = 0.01 * diag(3), 
                        v0=1, sigma_sq_0=1){
  sigma_sq <- sigma_prior_f(1, v0, sigma_sq_0)
  beta <- beta_prior_f(1, mu0, sigma_sq, omega_0_inv)
  return(beta)
}

beta_draws <- sapply(1:100, joint_prior)



# Transpose to make each row one set of coefficients
beta_draws_t <- t(beta_draws)

x_vals <- seq(1, 365, 1)/365

# Build data frame of predictions
curve_data <- lapply(1:nrow(beta_draws_t), function(i) {
  beta <- beta_draws_t[i, ]
  data.frame(
    x = x_vals,
    y = beta[1] + beta[2]*x_vals + beta[3]*x_vals^2,
    curve = paste0("Draw_", i)
  )
}) %>% bind_rows()

p <- ggplot(curve_data, aes(x = x, y = y, group = curve, color = curve)) +
  geom_line(alpha = 0.7) +
  labs(title = "Regression Curves from Prior Draws",
       x = "x", y = "y") +
  theme_minimal() +
  theme(legend.position = "none")  

print(p)


## 1.2.1

# Returns posterior parameters (mu_n, omega_n, sigma_sq_n, vn)
get_posterior_parameters <- function(mu0, omega0, sigma_sq_0, v0, y, X) {
  vn <- v0 + length(y)
  XtX <- t(X) %*% X
  omega_n <- XtX + omega0
  beta_hat <- solve(XtX) %*% t(X) %*% y
  mu_n <- solve(omega_n) %*% (XtX %*% beta_hat + omega0 %*% mu0)
  sigma_sq_n <- as.numeric(
    (v0 * sigma_sq_0 + t(y) %*% y + t(mu0) %*% omega0 %*% mu0 - t(mu_n) %*% omega_n %*% mu_n) / vn
  )
  
  return(list(
    mu_n = mu_n,
    omega_n = omega_n,
    sigma_sq_n = sigma_sq_n,
    vn = vn
  ))
}

sigma_posterior_f <- function(n = 1, posterior_params){
  with(posterior_params, {
    return(rinvchisq(n, df = vn, scale = sigma_sq_n))
  })
}

beta_posterior_f <- function(n = 1, posterior_params) {
  with(posterior_params, {
    # Sample sigma
    sigma_sq_samples <- sigma_posterior_f(n = n, posterior_params)
    
    # Sample beta conditional for each sigma
    beta_samples <- t(sapply(sigma_sq_samples, function(sigma_sq) {
      rmvnorm(1, mean = mu_n, sigma = sigma_sq * solve(omega_n))
    }))
    
    return(list(beta = beta_samples, sigma_sq = sigma_sq_samples))
  })
}

marginal_beta <- function(n, mu, Sigma, df) {
  p <- length(mu)
  chi_samples <- rchisq(n, df)
  Z <- matrix(rnorm(n * p), nrow = n)  # standard normal
  t_samples <- sweep(Z, 1, sqrt(chi_samples / df), FUN = "/") %*% chol(Sigma)
  sweep(t_samples, 2, mu, FUN = "+")
}

X <- data$time
y <- data$temp

X <- cbind(1, X, X^2)
posterior_params <- get_posterior_parameters(mu0 = c(20, -100, 100), 
                                             omega0 = 0.01 * diag(3), 
                                             sigma_sq_0 = 1, v0 = 1,
                                             y, X)

n_post = 1000
posterior_draws <- beta_posterior_f(n = n_post, posterior_params)
beta_posterior_draws <- posterior_draws$beta
sigma_posterior_draws <- posterior_draws$sigma

beta_df <- as.data.frame(beta_posterior_draws)
beta_df$beta0 <- beta_df$V1
beta_df$beta1 <- beta_df$V2
beta_df$beta2 <- beta_df$V3

beta_long <- tidyr::gather(beta_df, key = "Beta", value = "Value", beta0, beta1, beta2)

# Plot overlapping histograms
p <- ggplot(beta_long, aes(x = Value, fill = Beta)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 200) + 
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(title = "Overlapping Histograms of Betas",
       x = "Beta values", 
       y = "Frequency") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))
print(p)

p <- ggplot() +
  geom_histogram(aes(x = sigma_posterior_draws), alpha = 0.5, position = "identity", bins = 20, col = "red", fill = "red") +
  xlab("Sigma") + 
  ggtitle("Distribution of Posterior of Sigma") +
  theme_minimal()

print(p)


## 1.2.2

apply_regression <- function(beta_values, X){
  return(beta_values[1] + beta_values[2] * X[2] + beta_values[3] * X[3])
}

temp_matrix <- matrix(ncol = nrow(beta_posterior_draws), nrow = nrow(X))

for(i in 1:nrow(temp_matrix)){
  temp_matrix[i, ] <- apply(beta_posterior_draws, 1, function(betas) apply_regression(betas, X[i, ]))
}

median_temp <- apply(temp_matrix, 1, median)
fifth_perc <- apply(temp_matrix, 1, function(row) quantile(row, 0.05))
ninetyfifth_perc <- apply(temp_matrix, 1, function(row) quantile(row, 0.95))

temperature_frame <- data.frame(actual = y, median_posterior_temp = median_temp, 
                                lower_bound = fifth_perc, upper_bound = ninetyfifth_perc)

p <- ggplot(temperature_frame) +
  aes(x = 1:366) +
  geom_point(aes(y = actual, color = "Actual Temperature"), shape = 21, fill = "red") +
  geom_line(aes(y = median_posterior_temp, color = "Posterior Median"), linewidth = 1) +
  geom_line(aes(y = lower_bound), color = "blue", alpha = 0.2, linewidth = 1, show.legend = FALSE) + 
  geom_line(aes(y = upper_bound), color = "blue", alpha = 0.2, linewidth = 1, show.legend = FALSE) + 
  scale_color_manual(values = c("Actual Temperature" = "red", 
                                "Posterior Median" = "blue")) + 
  xlab("Number of Days since first Observation") + 
  ylab("Temperature in °Celsius") + 
  ggtitle("Median Regression Curve vs Actual Observations") + 
  theme_minimal() + 
  theme(legend.position = c(0.5, 0.8))

print(p)


## 1.3


# Posterior distribution of value that minimizes regression curve
posterior_x <- function(betas){
  return(-1/2 * (betas[2]/betas[3]))
}
x_s <- apply(beta_posterior_draws, 1, posterior_x) * 365


## 1.4

X <- data$time
X <- cbind(1, X, X^2, X^3, X^4, X^5, X^6, X^7, X^8, X^9, X^10)

m <- 10
n <- 1000
eta0 <- 10
lambda0 <- 1
mu0 <- c(20, -100, 100, -10, 10, -1, 1, -0.1, 0.1, -0.01, 0.01)
sigma_sq_0 <- 1 
v0 <- 1
lambda_prior_draws <- rinvchisq(n, eta0, lambda0) # prior of lambda
sigma_prior_draws <- rinvchisq(n, v0, sigma_sq_0) 
beta_prior_draws <- t(sapply((1:n), function(i) rmvnorm(n = 1, mean = mu0,
                                                        sigma = (sigma_prior_draws[i] * 1/lambda_prior_draws[i] * diag(m + 1)))))
temp_reg_values <- apply(beta_prior_draws, 1, function(betas) X %*% betas)

time_values <- data$time
temp_reg_df <- as.data.frame(temp_reg_values)
temp_reg_df$time <- time_values

# Convert to long format
long_df <- pivot_longer(temp_reg_df, cols = -time, names_to = "draw", values_to = "y")

# Assuming 'long_df' is structured like curve_data:
# long_df has columns: time (x-axis), y (y-axis), draw (equivalent to "curve")

p <- ggplot(long_df, aes(x = time, y = y, group = draw, color = draw)) +
  geom_line(alpha = 0.7) +  
  labs(title = "Posterior Regression Lines",
       x = "Time", y = "Predicted Value") +
  ylim(-40, 40) +
  theme_minimal() +
  theme(legend.position = "none")  

print(p)



## 2.1

library(mvtnorm)
data <- read.csv("Disease.csv")

Covs <- c(1:6)
Nobs <- nrow(data)
Npar <- ncol(data) -1
tau <- 2
mu <- as.matrix(rep(0,Npar)) # Prior mean vector
Sigma <- tau^2*diag(Npar) # Prior covariance matrix

y <- data$class_of_diagnosis
X <- data[,1:6]
Xnames <- colnames(X)
X <- as.matrix(X)
# Functions that returns the log posterior for the logistic and probit regression.
# First input argument of this function must be the parameters we optimize on, 
# i.e. the regression coefficients beta.

LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

initVal <- matrix(0,Npar,1)

OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

approxPostMode <- matrix(OptimRes$par,1, length(Covs)) # beta-tilde (MLE)
colnames(approxPostMode) <- Xnames 
approxPostStd <- sqrt(diag(solve(-OptimRes$hessian))) # Computing approximate standard deviations.
names(approxPostStd) <- Xnames 

Cred_int <- matrix(0,2,length(Covs)) # Create 95 % approximate credibility intervals for each coefficient
Cred_int[1,] <- approxPostMode - 1.96*approxPostStd # LCI
Cred_int[2,] <- approxPostMode + 1.96*approxPostStd # UCI
colnames(Cred_int) <- Xnames
rownames(Cred_int) <- c("LCI","UCI") # CI
j_theta <- -OptimRes$hessian
j_theta_inv <- solve(j_theta)
print(OptimRes$hessian) # Mode of Posterior of Beta
print(j_theta_inv) # J(theta)^-1


## 2.2

draw_betas <- function(n, mus, cov_matrix){
  return(rmvnorm(n, mean = mus, sigma = cov_matrix))
}

beta_draws <- draw_betas(1000, approxPostMode, solve((-OptimRes$hessian)))
test_subject <- c(1, 38, 1, 10, 0, 11000)

posterior_predictive_distribution <- 1/(1 + exp(beta_draws %*% test_subject))
p <- ggplot() +
  geom_histogram(aes(x = posterior_predictive_distribution), alpha = 0.5, position = "identity", bins = 20, col = "red", fill = "red") +
  xlab("Probability Value") + 
  ylab("Probabilty Density") +
  ggtitle("Distribution of posterior predictive distribution for a given patient")

print(p)