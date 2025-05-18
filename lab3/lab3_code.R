library(mvtnorm)
library(BayesLogit)
library(ggplot2)
library(tidyr)
library(dplyr)
library(rstan)

# Exercise 1

# Variables
X <- read.csv("Disease.csv")

# Sampling Functions
beta_prior <- function(n, b, B){
  rmvnorm(n, mean = b, sigma = B)
} 


# Ineffeciency factor relies on the assumption that all auto correlation values are positive
# Geyer's inefficiency therefore terminates the summation when faced with the first
# negative auto-correlation value
compute_if_geyer <- function(x, lag.max = 100) {
  acf_vals <- acf(x, lag.max = lag.max, plot = FALSE)$acf[-1] # compute auto correlations
  m <- floor(length(acf_vals) / 2)
  psum <- acf_vals[seq(1, 2 * m, by = 2)] + acf_vals[seq(2, 2 * m, by = 2)]
  k <- which(psum <= 0)
  if (length(k) > 0) {
    m <- k[1] - 1
  }
  pos_rho <- acf_vals[1:(2 * m)]
  IF <- 1 + 2 * sum(pos_rho)
  return(max(1, IF))
}

run_experiment <- function(X){
  d <- ncol(X)
  n <- nrow(X)
  # Split X and y, add y-intercept to X, type conversion
  y <- X$class_of_diagnosis; X$class_of_diagnosis <- NULL; X <- as.matrix(X); X <- cbind(1, X); y <- as.numeric(y)
  tau <- 3
  I <- diag(d)
  kappa <- y - rep(-1/2, n)
  b <- rep(0, d)
  B <- tau^2*I
  
  # Execution context
  beta_vec <- beta_prior(n = 1, b = b, B = B) # sampling from prior
  
  beta_matrix <- matrix(nrow = n, ncol = d)
  
  # Gibbs sampling
  for(i in 1:n){
    omega_vec <- rpg(n, 1, z = X[i,] %*% beta_vec)
    omega_matrix <- diag(omega_vec)
    V_w <- solve(t(X) %*% omega_matrix %*% X + solve(B))
    m_w <- V_w %*% (t(X) %*% kappa + solve(B) %*% b)
    beta_vec <- rmvnorm(1, mean = m_w, sigma = V_w)
    beta_matrix[i,] <- beta_vec
  }
  
  beta_frame <- data.frame(beta_matrix)
  colnames(beta_frame) <- paste0("Beta", seq_len(ncol(beta_frame)))
  beta_frame$Iteration <- 1:nrow(beta_frame)  # Add iteration column
  beta_long <- pivot_longer(beta_frame, 
                            cols = starts_with("Beta"), 
                            names_to = "Coefficient", 
                            values_to = "Value")
  
  # Plot
  p <- ggplot(beta_long, aes(x = Iteration, y = Value, color = Coefficient)) +
    geom_line(alpha = 0.9) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Gibbs Sampling Trace Plots for beta Coefficients",
      x = "Iteration",
      y = expression(beta),
      color = "Coefficient"
    ) +
    theme(legend.position = "right")
  print(p)
  print("Ineffciency Factors: ")
  print(sapply(beta_frame[1:7], compute_if_geyer))
}

for(m in c(313, 80, 40, 10)){
  X_current <- X[1:m, ]
  print(m)
  run_experiment(X_current)
}

################################################################################

# Exercise 2

# Load data
X <- read.table("eBayNumberOfBidderData_2025.dat", header = TRUE)
y <- X$nBids
X$nBids <- NULL
x <- as.matrix(X)

# Prior parameters
p <- ncol(x)
mu <- rep(0, p)
sigma <- 100 * solve(t(x) %*% x + diag(1e-6, p))  # regularize to avoid singularity




# task b


library(mvtnorm)
set.seed(12345)
X <- read.table("eBayNumberOfBidderData_2025.dat", header = TRUE)
y <- X$nBids
X$nBids <- NULL
x <- as.matrix(X)

# Prior parameters
p <- ncol(x)
mu <- rep(0, p)
sigma <- 100 * solve(t(x) %*% x + diag(1e-6, p)) # Ensure numeric stability

log_prior <- function(sigma, mu, beta){
  k <- nrow(sigma)
  
  return(
    log(2*pi) * k/2 -
      1/2 * log(det(sigma)) -
      1/2 * t(beta - mu) %*% solve(sigma) %*% (beta - mu)
  )
}

log_likelihood <- function(x, y, beta){
  n <- nrow(x)
  return(
    sum(
      sapply((1:n), function(i) t(x[i, ]) %*% beta * y[i]
             - exp(t(x[i, ]) %*% beta) 
             - log(factorial(y[i])))
    )
  )
}

log_posterior <- function(sigma, mu, x, y, beta){
  return(log_prior(sigma, mu, beta) + log_likelihood(x, y, beta))
}

# Objective to minimize (negative log-posterior)
neg_log_posterior <- function(beta) {
  -log_posterior(sigma, mu, x, y, beta)
}

# Optimization (to find posterior mode)
optres <- optim(par = rep(0, p), fn = neg_log_posterior, hessian = TRUE, method = "BFGS")

# Posterior mode
beta_mode <- optres$par

# Negative Hessian at mode (i.e. observed information)
neg_hessian <- optres$hessian

# Output results
cat("Posterior mode (beta):\n")
print(beta_mode)

cat("\nNegative Hessian at mode:\n")
print(neg_hessian)


# task c

metropolis_alg <- function(n, initial_vec, sample_func, c=1, 
                           sigma, ...){
  p <- length(initial_vec)
  sample_matrix <- matrix(nrow = n, ncol = length(initial_vec))
  x_t <- rmvnorm(1, mean = rep(0, p), sigma = sigma)
  f_x_t <- sample_func(sigma, mu, x, y, as.vector(x_t))
  for(i in 1:n){
    proposal <- rmvnorm(1, mean = x_t, sigma = c * sigma)
    f_proposal <- sample_func(sigma, mu, x, y, as.vector(proposal))
    acceptance_prob <- f_proposal - f_x_t # since we are in log-space, we subtract instead of dividing
    if(acceptance_prob >= 1 || (runif(1) <= acceptance_prob)){
      sample_matrix[i, ] <- proposal
      x_t <- proposal
      f_x_t <- f_proposal
    } else{
      sample_matrix[i, ] <- x_t
    }
  }
  return(sample_matrix)
}

samples <- metropolis_alg(n = 100, initial_vector=mp00, beta_mode, sample_func=log_posterior, c=0.0001, sigma=neg_hessian, mu = rep(0, p), x = x, y = y)


matplot(samples, type = "l", lty = 1, col = rainbow(ncol(samples)),
        main = "Trace Plots", ylab = "Parameter value", xlab = "Iteration")


#task d
x_new <- c(1,   # Intercept
           1,   # PowerSeller
           0,   # VerifyID
           1,   # Sealed
           0,   # MinBlem
           1,   # MajBlem
           0,   # LargNeg
           1.3, # LogBook
           0.7) # MinBidShare

after_burnin <- samples[4000:10000, ]

eta <- after_burnin %*% x_new
lambda <- exp(eta)
y_pred <- rpois(length(lambda), lambda)

# Predictive distribution
hist(y_pred, breaks = 30, probability = TRUE,
     main = "Posterior Predictive Distribution",
     xlab = "Number of Bidders")

# Probability of zero bidders
pr_0 <- mean(y_pred == 0)
print(pr_0)


################################################################################


# Exercise 3

# task a

simulate_ar1 <- function(mu, sigma2, phi, T){
  x <- numeric(T)
  x[1] <- mu # use mu as initial value
  
  for (t in 2:T){
    epsilon_t <- rnorm(1,mean=0, sd=sqrt(sigma2))
    x[t] <- mu + phi * (x[t-1] - mu) + epsilon_t
  }
  
  return (x)
}


mu <- 5
sigma2 <- 9
T <- 300
phi_values <- c(-0.9, -0.5, 0, 0.5, 0.9) # arbitrary values from -1 to 1





par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

for (phi in phi_values) {
  x <- simulate_ar1(mu=mu, sigma2=sigma2, phi=phi, T=T)
  plot(x, type = "l", main = paste("AR(1) Simulation, phi =", phi),
       xlab = "Time", ylab = expression(x[t]))
  abline(h = mu, col = "red", lty = 2)
}




# task b


set.seed(42)
x1 <- simulate_ar1(mu=mu, sigma2=sigma2, phi=0.4, T=T)
x2 <- simulate_ar1(mu=mu, sigma2=sigma2, phi=0.98, T=T)

stan_data_x1 <- list(T = T, x = x1)
stan_data_x2 <- list(T = T, x = x2)

stan_model <- "
  data{
    int<lower=2> T;
    vector[T] x;
  }
  
  parameters{
    real mu;
    real<lower=-1, upper=1> phi;
    real<lower=0> sigma;
  }
  
  model {
    mu ~ normal(5,10);
    phi ~ uniform(-1,1);
    sigma ~ cauchy(0,5);
    
    for(t in 2:T) {
      x[t] ~ normal(mu + phi * (x[t-1] - mu), sigma);
    }
  }
"
fit_x1 <- stan(model_code=stan_model, data = stan_data_x1, chains = 4, iter = 2000)
fit_x2 <- stan(model_code=stan_model, data = stan_data_x2, chains = 4, iter = 2000)


print(fit_x1, pars = c("mu", "phi", "sigma"), probs = c(0.025, 0.975))
print(fit_x2, pars = c("mu", "phi", "sigma"), probs = c(0.025, 0.975))


traceplot(fit_x1, pars = c("mu", "phi", "sigma"))
traceplot(fit_x2, pars = c("mu", "phi", "sigma"))



posterior_samples_x1 <- as.data.frame(rstan::extract(fit_x1, pars = c("mu", "phi")))
posterior_samples_x2 <- as.data.frame(rstan::extract(fit_x2, pars = c("mu", "phi")))

ggplot(posterior_samples_x1, aes(x = mu, y = phi)) +
  geom_point(alpha = 0.2) +
  labs(title = "Joint Posterior (Phi = 0.4)", x = expression(mu), y = expression(phi)) +
  theme_minimal()

ggplot(posterior_samples_x2, aes(x = mu, y = phi)) +
  geom_point(alpha = 0.2) +
  labs(title = "Joint Posterior (Phi = 0.98)", x = expression(mu), y = expression(phi)) +
  theme_minimal()
