# Imports
library(mvtnorm)
library(BayesLogit)
library(ggplot2)
library(tidyr)
library(dplyr)

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
      title = "Gibbs Sampling Trace Plots for β Coefficients",
      x = "Iteration",
      y = expression(beta),
      color = "Coefficient"
    ) +
    theme(legend.position = "right")
  print(p)
  
  print(sapply(beta_frame[1:7], compute_if_geyer))
}

for(m in c(313, 10, 40, 80)){
  X_current <- X[1:m, ]
  print(m)
  run_experiment(X_current)
}

