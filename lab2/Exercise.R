library(asbio)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(matlib)

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
  theme(legend.position = "none")  # optional, removes legend

print(p)

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
  theme(legend.title = element_blank())
print(p)

p <- hist(sigma_posterior_draws, col =  "red")

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
  ggtitle("Median Regression Curve vs Actual Observations")

print(p)

# Poster distribution of value that minimizes regression curve
posterior_x <- function(betas){
  return(-1/2 * (betas[2]/betas[3]))
}

x_s <- apply(beta_posterior_draws, 1, posterior_x) * 365

get_posterior_parameters_reg <- function(mu0, omega0, sigma_sq_0, v0, y, X) {
  vn <- v0 + length(y)
  XtX <- t(X) %*% X
  print(XtX)
  print(omega0)
  omega_n <- XtX + omega0
  beta_hat <- solve(XtX) %*% t(X) %*% y
  mu_n <- solve(omega_n) %*% t(X) %*% y
  sigma_sq_n <- as.numeric(
    (v0 * sigma_sq_0 + t(y) %*% y - t(mu_n) %*% omega_n %*% mu_n) / vn
  )
  
  return(list(
    mu_n = mu_n,
    omega_n = omega_n,
    sigma_sq_n = sigma_sq_n,
    vn = vn
  ))
}




lambda_posterior_f <- function(n, m, X, y, eta0, lambda0){
  lambda_samples <- rinvchisq(n, eta0, lambda0)
  mu0 <-  rep(0, m + 1)
  sigma_sq_0 <- 1
  v0 <- 1
  
  for(lambda in lambda_samples){
    omega0 <- lambda * diag(m + 1)
    posterior_params <- get_posterior_parameters_reg(mu0 = mu0, 
                                                     omega0 = omega0, 
                                                     sigma_sq_0 =  sigma_sq_0,
                                                     v0 = v0, 
                                                     X = X, 
                                                     y = y)
    lambda_inv <- (eta0 * lambda0)/lambda
    with(posterior_params,
         sqrt(
           abs(omega0)/abs(t(X) %*% X + omega0)
         ) * ((vn * sigma_sq_n)/2)^(-vn/2) * dchisq(lambda_inv,  eta0, lambda0)
  )
  }
}
X <- data$time
X <- cbind(1, X, X^2, X^3, X^4, X^5, X^6, X^7, X^8, X^9, X^10)


  