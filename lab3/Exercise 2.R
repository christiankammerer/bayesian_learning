library(mvtnorm)

# Load data
X <- read.table("eBayNumberOfBidderData_2025.dat", header = TRUE)
y <- X$nBids
X$nBids <- NULL
x <- as.matrix(X)

# Prior parameters
p <- ncol(x)
mu <- rep(0, p)
sigma <- 100 * solve(t(x) %*% x + diag(1e-6, p))  # regularize to avoid singularity

# Stable log prior
log_prior <- function(sigma, mu, beta){
  k <- length(mu)
  term1 <- - (k / 2) * log(2 * pi)
  term2 <- - 0.5 * log(det(sigma))
  term3 <- - 0.5 * t(beta - mu) %*% solve(sigma) %*% (beta - mu)
  return(as.numeric(term1 + term2 + term3))
}

# Stable log likelihood
log_likelihood <- function(x, y, beta){
  eta <- x %*% beta
  eta <- pmin(pmax(eta, -100), 100)  # clip to avoid exp overflow
  return(sum(eta * y - exp(eta) - lgamma(y + 1)))
}

# Log posterior
log_posterior <- function(sigma, mu, x, y, beta){
  return(log_prior(sigma, mu, beta) + log_likelihood(x, y, beta))
}

# Optimization (posterior mode and Hessian)
neg_log_posterior <- function(beta) {
  -log_posterior(sigma, mu, x, y, beta)
}

optres <- optim(par = rep(0, p), fn = neg_log_posterior, hessian = TRUE, method = "BFGS")
beta_mode <- optres$par
neg_hessian <- optres$hessian

# Metropolis algorithm
metropolis_alg <- function(n, initial_vec, sample_func, c = 1, proposal_cov, ...) {
  p <- length(initial_vec)
  samples <- matrix(NA, nrow = n, ncol = p)
  current <- initial_vec
  current_log_post <- sample_func(current, ...)
  samples[1, ] <- current
  
  for (i in 2:n) {
    proposal <- mvtnorm::rmvnorm(1, mean = current, sigma = c * proposal_cov)
    proposal <- as.vector(proposal)
    proposal_log_post <- sample_func(proposal, ...)
    
    log_ratio <- proposal_log_post - current_log_post
    
    if (is.finite(log_ratio) && (log_ratio >= 0 || log(runif(1)) <= log_ratio)) {
      current <- proposal
      current_log_post <- proposal_log_post
    }
    samples[i, ] <- current
  }
  return(samples)
}

# Run sampler
initial_vec <- rep(0, p)
samples <- metropolis_alg(
  n = 10000,
  initial_vec = initial_vec,
  sample_func = log_posterior,
  c = 1,
  proposal_cov = solve(neg_hessian),
  sigma = sigma,
  mu = mu,
  x = x,
  y = y
)
