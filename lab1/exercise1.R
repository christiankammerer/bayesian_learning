n <- 78
f <- 35
s <- n - f
alpha <- 7 
beta <- 7
drawsize <- 10000

posterior_sample <- function(n, alpha, beta, f, s){
  return(rbeta(n, alpha + s, beta + f))
}

n_draws <- posterior_sample(drawsize, alpha, beta, f, s)
n_means <- numeric(drawsize)
n_sds <- numeric(drawsize)

for(i in 1:drawsize){
  n_means[[i]] <- mean(n_draws[1:i])
  n_sds[[i]] <- sd(n_draws[1:i])
}

n_sds[[0]] <- 0