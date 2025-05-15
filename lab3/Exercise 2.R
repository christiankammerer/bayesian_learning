X <- read.table("eBayNumberOfBidderData_2025.dat", header = TRUE)
# Exercise a)
lm <- glm(nBids ~ . - Const, X, family = "poisson")
print(summary(lm))

# Exercise b)
