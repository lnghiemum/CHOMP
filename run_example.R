# An example of fitting CHOMP and Adaptive CHOMP using PIC to select the tuning parameter
library(MASS)
n = 500; p = 100
# Generating X ~ N_p(0, Sigma)
set.seed(2000)
cov_sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
omega <- diag(sqrt(runif(p, 0.5, 4)))
SIGMA <- omega %*% cov_sigma %*%  omega
X = mvrnorm(n, mu = rep(0,p), Sigma = SIGMA)

# Example 1: Single index model
s <- 5
beta <- rep(0, p)
beta[1:s] <- sample(c(-1,1), s, replace=TRUE)*runif(s, 1, 1.5)
y <- exp(X %*% beta) + rnorm(n)
fit.AdaptCHOMP_WeightPower1 <- AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 1)
fit.AdaptCHOMP_WeightPower2 <- AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 2)

# Example 2: Double index model
# Generating sparse beta and outcome from a double index model
set.seed(200)
s <- list(s1 = 1:4, s2 = 3:6) # each element of s corresponds to the index of non-zero coefficient in each dimension
beta = matrix(0, p, length(s))
for (j in 1:length(s)){
  beta[s[[j]], j] = sample(c(-1,1), length(s[[j]]), replace=TRUE)*runif(length(s[[j]]), 1, 1.5)
}
Z = X %*% beta
y = (Z[, 1]) * exp(Z[,2] + rnorm(n))   

fit.AdaptChOMP <-  AdaptCHOMPwithPIC(X, y, d = 2)