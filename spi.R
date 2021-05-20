# Gibbs sampler for semi-parametric inference
library(mvtnorm)
Rcpp::sourceCpp("./spi.cpp")

# test data
theta.true <- matrix(2, ncol=1)
f.true <- Vectorize(function(x) ifelse(x > 5 && x <= 7, 2 * x^1.2, 0), "x")

n <- 50
x <- as.matrix(seq(0, 10, length.out=n))
y <- x %*% theta.true + f.true(x) + rnorm(n, 0, 1)

plot(x, y)
x <- t(x)

# Parameters
s <- 1
kern.s <- 0.8
B <- 0.001
b <- 2

# Semi-parametric regression
K <- rbf(x, x, kern.s)
K.y <- (K + s * diag(rep(1, n)))
x.star <- x
K.sx <- rbf(x.star, x, kern.s)
K.xs <- rbf(x, x.star, kern.s)
K.ss <- rbf(x.star, x.star, kern.s)


theta <- solve(solve(B) + x %*% solve(K.y) %*% t(x)) %*% 
	(x %*% solve(K.y) %*% y + solve(B) %*% b)
R <- x.star - x %*% solve(K.y) %*% K.sx
M.f <- t(x.star) %*% theta + K.sx %*% solve(K.y) %*% (y - t(x) %*% theta)
C.f <- K.ss - K.sx %*% solve(K.y) %*% K.xs + t(R) %*% 
    solve(solve(B) + x %*% solve(K.y) %*% t(x)) %*% R

f.s <- rmvnorm(1000, M.f, C.f)
f.m <- apply(f.s, 2, mean)
f.u <- apply(f.s, 2, function(x) quantile(x, probs=0.95))
f.l <- apply(f.s, 2, function(x) quantile(x, probs=0.05))

plot(x, y)
lines(x, f.m)
lines(x, f.u)
lines(x, f.l)

