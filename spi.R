# Gibbs sampler for semi-parametric inference
library(mvtnorm)
Rcpp::sourceCpp("./spi.cpp")

# test data
theta.true <- matrix(2, ncol=1)
f.true <- Vectorize(function(x) ifelse(x > 5 && x <= 7, 2 * x, 0), "x")

n <- 50
x <- as.matrix(seq(0, 10, length.out=n))
y <- x %*% theta.true + f.true(x) + rnorm(n, 0, 1)

plot(x, y)
x <- t(x)

# Parameters
s <- 1
kern.s <- 0.8
B <- 5
b <- 0

# Semi-parametric regression
K <- rbf(x, x, kern.s)
K.y <- (K + s * diag(rep(1, n)))
x.star <- seq(-1, 11, length.out=100)
K.sx <- rbf(x.star, x, kern.s)
K.xs <- rbf(x, x.star, kern.s)
K.ss <- rbf(x.star, x.star, kern.s)


theta <- solve(solve(B) + x %*% solve(K.y) %*% t(x)) %*% 
	(x %*% solve(K.y) %*% y + solve(B) %*% b)
R <- x.star - x %*% solve(K.y) %*% K.xs
M.f <- x.star %*% theta + K.sx %*% solve(K.y) %*% (y - t(x) %*% theta)
C.f <- K.ss - K.sx %*% solve(K.y) %*% K.xs + t(R) %*% 
    solve(solve(B) + x %*% solve(K.y) %*% t(x)) %*% R

f.s <- rmvnorm(1000, M.f, C.f)
f.m <- apply(f.s, 2, mean)
f.u <- apply(f.s, 2, function(x) quantile(x, probs=0.995))
f.l <- apply(f.s, 2, function(x) quantile(x, probs=0.005))

pdf("spi.pdf", width=6, height=5)
plot(x, y, pch=20, col="darkorchid1", las=1)
lines(x.star, f.m, lwd=1.2, lty=1)
polygon(c(x.star, rev(x.star)), c(f.l, rev(f.u)), col=rgb(.75, .25, 1, alpha=0.2), border=NA)
dev.off()

