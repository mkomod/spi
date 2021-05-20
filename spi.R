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


# prior parameters
# likelihood
s <- 1
# beta
m.0 <- 0
s.0 <- 1
# f
K <- rbf(x, x, 0.1)
K.y <- K + s * diag(length(x))
C.f <- K - K %*% solve(K.y) %*% K

# init values
theta <- rnorm(1)
f <- t(rmvnorm(1, rep(0, n), K))

N <- 1e3
THETA <- matrix(0, nrow=1, ncol=N)
FUN <- matrix(0, nrow=n, ncol=N)

# sample beta
for (iter in 1:N) {
    # update theta
    x.bar <- mean(x %*% theta + f)
    s.1 <- (1/s.0^2 + n/s^2)^(-1/2)
    m.1 <- s.1^2 * (m.0/s.0^2 + n*x.bar/s^2)
    theta <- rnorm(1, m.1, s.1)

    # update f
    # M.f <- (x %*% theta + f) + K %*% solve(K.y) %*% (y - (x %*% theta + f))
    for (i in 1:N) {
	M.f <- x.bar[i] %*% theta 
	+ K[ , i] %*% solve(K.y) %*% y
	f[i] <- t(rmvnorm(1, M.f, C.f))
    }

    # save
    THETA[ , iter] <- theta
    FUN[ , iter] <- f
}

plot(density(THETA[1, ]))
mean(THETA[1, ])
matplot(FUN, type="l")

