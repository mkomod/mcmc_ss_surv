# spike and slab gibbs sampling
library(Rcpp)
sourceCpp("./lpl.cpp")

# test data
set.seed(1)
n <- 250; omega <- 1; censoring_lvl <- 0.2
b <- c(0.5, 2, rep(0, 28)); p <- length(b)
X <- matrix(rnorm(n * p), nrow=n)
y <- runif(nrow(X))
Y <- log(1 - y) / - (exp(X %*% b) * omega)
delta  <- runif(n) > censoring_lvl   	# 0: censored, 1: uncensored
Y[!delta] <- Y[!delta] * runif(sum(!delta))

# set-up indices of failures
Y_sorted <- order(Y) - 1                # we -1 as C++ indices start at 0
Y_failure <- which(delta[Y_sorted] == 1) - 1
b.true <- b                             # save the true beta

# testing cpp funcs
log_PL(X, b, Y_sorted, Y_failure) 	# should be larger
log_PL(X, rnorm(p, 0, 0.2), Y_sorted, Y_failure)
xs <- seq(-4, 4, by=0.1)
ys <- sapply(xs, function(x) log_Laplace(x, 1))
plot(xs, exp(ys))


# Gibbs sampler
lambda <- 1
kernel.sd <- 0.3

# initial values
z <- matrix(sample(c(0, 1), p, replace=T), nrow=p)
b <- matrix(rnorm(p, 0, 0.2), nrow=p)  # init with random beta 
w <- runif(1)

# store results
Z <- matrix(0, nrow=p, ncol=1e4)
B <- matrix(0, nrow=p, ncol=1e4)
W <- matrix(0, nrow=1, ncol=1e4)


system.time(
for (iter in 1:ncol(B)) {

# update w
w <- runif(1)
pz <- matrix(c(w, 1-w), nrow=2)

# update z
for (i in 1:p) {
    z[i] <- 0
    p0 <- log(pz[1]) + log_PL(X, b * z, Y_sorted, Y_failure)
    z[i] <- 1
    p1 <- log(pz[2]) + log_PL(X, b * z, Y_sorted, Y_failure)

    prob0 <- sigmoid(p0 - p1)
    prob1 <- sigmoid(p1 - p0)

    z[i] <- sample(c(0, 1), 1, prob=c(prob0, prob1))
}


# update beta, we use MH
for (i in 1:p) {
    b.old <- b[i]
    b.new <- rnorm(1, b.old, kernel.sd * 10^(1 - z[i])) 	# cancels out by symmetry
    
    # MH numerator
    b[i] <- b.new
    p1 <- log_PL(X, b * z, Y_sorted, Y_failure) + log_Laplace(b.new, lambda) 
	# + dnorm(b.old, b.new, kernel.sd, log=T)
    
    # MH denominator
    b[i] <- b.old
    p2 <- log_PL(X, b * z, Y_sorted, Y_failure) + log_Laplace(b.old, lambda) 
	# + dnorm(b.new, b.old, kernel.sd, log=T)

    A <- min(1, exp(p1 - p2))
    b[i] <- ifelse(A > runif(1), b.new, b.old)
}

# save samples
Z[ , iter] <- z
B[ , iter] <- b
W[ , iter] <- w

} # outer loop
) # timer


# results take 1e3 - 1e4 (burn in of 1e3)
par(mfrow=c(2,4))
plot(B[1, 1e3:1e4], type="l", ylab=expression(beta[1])); plot(density(B[1, 1e3:1e4]), main="")
plot(B[2, 1e3:1e4], type="l", ylab=expression(beta[2])); plot(density(B[2, 1e3:1e4]), main="")
plot(B[3, 1e3:1e4], type="l", ylab=expression(beta[3])); plot(density(B[3, 1e3:1e4]), main="")
plot(B[4, 1e3:1e4], type="l", ylab=expression(beta[4])); plot(density(B[4, 1e3:1e4]), main="")

plot(Z[1, 1e3:1e4], type="l", ylab=expression(beta[1])); hist(Z[1, 1e3:1e4], main="")
plot(Z[2, 1e3:1e4], type="l", ylab=expression(beta[2])); hist(Z[2, 1e3:1e4], main="")
plot(Z[3, 1e3:1e4], type="l", ylab=expression(beta[3])); hist(Z[3, 1e3:1e4], main="")
plot(Z[4, 1e3:1e4], type="l", ylab=expression(beta[4])); hist(Z[4, 1e3:1e4], main="")

# posterior means
apply(B[ , 1e3:1e4], 1, mean)
apply(Z[ , 1e3:1e4], 1, mean)
mean(W)

# posterior sd
apply(B[ , 1e3:1e4], 1, sd)
apply(Z[ , 1e3:1e4], 1, sd)
sd(W)


for(i in 1:8)
    acf(B[i, ])

