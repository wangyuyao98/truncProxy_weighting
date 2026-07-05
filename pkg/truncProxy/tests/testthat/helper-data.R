paper_like_parameters <- function() {
  list(
    mu_Z = 0.6,
    sigma_Z = 0.45,
    mu_U = 0.6,
    sigma_U = 0.45,
    mu_W1 = c(1.4, 0.3, -0.9),
    sigma_W1 = 0.25,
    mu_W2 = c(0.6, -0.2, 0.5),
    sigma_W2 = 0.25,
    mu_Q = c(0.1, 0.25, 1),
    mu_TT = c(0.25, 0.3, 0.6),
    T.min = 0,
    Q.max = 2,
    shape_D = 2,
    scale_D = 2
  )
}

make_example_data <- function(n = 120L, multi = 20L) {
  set.seed(1)
  para_set <- paper_like_parameters()

  Z <- pmax(0, para_set$mu_Z + rnorm(multi * n, 0, para_set$sigma_Z))
  U <- pmax(0, para_set$mu_U + rnorm(multi * n, 0, para_set$sigma_U))

  W1 <- cbind(1, Z, U) %*% para_set$mu_W1 + rnorm(multi * n, 0, para_set$sigma_W1)
  W2 <- cbind(1, Z, U) %*% para_set$mu_W2 + rnorm(multi * n, 0, para_set$sigma_W2)

  TT <- para_set$T.min + rexp(multi * n, cbind(1, Z, U) %*% para_set$mu_TT)

  tau <- para_set$Q.max
  Q2 <- rexp(multi * n, cbind(1, Z, U) %*% para_set$mu_Q)
  Q2 <- pmin(Q2, tau)
  Q <- tau - Q2

  D <- rweibull(multi * n, shape = para_set$shape_D, scale = para_set$scale_D)
  C <- Q + D
  X <- pmin(TT, C)
  delta <- as.integer(TT < C)

  dat_full <- data.frame(X = X, TT = TT, delta = delta, Q = Q, W1 = W1, W2 = W2, Z = Z)
  dat_obs <- dat_full[dat_full$Q < dat_full$TT, , drop = FALSE]
  if (nrow(dat_obs) < n) {
    stop("Truncation rate is high. Increase `multi` in the test helper.")
  }

  dat_obs[seq_len(n), , drop = FALSE]
}
