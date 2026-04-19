# truncProxy

`truncProxy` implements proximal weighting estimators for the expectation of an
arbitrarily transformed event time under dependent left truncation, with
optional inverse probability of censoring weighting to handle right censoring.

The current package exports:

- `PQB_estimator()`
- `PQB_IPCW_estimator()`

The associated paper is available at <https://arxiv.org/pdf/2512.21283>.

## Installation

Once the package is published on CRAN, it can be installed with:

```r
install.packages("truncProxy")
```

During development, it can also be installed from GitHub with:

```r
devtools::install_github("wangyuyao98/truncProxy_weighting", subdir = "pkg/truncProxy")
```

During development from this repository root, it can be installed with:

```r
devtools::install("pkg/truncProxy")
```

## Simulation-Based Example

```r
library(truncProxy)

simulate_truncproxy_data <- function(n = 300, multi = 20) {
  para_set <- list(
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

  Z <- pmax(0, para_set$mu_Z + rnorm(multi * n, 0, para_set$sigma_Z))
  U <- pmax(0, para_set$mu_U + rnorm(multi * n, 0, para_set$sigma_U))
  W1 <- cbind(1, Z, U) %*% para_set$mu_W1 + rnorm(multi * n, 0, para_set$sigma_W1)
  W2 <- cbind(1, Z, U) %*% para_set$mu_W2 + rnorm(multi * n, 0, para_set$sigma_W2)
  TT <- para_set$T.min + rexp(multi * n, cbind(1, Z, U) %*% para_set$mu_TT)

  tau <- para_set$Q.max
  Q2 <- rexp(multi * n, cbind(1, Z, U) %*% para_set$mu_Q)
  Q2 <- pmin(Q2, tau)
  Q <- tau - Q2

  D <- rweibull(n, shape = para_set$shape_D, scale = para_set$scale_D)
  C <- Q + D
  X <- pmin(TT, C)
  delta <- as.integer(TT < C)

  dat_full <- data.frame(X = X, TT = TT, delta = delta, Q = Q, W1 = W1, W2 = W2, Z = Z)
  dat_obs <- dat_full[dat_full$Q < dat_full$TT, , drop = FALSE]
  dat_obs[seq_len(n), , drop = FALSE]
}

set.seed(1)
dat <- simulate_truncproxy_data()
nu <- function(t) as.numeric(t > 1)

PQB_estimator(
  nu = nu,
  dat = dat,
  time.name = "TT",
  Q.name = "Q",
  W1.name = "W1",
  W2.name = "W2",
  Z.name = "Z"
)

PQB_IPCW_estimator(
  nu = nu,
  t0 = 1,
  dat = dat,
  time.name = "X",
  Q.name = "Q",
  event.name = "delta",
  W1.name = "W1",
  W2.name = "W2",
  Z.name = "Z",
  IPCW_time_varying = TRUE
)
```

## Notes

In the current implementation, the IPCW weights are computed from a weighted
Kaplan-Meier estimator on the residual time scale `time - Q`.
