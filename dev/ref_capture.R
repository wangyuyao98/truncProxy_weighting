# Capture reference outputs from the implementation currently in the working tree.
# Run from the repo root:  Rscript <path>/ref_capture.R <output.rda>
suppressMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(survival)
})

out_file <- commandArgs(trailingOnly = TRUE)[1]

source("src/data_generating.R")
source("src/truncProxy.R")
source("src/est_nuisance.R")
Rcpp::sourceCpp("src/est_bridge_safe.cpp")

para_set <- para_set_simu3
t00 <- 1
nu <- function(t, t0 = t00) as.numeric(t > t0)
trim.IPCW <- 0.05

set.seed(20260705)
dat.list <- data_gen(1000, para_set)
dat <- dat.list$dat
N <- nrow(dat)

set.seed(999)
w_boot <- rexp(N)
w_boot <- w_boot / mean(w_boot)
w_unif <- rep(1, N)

ref <- list()

for (wn in c("unif", "boot")) {
  weights <- if (wn == "unif") w_unif else w_boot

  ref[[paste0("PQB_noC_", wn)]] <- PQB_estimator(
    nu = nu, dat = dat, time.name = "TT", Q.name = "Q",
    W1.name = "W1", W2.name = "W2", Z.name = "Z", weights = weights)

  ref[[paste0("PQB_tv_", wn)]] <- PQB_IPCW_estimator(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    W1.name = "W1", W2.name = "W2", Z.name = "Z", weights = weights,
    IPCW_time_varying = TRUE, trim.IPCW = trim.IPCW)

  ref[[paste0("PQB_ti_", wn)]] <- PQB_IPCW_estimator(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    W1.name = "W1", W2.name = "W2", Z.name = "Z", weights = weights,
    IPCW_time_varying = FALSE, trim.IPCW = trim.IPCW)

  ref[[paste0("IPQW_tv_", wn)]] <- IPQW_IPCW_estimator_aalen(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    Z.name = c("Z", "W1", "W2"), weights = weights,
    IPCW_time_varying = TRUE, trim.IPCW = trim.IPCW)

  ref[[paste0("IPQW_ti_", wn)]] <- IPQW_IPCW_estimator_aalen(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    Z.name = c("Z", "W1", "W2"), weights = weights,
    IPCW_time_varying = FALSE, trim.IPCW = trim.IPCW)

  ref[[paste0("IPQWu_tv_", wn)]] <- IPQW_IPCW_estimator_aalen(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    Z.name = c("Z", "U"), weights = weights,
    IPCW_time_varying = TRUE, trim.IPCW = trim.IPCW)

  ref[[paste0("IPQWo_", wn)]] <- IPQW_oracle_estimator(
    nu = nu, dat = dat, time.name = "TT", weights = weights)
}

## Raw bridge outputs (bootstrap weights, same inputs PQB_estimator builds)
TT <- dat$TT
Q <- dat$Q
q_Z <- as.matrix(cbind(1, dat$W1, dat$Z))
q_W <- as.matrix(cbind(1, dat$W2, dat$Z))
v <- sort(unique(Q))
qtime <- v[min(TT) <= v & v <= max(TT)]

qfit <- qbridge(TT = TT, Q = Q, qtime = qtime, Z = q_Z, W = q_W, weights = w_boot)
ref$qbridge_B_int <- qfit$B_int
ref$qbridge_dB <- qfit$dB
ref$qbridge_B <- qfit$B
ref$B_int_eval <- B_int_eval(TT, qfit$B_int, qfit$qtime)

## Raw qbridge_IPCW outputs (censored time scale, as PQB_IPCW_estimator builds them)
X <- dat$X
v2 <- sort(unique(Q))
qtime2 <- v2[min(X) <= v2 & v2 <= max(X)]
fitD <- survfit(Surv(X - Q, 1 - dat$delta) ~ 1, type = "kaplan-meier", weights = w_boot)
stepfunD <- stepfun(fitD$time, c(1, fitD$surv))
IPCW_mx <- matrix(nrow = N, ncol = length(qtime2))
for (j in seq_along(qtime2)) {
  d_j <- as.numeric(dat$delta == 1 | X > qtime2[j])
  IPCW_mx[, j] <- d_j / pmax(stepfunD(pmin(X, qtime2[j]) - Q), trim.IPCW)
}
qfit2 <- qbridge_IPCW(TT = X, Q = Q, qtime = qtime2, Z = q_Z, W = q_W,
                      IPCW_weights = IPCW_mx, weights = w_boot)
ref$qbridge_IPCW_B_int <- qfit2$B_int
ref$qbridge_IPCW_dB <- qfit2$dB

## est_G + CDF_eval + baseS chain (Cox model), via IPQW_estimator with details
ipqw_cox <- IPQW_estimator(nu = nu, dat = dat, time.name = "TT", Q.name = "Q",
                           cov.names = c("Z", "W1"), model = "Cox",
                           weights = w_unif, details = TRUE)
ref$IPQW_cox_est <- ipqw_cox$estimator
ref$IPQW_cox_GTZ <- ipqw_cox$GTZ

## quick timing baseline for one expensive call
ref$time_PQB_tv_sec <- system.time(
  PQB_IPCW_estimator(nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q",
                     event.name = "delta", W1.name = "W1", W2.name = "W2", Z.name = "Z",
                     weights = w_boot, IPCW_time_varying = TRUE, trim.IPCW = trim.IPCW)
)[["elapsed"]]

save(dat, w_boot, w_unif, ref, file = out_file)
cat("Saved", length(ref), "reference objects to", out_file, "\n")
cat("Scalar estimates:\n")
for (nm in names(ref)) {
  if (is.numeric(ref[[nm]]) && length(ref[[nm]]) == 1) {
    cat(sprintf("  %-16s %.12f\n", nm, ref[[nm]]))
  }
}
