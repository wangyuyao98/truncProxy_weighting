# Re-run all estimators with the CURRENT working-tree code on the saved
# reference dataset and compare against the frozen baseline.
# Run from the repo root:  Rscript <path>/check_equivalence.R <ref_outputs_v0.rda>
suppressMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(survival)
})

ref_file <- commandArgs(trailingOnly = TRUE)[1]
load(ref_file)  # dat, w_boot, w_unif, ref

source("src/data_generating.R")
source("src/truncProxy.R")
source("src/est_nuisance.R")
Rcpp::sourceCpp("pkg/truncProxy/src/est_bridge_safe.cpp")

para_set <- para_set_simu3
t00 <- 1
nu <- function(t, t0 = t00) as.numeric(t > t0)
trim.IPCW <- 0.05
N <- nrow(dat)

TOL <- 1e-8
n_fail <- 0
check <- function(name, new_val) {
  old_val <- ref[[name]]
  d <- max(abs(as.numeric(new_val) - as.numeric(old_val)))
  ok <- is.finite(d) && d <= TOL
  if (!ok) n_fail <<- n_fail + 1
  cat(sprintf("  %-22s max|diff| = %.3e  %s\n", name, d, if (ok) "OK" else "** FAIL **"))
}

# IPQW_oracle_estimator gained an explicit para_set argument on this branch;
# fall back to the old signature if running against pre-fix code.
oracle <- function(weights) {
  if ("para_set" %in% names(formals(IPQW_oracle_estimator))) {
    IPQW_oracle_estimator(nu = nu, dat = dat, time.name = "TT",
                          para_set = para_set, weights = weights)
  } else {
    IPQW_oracle_estimator(nu = nu, dat = dat, time.name = "TT", weights = weights)
  }
}

cat("Scalar estimators:\n")
for (wn in c("unif", "boot")) {
  weights <- if (wn == "unif") w_unif else w_boot

  check(paste0("PQB_noC_", wn), PQB_estimator(
    nu = nu, dat = dat, time.name = "TT", Q.name = "Q",
    W1.name = "W1", W2.name = "W2", Z.name = "Z", weights = weights))

  check(paste0("PQB_tv_", wn), PQB_IPCW_estimator(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    W1.name = "W1", W2.name = "W2", Z.name = "Z", weights = weights,
    IPCW_time_varying = TRUE, trim.IPCW = trim.IPCW))

  check(paste0("PQB_ti_", wn), PQB_IPCW_estimator(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    W1.name = "W1", W2.name = "W2", Z.name = "Z", weights = weights,
    IPCW_time_varying = FALSE, trim.IPCW = trim.IPCW))

  check(paste0("IPQW_tv_", wn), IPQW_IPCW_estimator_aalen(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    Z.name = c("Z", "W1", "W2"), weights = weights,
    IPCW_time_varying = TRUE, trim.IPCW = trim.IPCW))

  check(paste0("IPQW_ti_", wn), IPQW_IPCW_estimator_aalen(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    Z.name = c("Z", "W1", "W2"), weights = weights,
    IPCW_time_varying = FALSE, trim.IPCW = trim.IPCW))

  check(paste0("IPQWu_tv_", wn), IPQW_IPCW_estimator_aalen(
    nu = nu, t0 = t00, dat = dat, time.name = "X", Q.name = "Q", event.name = "delta",
    Z.name = c("Z", "U"), weights = weights,
    IPCW_time_varying = TRUE, trim.IPCW = trim.IPCW))

  check(paste0("IPQWo_", wn), oracle(weights))
}

cat("Raw bridge outputs:\n")
TT <- dat$TT
Q <- dat$Q
q_Z <- as.matrix(cbind(1, dat$W1, dat$Z))
q_W <- as.matrix(cbind(1, dat$W2, dat$Z))
v <- sort(unique(Q))
qtime <- v[min(TT) <= v & v <= max(TT)]

qfit <- qbridge(TT = TT, Q = Q, qtime = qtime, Z = q_Z, W = q_W, weights = w_boot)
check("qbridge_B_int", qfit$B_int)
check("qbridge_dB", qfit$dB)
check("qbridge_B", qfit$B)
check("B_int_eval", B_int_eval(TT, qfit$B_int, qfit$qtime))

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
check("qbridge_IPCW_B_int", qfit2$B_int)
check("qbridge_IPCW_dB", qfit2$dB)

cat("est_G / CDF_eval chain:\n")
ipqw_cox <- IPQW_estimator(nu = nu, dat = dat, time.name = "TT", Q.name = "Q",
                           cov.names = c("Z", "W1"), model = "Cox",
                           weights = w_unif, details = TRUE)
check("IPQW_cox_est", ipqw_cox$estimator)
check("IPQW_cox_GTZ", ipqw_cox$GTZ)

if (n_fail > 0) {
  cat(sprintf("\n%d CHECK(S) FAILED (tolerance %.0e)\n", n_fail, TOL))
  quit(status = 1)
} else {
  cat(sprintf("\nAll checks passed at tolerance %.0e\n", TOL))
}
