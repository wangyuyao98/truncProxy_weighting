# Main simulation for setting simu3.
# Run from the repository root, e.g.:  Rscript main_simu3.R

library(Rcpp)
library(RcppArmadillo)
library(survival)
library(parallel)

if (!file.exists("src/truncProxy.R")) {
  stop("Please run this script from the repository root.")
}
source("src/data_generating.R")
source("src/truncProxy.R")

Rcpp::sourceCpp("pkg/truncProxy/src/est_bridge_safe.cpp")

setting = "simu3"

if(setting == "simu1"){
  para_set = para_set_simu1
}else if(setting == "simu2"){
  para_set = para_set_simu2
}else if(setting == "simu3"){
  para_set = para_set_simu3
}

# transformation nu
t00 = 1
nu <- function(t,t0=t00){
  # indicator function
  result = as.numeric(t>t0)

  return(result)
}

alpha = 0.05
qz = qnorm(alpha/2, lower = F)


# # The truth
# data.list <- data_gen(10^5, para_set, multi = 100)
# df = data.list$dat.all
# truth <- mean(nu(df$TT))
# truth  # simu1: truth = 0.54 with t00 = 1.5; simu2: truth = 0.463 with t00 = 1
# TR = mean(df$Q >= df$TT)
# TR
# CR.full = 1 - mean(df$delta)
# CR = 1 - mean(df$delta[df$Q < df$TT])
# rm(data.list, df)
load(paste0("results/", setting, "_truth.rda"))
TR
CR
CR.full



time.name = "X"
event.name = "delta"
Q.name = "Q"
W1.name = "W1"
W2.name = "W2"
Z.name = "Z"
U.name = "U"

trim.IPCW = 0.05
IPCW_time_varying = TRUE   # Whether to use time-varying IPCW weights or case weights --  TRUE: time varying weights, FALSE: case weights

N <- 1000
rep_num <- 5 # number of simulation runs
B <- 200 # number of bootstrap

# Number of worker processes for the simulation reps. Each rep re-seeds its
# own RNG stream, so the results are identical for any n_cores (including 1).
# mclapply forks are not available on Windows, hence the fallback.
n_cores <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)
n_cores <- min(n_cores, rep_num)

methods <- c("PQB", "IPQW", "IPQW2", "IPQWu", "PQB_noC", "IPQWo", "naive", "PL", "KM")

# Compute all estimators for one dataset and one weight vector.
compute_all_estimators <- function(dat, weights) {

  # survfit() looks the weights argument up in the formula environment, so
  # the formulas must be created here where `weights` is in scope.
  formula.PL = as.formula(paste0("Surv(", Q.name, ",", time.name, ",", event.name, ") ~ 1"))
  formula.KM = as.formula(paste0("Surv(", time.name, ",", event.name, ") ~ 1"))

  PQB_est <- PQB_IPCW_estimator(nu = nu, t0 = t00, dat = dat,
                                time.name = time.name, Q.name = Q.name, event.name = event.name,
                                W1.name = W1.name, W2.name = W2.name, Z.name = Z.name,
                                weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)

  IPQW_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                        dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                        Z.name = c(Z.name, W1.name, W2.name),
                                        weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)

  IPQW2_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                         dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                         Z.name = c(Z.name, W1.name),
                                         weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)

  IPQWu_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                         dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                         Z.name = c(Z.name, U.name),
                                         weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)

  # Oracle estimators
  PQB_noC_est = PQB_estimator(nu = nu, dat = dat,
                              time.name = "TT", Q.name = Q.name,
                              W1.name = W1.name, W2.name = W2.name, Z.name = Z.name,
                              weights = weights)

  IPQWo_est <- IPQW_oracle_estimator(nu = nu, dat = dat, time.name = "TT",
                                     para_set = para_set, weights = weights)

  # Naive estimator
  naive_est = mean(nu(dat[,time.name]) * weights)

  # PL estimator that assumes random LTRC
  PLfit = survfit(formula.PL, data = dat, weights = weights)
  surv_PL = stepfun(PLfit$time, c(1, PLfit$surv))
  PL_est <- surv_PL(t00)

  # KM estimator - KM that account for random censoring, but ignore left truncation
  KMfit = survfit(formula.KM, data = dat, weights = weights)
  surv_KM = stepfun(KMfit$time, c(1, KMfit$surv))
  KM_est <- surv_KM(t00)

  c(PQB = PQB_est, IPQW = IPQW_est, IPQW2 = IPQW2_est, IPQWu = IPQWu_est,
    PQB_noC = PQB_noC_est, IPQWo = IPQWo_est, naive = naive_est,
    PL = PL_est, KM = KM_est)
}

# One simulation rep: generate data, point estimates, Bayesian bootstrap SE.
run_one_rep <- function(rep_seed) {
  set.seed(rep_seed)
  dat.list <- data_gen(N, para_set)
  dat <- dat.list$dat

  est <- compute_all_estimators(dat, rep(1, N))

  boot_est <- matrix(NA_real_, nrow = B, ncol = length(methods),
                     dimnames = list(NULL, methods))
  for (b in 1:B) {
    weights <- rexp(N)
    weights <- weights / mean(weights)
    boot_est[b, ] <- compute_all_estimators(dat, weights)
  }

  list(est = est, bootSE = apply(boot_est, 2, sd))
}

set.seed(213)
seeds = sample(10^6, rep_num, replace = FALSE)

start_time <- Sys.time()
rep_results <- mclapply(seeds, run_one_rep, mc.cores = n_cores)
end_time <- Sys.time()

failed <- vapply(rep_results, function(r) inherits(r, "try-error") || is.null(r$est), logical(1))
if (any(failed)) {
  stop("Simulation rep(s) failed: ", paste(which(failed), collapse = ", "), "\n",
       paste(unlist(rep_results[failed]), collapse = "\n"))
}

est_results <- do.call(rbind, lapply(rep_results, `[[`, "est"))        # rep_num x methods
bootSE_results <- do.call(rbind, lapply(rep_results, `[[`, "bootSE"))  # rep_num x methods

elapsed_time <- end_time - start_time
print(elapsed_time)

cover_results <- (est_results - qz * bootSE_results <= truth) &
  (truth <= est_results + qz * bootSE_results)

summary_df <- data.frame(
  method = methods,
  bias = colMeans(est_results) - truth,
  empSE = apply(est_results, 2, sd),
  bootSE = colMeans(bootSE_results),
  coverage = colMeans(cover_results)
)
print(summary_df, digits = 3)


# save the results
# save(est_results, bootSE_results, cover_results, summary_df, truth,
#      N, B, rep_num, seeds, IPCW_time_varying, t00, trim.IPCW,
#      file = paste0("results/", setting, "_n", N, "_B", B, "_tv", IPCW_time_varying, "_t0_", t00, ".rda"))
