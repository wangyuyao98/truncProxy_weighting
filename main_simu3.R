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
# Settings used for Table 1 of the paper: 500 simulation replicates, B = 200.
# SIMU_REPS overrides the replicate count for a quick smoke test, e.g.
#   SIMU_REPS=5 Rscript main_simu3.R
rep_num <- as.integer(Sys.getenv("SIMU_REPS", "500")) # number of simulation runs
B <- 200 # number of bootstrap

# Number of worker processes. Set this to control CPU usage; the results are
# identical for any value (each rep re-seeds its own RNG stream and all
# bootstrap weights are pre-drawn from it). Forked parallelism is not
# available on Windows, hence the fallback to 1.
n_cores <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)

# Bootstrap iterations per parallel work unit; ~10 balances scheduling
# overhead against load balancing.
boot_chunk_size <- 10

methods <- c("PQB", "IPQW", "IPQW2", "IPQWu", "PQB_noC", "IPQWo", "naive", "PL", "KM")

# Compute all estimators for one dataset and one weight vector.
compute_all_estimators <- function(dat, weights) {

  # survfit() looks the weights argument up in the formula environment, so
  # the formulas must be created here where `weights` is in scope.
  formula.PL = as.formula(paste0("Surv(", Q.name, ",", time.name, ",", event.name, ") ~ 1"))
  formula.KM = as.formula(paste0("Surv(", time.name, ",", event.name, ") ~ 1"))

  # The censoring KM fit and time-varying IPCW weight matrix are the same for
  # all four IPCW estimators below, so compute them once per weight draw.
  IPCW_prep <- prep_IPCW(dat = dat, time.name = time.name, Q.name = Q.name,
                         event.name = event.name, weights = weights,
                         trim.IPCW = trim.IPCW, time_varying = IPCW_time_varying)

  PQB_est <- PQB_IPCW_estimator(nu = nu, t0 = t00, dat = dat,
                                time.name = time.name, Q.name = Q.name, event.name = event.name,
                                W1.name = W1.name, W2.name = W2.name, Z.name = Z.name,
                                weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW,
                                IPCW_prep = IPCW_prep)

  IPQW_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                        dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                        Z.name = c(Z.name, W1.name, W2.name),
                                        weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW,
                                        IPCW_prep = IPCW_prep)

  IPQW2_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                         dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                         Z.name = c(Z.name, W1.name),
                                         weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW,
                                         IPCW_prep = IPCW_prep)

  IPQWu_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                         dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                         Z.name = c(Z.name, U.name),
                                         weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW,
                                         IPCW_prep = IPCW_prep)

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

# Phase 1 of a rep: generate the data, pre-draw all bootstrap weights from
# the rep's seeded RNG stream (same draw order as a sequential rexp(N) loop,
# so the weights are identical to a serial run), and compute the point
# estimates.
gen_one_rep <- function(rep_seed) {
  set.seed(rep_seed)
  dat.list <- data_gen(N, para_set)
  dat <- dat.list$dat

  weights_mat <- matrix(rexp(B * N), nrow = B, byrow = TRUE)

  list(dat = dat, weights_mat = weights_mat,
       est = compute_all_estimators(dat, rep(1, N)))
}

# Phase 2 work unit: one chunk of bootstrap iterations for one rep.
run_boot_chunk <- function(unit) {
  rep_data <- rep_setups[[unit$rep]]
  vapply(unit$b_idx, function(b) {
    weights <- rep_data$weights_mat[b, ]
    weights <- weights / mean(weights)
    compute_all_estimators(rep_data$dat, weights)
  }, numeric(length(methods)))  # methods x length(b_idx)
}

set.seed(213)
seeds = sample(10^6, rep_num, replace = FALSE)

start_time <- Sys.time()

rep_setups <- mclapply(seeds, gen_one_rep, mc.cores = min(n_cores, rep_num))
failed <- vapply(rep_setups, function(r) inherits(r, "try-error") || is.null(r$est), logical(1))
if (any(failed)) {
  stop("Simulation rep(s) failed: ", paste(which(failed), collapse = ", "), "\n",
       paste(unlist(rep_setups[failed]), collapse = "\n"))
}

# Flatten (rep x bootstrap-chunk) into one work list so all cores stay busy
# even when rep_num is smaller than n_cores.
work_units <- do.call(c, lapply(seq_len(rep_num), function(r) {
  lapply(split(seq_len(B), ceiling(seq_len(B) / boot_chunk_size)),
         function(bs) list(rep = r, b_idx = bs))
}))

chunk_results <- mclapply(work_units, run_boot_chunk, mc.cores = n_cores)
failed <- vapply(chunk_results, function(r) inherits(r, "try-error") || !is.numeric(r), logical(1))
if (any(failed)) {
  bad <- which(failed)
  bad_desc <- vapply(bad, function(i) {
    sprintf("rep %d, bootstrap %s (%s)",
            work_units[[i]]$rep,
            paste(range(work_units[[i]]$b_idx), collapse = "-"),
            paste(as.character(chunk_results[[i]]), collapse = " "))
  }, character(1))
  stop("Bootstrap chunk(s) failed:\n", paste(bad_desc, collapse = "\n"))
}

end_time <- Sys.time()

est_results <- do.call(rbind, lapply(rep_setups, `[[`, "est"))  # rep_num x methods

# Reassemble the bootstrap estimates per rep and compute the bootstrap SEs.
unit_reps <- vapply(work_units, `[[`, integer(1), "rep")
bootSE_results <- do.call(rbind, lapply(seq_len(rep_num), function(r) {
  boot_est <- t(do.call(cbind, chunk_results[unit_reps == r]))  # B x methods
  apply(boot_est, 2, sd)
}))
colnames(bootSE_results) <- methods

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
