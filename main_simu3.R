rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(Rcpp)
library(RcppArmadillo)
library(survival)  
source("src/data_generating.R")
source("src/truncProxy.R")

Rcpp::sourceCpp("src/est_bridge_safe.cpp")

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

alpha = 0.05
qz = qnorm(alpha/2, lower = F)

N <- 1000
# rep = 1
rep_num <- 5 # number of simulation runs
B <- 200 # number of bootstrap

PQB_est_results <- rep(NA, rep_num)
IPQW_est_results <- rep(NA, rep_num)
IPQW2_est_results <- rep(NA, rep_num)
IPQWu_est_results <- rep(NA, rep_num)
PQB_noC_est_results <- rep(NA, rep_num)
IPQWo_est_results <- rep(NA, rep_num)
naive_est_results <- rep(NA, rep_num)
PL_est_results <- rep(NA, rep_num)
KM_est_results <- rep(NA, rep_num)

PQB_bootSE_results <- rep(NA, rep_num)
IPQW_bootSE_results <- rep(NA, rep_num)
IPQW2_bootSE_results <- rep(NA, rep_num)
IPQWu_bootSE_results <- rep(NA, rep_num)
PQB_noC_bootSE_results <- rep(NA, rep_num)
IPQWo_bootSE_results <- rep(NA, rep_num)
naive_bootSE_results <- rep(NA, rep_num)
PL_bootSE_results <- rep(NA, rep_num)
KM_bootSE_results <- rep(NA, rep_num)

set.seed(213)
seeds = sample(10^6, rep_num, replace = FALSE)

start_time <- Sys.time()

for(rep in 1:rep_num) {
  set.seed(seeds[rep])
  dat.list <- data_gen(N, para_set)
  dat <- dat.list$dat
  dat.full = dat.list$dat.full
  
  # Compute the estimators
  weights <- rep(1, N)
  weights <- weights/mean(weights)
  
  PQB_est <- PQB_IPCW_estimator(nu = nu, t0 = t00, dat = dat,
                             time.name = time.name, Q.name = Q.name, event.name = event.name,
                              W1.name = W1.name, W2.name = W2.name, Z.name = Z.name,
                              weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
  PQB_est_results[rep] <- PQB_est


  IPQW_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                        dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                        Z.name = c(Z.name, W1.name, W2.name),
                                        weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
  IPQW_est_results[rep] <-  IPQW_est

  IPQW2_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                         dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                         Z.name = c(Z.name, W1.name),
                                         weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
  IPQW2_est_results[rep] <- IPQW2_est

  IPQWu_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                         dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                         Z.name = c(Z.name, U.name),
                                         weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
  IPQWu_est_results[rep] <- IPQWu_est


  # # Oracle estimators
  PQB_noC_est = PQB_estimator(nu = nu, dat = dat,
                              time.name = "TT", Q.name = Q.name,
                              W1.name = W1.name, W2.name = W2.name, Z.name = Z.name,
                              weights = weights)
  PQB_noC_est_results[rep] <- PQB_noC_est

  IPQWo_est <- IPQW_oracle_estimator(nu = nu, dat = dat, time.name = "TT",
                                     weights = weights)
  IPQWo_est_results[rep] <- IPQWo_est

  # Naive estimator
  naive_est = mean(nu(dat[,time.name]) * weights)
  naive_est_results[rep] <- naive_est
  
  # PL estimator that assumes random LTRC
  formula.PL = as.formula(paste0("Surv(", Q.name,",",time.name,",", event.name, ") ~ 1"))
  PLfit = survfit(formula.PL, data = dat, weights = weights)
  surv_PL = stepfun(PLfit$time, c(1, PLfit$surv))
  PL_est_results[rep]  <-  surv_PL(t00)
  
  # KM estimator - KM that account for random censoring, but ignore left truncation
  formula.KM = as.formula(paste0("Surv(", time.name, ",", event.name, ") ~ 1"))
  KMfit = survfit(formula.KM, data = dat, weights = weights)
  surv_KM = stepfun(KMfit$time, c(1, KMfit$surv))
  KM_est_results[rep]  <-  surv_KM(t00)
  
  
  # Bayesian boostrap
  boot_PQB_est_results <- rep(NA, B)
  boot_IPQW_est_results <- rep(NA, B)
  boot_IPQW2_est_results <- rep(NA, B)
  boot_IPQWu_est_results <- rep(NA, B)
  boot_PQB_noC_est_results <- rep(NA, B)
  boot_IPQWo_est_results <- rep(NA, B)
  boot_naive_est_results <- rep(NA, B)
  boot_PL_est_results <- rep(NA, B)
  boot_KM_est_results <- rep(NA, B)
  for(b in 1:B){
      weights <- rexp(N)
      weights <- weights/mean(weights)

      PQB_est <- PQB_IPCW_estimator(nu = nu, t0 = t00, dat = dat,
                                    time.name = time.name, Q.name = Q.name, event.name = event.name,
                                    W1.name = W1.name, W2.name = W2.name, Z.name = Z.name,
                                    weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
      boot_PQB_est_results[b] <- PQB_est


      IPQW_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                            dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                            Z.name = c(Z.name, W1.name, W2.name),
                                            weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
      boot_IPQW_est_results[b] <-  IPQW_est

      IPQW2_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                             dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                             Z.name = c(Z.name, W1.name),
                                             weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
      boot_IPQW2_est_results[b] <- IPQW2_est

      IPQWu_est <- IPQW_IPCW_estimator_aalen(nu = nu, t0 = t00,
                                             dat = dat, time.name = time.name, Q.name = Q.name, event.name = event.name,
                                             Z.name = c(Z.name, U.name),
                                             weights = weights, IPCW_time_varying = IPCW_time_varying, trim.IPCW = trim.IPCW)
      boot_IPQWu_est_results[b] <- IPQWu_est
    
      
      # Oracle estimators
      PQB_noC_est = PQB_estimator(nu = nu, dat = dat,
                                  time.name = "TT", Q.name = Q.name,
                                  W1.name = W1.name, W2.name = W2.name, Z.name = Z.name,
                                  weights = weights)
      boot_PQB_noC_est_results[b] <- PQB_noC_est

      IPQWo_est <- IPQW_oracle_estimator(nu = nu, dat = dat, time.name = "TT",
                                         weights = weights)
      boot_IPQWo_est_results[b] <- IPQWo_est

      # Naive estimator
      naive_est = mean(nu(dat[,time.name]) * weights)
      boot_naive_est_results[b] <- naive_est
      
      # PL estimator that assumes random LTRC
      PLfit = survfit(formula.PL, data = dat, weights = weights)
      surv_PL = stepfun(PLfit$time, c(1, PLfit$surv))
      boot_PL_est_results[b]  <-  surv_PL(t00)
      
      # KM estimator - KM that account for random censoring, but ignore left truncation
      formula.KM = as.formula(paste0("Surv(", time.name, ",", event.name, ") ~ 1"))
      KMfit = survfit(formula.KM, data = dat, weights = weights)
      surv_KM = stepfun(KMfit$time, c(1, KMfit$surv))
      boot_KM_est_results[b]  <-  surv_KM(t00)
  }

  PQB_bootSE_results[rep] <- sd(boot_PQB_est_results)
  # PQB_tv_bootSE_results[rep] <- sd(boot_PQB_tv_est_results)
  IPQW_bootSE_results[rep] <- sd(boot_IPQW_est_results)
  IPQW2_bootSE_results[rep] <- sd(boot_IPQW2_est_results)
  IPQWu_bootSE_results[rep] <- sd(boot_IPQWu_est_results)
  PQB_noC_bootSE_results[rep] <- sd(boot_PQB_noC_est_results)
  IPQWo_bootSE_results[rep] <- sd(boot_IPQWo_est_results)
  naive_bootSE_results[rep] <- sd(boot_naive_est_results)
  PL_bootSE_results[rep] <- sd(boot_PL_est_results)
  KM_bootSE_results[rep] <- sd(boot_KM_est_results)
  
  print(rep)
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


PQB_cover_results = (PQB_est_results - qz * PQB_bootSE_results <= truth & truth <= PQB_est_results + qz * PQB_bootSE_results)
IPQW_cover_results = (IPQW_est_results - qz * IPQW_bootSE_results <= truth & truth <= IPQW_est_results + qz * IPQW_bootSE_results)
IPQW2_cover_results = (IPQW2_est_results - qz * IPQW2_bootSE_results <= truth & truth <= IPQW2_est_results + qz * IPQW2_bootSE_results)
IPQWu_cover_results = (IPQWu_est_results - qz * IPQWu_bootSE_results <= truth & truth <= IPQWu_est_results + qz * IPQWu_bootSE_results)

PQB_noC_cover_results = (PQB_noC_est_results - qz * PQB_noC_bootSE_results <= truth & truth <= PQB_noC_est_results + qz * PQB_noC_bootSE_results)
IPQWo_cover_results = (IPQWo_est_results - qz * IPQWo_bootSE_results <= truth & truth <= IPQWo_est_results + qz * IPQWo_bootSE_results)
naive_cover_results = (naive_est_results - qz * naive_bootSE_results <= truth & truth <= naive_est_results + qz * naive_bootSE_results)
PL_cover_results = (PL_est_results - qz * PL_bootSE_results <= truth & truth <= PL_est_results + qz * PL_bootSE_results)
KM_cover_results = (KM_est_results - qz * KM_bootSE_results <= truth & truth <= KM_est_results + qz * KM_bootSE_results)


# save the results
# save(list = ls(),
#      file = paste0("results/", setting, "_n", N, "_B", B, "_tv", IPCW_time_varying, "_t0_", t00, ".rda"))

# save(list = ls(),
#      file = paste0("results/PQB_", setting, "_n", N, "_B", B, "_tv", IPCW_time_varying, ".rda"))

# save(list = ls(),
#      file = paste0("results/PQB_", setting, "_n", N, "_B", B, "_PLKM.rda"))



##### Organize the results into tables -----------------------------------------------------

# List of methods you used
methods <- c("PQB", "IPQW", "IPQW2", "IPQWu", 
             "PQB_noC", "IPQWo", "naive", "PL", "KM")  # Add any additional methods here

# Initialize storage for summary table
summary_table <- data.frame(
  Method = character(),
  Bias = numeric(),
  SD = numeric(),
  BootSE = numeric(),
  CP = numeric(),
  stringsAsFactors = FALSE
)

# Loop through methods and compute summaries
for (method in methods) {
  
  # Load estimates and bootSE
  est <- get(paste0(method, "_est_results"))
  bootSE <- get(paste0(method, "_bootSE_results"))
  cover <- get(paste0(method, "_cover_results"))
  
  # Calculate summary statistics
  bias <- mean(est) - truth
  empirical_sd <- sd(est)
  mean_bootse <- mean(bootSE)
  CP = mean(cover)
  
  
  # Append to summary table
  summary_table <- rbind(summary_table, data.frame(
    Method = method,
    Bias = bias,
    SD = empirical_sd,
    bootSE = mean_bootse,
    CP = CP
  ))
}

# summary_table

# function for round the results into certain digits
format_round <- function(X){
  X[,2] = sprintf("%.4f", X[,2])
  X[,3] = sprintf("%.4f", X[,3])
  X[,4] = sprintf("%.4f", X[,4])
  X[,5] = sprintf("%.3f", X[,5])
  
  return(X)
}

summary_table = format_round(summary_table)
# summary_table

library(xtable)
print(xtable(summary_table,
             type = "latex", 
             align = rep("c", 6),
             file = "tab.tex"),
      include.rownames = FALSE)




