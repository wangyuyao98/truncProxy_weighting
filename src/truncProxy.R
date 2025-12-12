
# The PQB estimator for marginal survival probability 
# nu the transformation
PQB_estimator <- function(nu, dat, time.name, Q.name, W1.name, W2.name, Z.name, 
                          weights = rep(1, nrow(dat)),
                          trim.min = 1e-7, trim.max = 1e7) {
  
  TT = dat[, time.name]
  Q = dat[, Q.name]
  q_Z = as.matrix(cbind(1, dat[,W1.name], dat[,Z.name]))
  colnames(q_Z) <- c(1, W1.name, Z.name)
  q_W = as.matrix(cbind(1, dat[,W2.name], dat[,Z.name]))
  colnames(q_W) <- c(1, W2.name, Z.name)
  
  v = sort(unique(Q))
  qtime.id = (min(TT) <= v & v <= max(TT))
  qtime = v[qtime.id]
  
  
  # estimate the bridge process b
  qfit = qbridge(TT = TT, Q = Q, qtime = qtime, Z = q_Z, W = q_W, weights = weights)
  B_int = qfit$B_int
  qtime = qfit$qtime
  # B = qfit$B
  # dB = qfit$dB
  # plot(qtime2, dB[,3])
  
  # colnames(B_int) <- qfit$qtime
  
  bT_hat = exp(B_int_eval(TT, B_int, qtime))    ## ?? Some of the weights are less than 1. Is this an issue? -- maybe not
  bT_hat = pmin(pmax(bT_hat, trim.min), trim.max)  # bound the estimated weights to improve stability
  
  numer = mean(weights * bT_hat * nu(TT))
  denom = mean(weights * bT_hat)
  
  return(numer/denom)
}



# t0 is the minimum time poit such that $\nu(T) = \nu(min(T,t_0))$
PQB_IPCW_estimator <- function(nu, t0, dat, time.name, Q.name, event.name = NULL, W1.name, W2.name, Z.name, 
                          weights = rep(1, nrow(dat)),
                          trim.min = 1e-7, trim.max = 1e7,
                          IPCW_time_varying = FALSE, trim.IPCW = 1e-7) {
  
  names = c(time.name, Q.name, W1.name, W2.name, Z.name)
  if(sum(names == "weights")>0){
    stop("The names of the variables cannot be 'weights', which is used in the middle of computation")
  }else if(sum(names == "IPCW_weights")>0){
    stop("The names of the variables cannot be 'IPCW_weights', which is used in the middle of computation")
  }
  
  if(is.null(event.name)){
    if(sum(names == "delta1")>0){
      stop("The names of the variables cannot be 'delta1', which is used in the middle of computation")
    }
    event.name = "delta1"
    dat$delta1 = rep(1, nrow(dat))
    IPCW_weights = rep(1, nrow(dat))
  }
  
  
  # estimate the IPCW weights by fitting KM for D
  Y = dat[, time.name] - dat[, Q.name]
  delta_C = 1 - dat[, event.name]
  fitD = survfit(Surv(Y, delta_C) ~1, type = "kaplan-meier", weights = weights)
  stepfunD = stepfun(fitD$time, c(1, fitD$surv))
  # # Check the KM fit -- it is close to the truth
  # plot(fitD)
  # xx = seq(0, 4, by = 0.1)
  # lines(xx, Sd_truth(xx), col = 2)
  
  if(IPCW_time_varying){
    
    X = dat[, time.name]
    Q = dat[, Q.name]
    q_Z = as.matrix(cbind(1, dat[,W1.name], dat[,Z.name]))
    colnames(q_Z) <- c(1, W1.name, Z.name)
    q_W = as.matrix(cbind(1, dat[,W2.name], dat[,Z.name]))
    colnames(q_W) <- c(1, W2.name, Z.name)
    
    v = sort(unique(Q))
    qtime.id = (min(X) <= v & v <= max(X))
    qtime = v[qtime.id]
    
    # prepare the IPCW weights and the delta_t matrices
    IPCW_weights_mx = matrix(nrow = nrow(dat), ncol = length(qtime))
    delta_t_mx = matrix(nrow = nrow(dat), ncol = length(qtime))
    for(j in 1:length(qtime)){
      temp_time = qtime[j]
      delta_t_mx[,j] = as.numeric(dat[,event.name] == 1 | X > temp_time)
      IPCW_weights_mx[,j] = delta_t_mx[,j]/pmax(stepfunD(pmin(X,temp_time)-Q), trim.IPCW)
    }
    
    
    # estimate the bridge process b
    qfit = qbridge_IPCW(TT = X, Q = Q, qtime = qtime, Z = q_Z, W = q_W, 
                   IPCW_weights = IPCW_weights_mx, weights = weights)
    B_int = qfit$B_int
    qtime = qfit$qtime
    # B = qfit$B
    # dB = qfit$dB
    # plot(qtime2, dB[,3])
    
    # colnames(B_int) <- qfit$qtime
    
    bT_hat = exp(B_int_eval(X, B_int, qtime))    ## ?? Some of the weights are less than 1. Is this an issue? -- maybe not
    bT_hat = pmin(pmax(bT_hat, trim.min), trim.max)  # bound the estimated weights to improve stability
    
    
    # Compute the IPCW weights for the estimator
    
    # if(sum(names == "delta_t0")>0){
    #   stop("The names of the variables cannot be 'delta_t0', which is used in the middle of computation")
    # }
    # delta_t0 = as.numeric(dat[,event.name] == 1 | dat[,time.name] > t0)
    # IPCW_weights_t0 = delta_t0/pmax(stepfunD(pmin(X,t0)-Q), trim.IPCW)
    # 
    # numer = mean(weights * IPCW_weights_t0 * bT_hat * nu(X))
    # denom = mean(weights * IPCW_weights_t0 * bT_hat)
    # est = numer/denom
    
    Q_max = max(Q)
    tau_delta = max(Q_max, t0)
    delta_t = as.numeric(dat[,event.name] == 1 | X > tau_delta)
    IPCW_weights = delta_t/pmax(stepfunD(pmin(X, tau_delta) - Q), trim.IPCW)
    numer = mean(weights * IPCW_weights * bT_hat * nu(X))
    denom = mean(weights * IPCW_weights * bT_hat)
    
    est = numer/denom
    
    
  }else{  # with time invariant IPCW weights
    
    X = dat[, time.name]
    Q = dat[, Q.name]
    
    Q_max = max(Q)
    tau_delta = max(Q_max, t0)
    delta_t = as.numeric(dat[,event.name] == 1 | X > tau_delta)
    IPCW_weights = delta_t/pmax(stepfunD(pmin(X, tau_delta) - Q), trim.IPCW)
    
    est = PQB_estimator(nu = nu, dat = dat, 
                        time.name = time.name, Q.name = Q.name, 
                        W1.name = W1.name, W2.name = W2.name, Z.name = Z.name, 
                        weights = weights * IPCW_weights,
                        trim.min = trim.min, trim.max = trim.max) 
  }
  
  
  return(est)
}




IPQW_estimator_aalen <- function(nu, dat, time.name, Q.name, Z.name, 
                                 weights = rep(1, nrow(dat)),
                                 trim.min = 1e-7, trim.max = 1e7){
  
  TT = dat[, time.name]
  Q = dat[, Q.name]
  q_Z = as.matrix(cbind(1, dat[,Z.name]))
  colnames(q_Z) <- c(1, Z.name)
  # q_W = as.matrix(cbind(1, dat[,Z.name]))
  # colnames(q_W) <- c(1, Z.name)
  
  v = sort(unique(Q))
  qtime.id = (min(TT) <= v & v <= max(TT))
  qtime = v[qtime.id]
  
  
  # estimate the bridge process b
  qfit = qbridge(TT = TT, Q = Q, qtime = qtime, Z = q_Z, W = q_Z, weights = weights)
  B_int = qfit$B_int
  qtime = qfit$qtime
  # B = qfit$B
  # dB = qfit$dB
  # plot(qtime2, dB[,3])
  
  # colnames(B_int) <- qfit$qtime
  
  bT_hat = exp(B_int_eval(TT, B_int, qtime))    ## ?? Some of the weights are less than 1. Is this an issue? -- maybe not
  bT_hat = pmin(pmax(bT_hat, trim.min), trim.max)  # bound the estimated weights to improve stability
  
  numer = mean(weights * bT_hat * nu(TT))
  denom = mean(weights * bT_hat)
  
  return(numer/denom)
}

# t0: the minimum time point such that nu(T) = mu(min(t,t0)). It can take Inf.
IPQW_IPCW_estimator_aalen <- function(nu, t0 = NULL, dat, time.name, Q.name, event.name, Z.name, 
                                 weights = rep(1, nrow(dat)),
                                 trim.min = 1e-7, trim.max = 1e7,
                                 IPCW_time_varying = FALSE, trim.IPCW = 1e-7){
  
  # names = c(time.name, Q.name, W1.name, W2.name, Z.name)
  # if(sum(names == "weights")>0){
  #   stop("The names of the variables cannot be 'weights', which is used in the middle of computation")
  # }else if(sum(names == "IPCW_weights")>0){
  #   stop("The names of the variables cannot be 'IPCW_weights', which is used in the middle of computation")
  # }
  
  # estimate the IPCW weights by fitting KM for D
  Y = dat[, time.name] - dat[, Q.name]
  delta_C = 1 - dat[, event.name]
  fitD = survfit(Surv(Y, delta_C) ~1, type = "kaplan-meier", weights = weights)
  stepfunD = stepfun(fitD$time, c(1, fitD$surv))

  
  if(IPCW_time_varying){
    
    TT = dat[, time.name]
    Q = dat[, Q.name]
    q_Z = as.matrix(cbind(1, dat[,Z.name]))
    colnames(q_Z) <- c(1, Z.name)
    
    v = sort(unique(Q))
    qtime.id = (min(TT) <= v & v <= max(TT))
    qtime = v[qtime.id]
    
    # prepare the IPCW weights and the delta_t matrices
    IPCW_weights_mx = matrix(nrow = nrow(dat), ncol = length(qtime))
    delta_t_mx = matrix(nrow = nrow(dat), ncol = length(qtime))
    for(j in 1:length(qtime)){
      temp_time = qtime[j]
      delta_t_mx[,j] = as.numeric(dat[,event.name] == 1 | TT > temp_time)
      IPCW_weights_mx[,j] = delta_t_mx[,j]/pmax(stepfunD(pmin(TT,temp_time)-Q), trim.IPCW)
    }
    
    # estimate the bridge process b
    qfit = qbridge_IPCW(TT = TT, Q = Q, qtime = qtime, Z = q_Z, W = q_Z, 
                        IPCW_weights = IPCW_weights_mx, weights = weights)
    B_int = qfit$B_int
    qtime = qfit$qtime
    # B = qfit$B
    # dB = qfit$dB
    # plot(qtime2, dB[,3])
    
    # colnames(B_int) <- qfit$qtime
    
    bT_hat = exp(B_int_eval(TT, B_int, qtime))  
    bT_hat = pmin(pmax(bT_hat, trim.min), trim.max)  # bound the estimated weights to improve stability
    
    
    # Compute the IPCW weights for the estimator

    # names = c(time.name, Q.name, W1.name, W2.name, Z.name)
    # if(sum(names == "delta_t0")>0){
    #   stop("The names of the variables cannot be 'delta_t0', which is used in the middle of computation")
    # }
    # delta_t0 = as.numeric(dat[,event.name] == 1 | dat[,time.name] > t0)
    # IPCW_weights_t0 = delta_t0/pmax(stepfunD(pmin(TT,t0)-Q), trim.IPCW)
    # 
    # numer = mean(weights * IPCW_weights_t0 * bT_hat * nu(TT))
    # denom = mean(weights * IPCW_weights_t0 * bT_hat)

    Q_max = max(Q)
    tau_delta = max(Q_max, t0)
    delta_t = as.numeric(dat[,event.name] == 1 | TT > tau_delta)
    IPCW_weights = delta_t/pmax(stepfunD(pmin(TT, tau_delta) - Q), trim.IPCW)
    
    numer = mean(weights * IPCW_weights * bT_hat * nu(TT))
    denom = mean(weights * IPCW_weights * bT_hat)

    est = numer/denom
    
    
  }else{ # with time invariant IPCW weights
    
    X = dat[, time.name]
    Q = dat[, Q.name]
    
    Q_max = max(Q)
    tau_delta = max(Q_max, t0)
    delta_t = as.numeric(dat[,event.name] == 1 | X > tau_delta)
    IPCW_weights = delta_t/pmax(stepfunD(pmin(X, tau_delta) - Q), trim.IPCW)
    
    est = IPQW_estimator_aalen(nu = nu, dat = dat, 
                               time.name = time.name, Q.name = Q.name, Z.name = Z.name, 
                               weights = weights * IPCW_weights,
                               trim.min = trim.min, trim.max = trim.max)
    
  }
  
 
  
  return(est)
}




IPQW_estimator <- function(nu, dat, time.name, Q.name, cov.names,
                           model = "Cox", trim = 1e-7, weights, details = FALSE) {
  
  names = c(time.name, Q.name, cov.names)
  if(sum(names == "delta.1")){
    stop("The names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
  }
  
  # compute the IPQW weights G(T|Z)
  dat$delta.1 = rep(1, nrow(dat))
  jumps.Q = sort(unique(dat[,Q.name]))
  v = c(jumps.Q, max(jumps.Q) + 1e-7)
  Gvz.mx = est_G(dat.fit = dat, dat.est = dat, time.eval = v, model = model, 
                 time.name = time.name, Q.name = Q.name, event.name = "delta.1", cov.names = cov.names,
                 weights = weights, trunc = TRUE)
  
  GTZ = CDF_eval(dat[,time.name], Gvz.mx)
  
  IPQW_weights = 1/pmax(GTZ, trim)
  
  numer <- mean(weights * IPQW_weights * nu(dat[,time.name]))
  denom <- mean(weights * IPQW_weights)
  
  est_IPQW = numer/denom
  
  if(details){
    results = list(estimator = est_IPQW,
                   GTZ = GTZ)
  }else{
    results = est_IPQW
  }
  
  return(results)
}



IPQW_oracle_estimator <- function(nu, dat, time.name, trim = 1e-7, 
                                  weights = rep(1, nrow(dat))) {
  
  
  GTZ = G_truth(dat$Z, dat$U, dat[, time.name])
  
  IPQW_weights = 1/pmax(GTZ, trim)
  
  numer <- mean(weights * IPQW_weights * nu(dat[,time.name]))
  denom <- mean(weights * IPQW_weights)
  
  est_IPQW = numer/denom
  
  return(est_IPQW)
}




## functions for evaluating B_int at given times ----------------------------------
## B_int = B_0(t) + B_1(t)W_1 + B_z(t)Z 

# Return a vector of CDF(time.eval_i|Z_i)
B_int_eval <- function(time.eval, B_int, qtime){
  if(length(time.eval) != nrow(B_int)){
    stop("The number of time points does not equal the number of subjects!")
  }
  
  if(ncol(B_int) != (length(qtime)+1) ){
    stop("ncol(B_int) should equal to length(qtime) + 1, where the last column of B_int correspond to the time in the initial condition.")
  }
  
  probs = rep(NA, length(time.eval))
  for(i in 1:length(time.eval)){
    id = findInterval(time.eval[i], c(0, qtime, Inf))
    probs[i] = B_int[i,id]
  }
  
  return(probs)
}




## functions for computing the CDF at given times ----------------------------------

# Return a vector of CDF(time.eval_i|Z_i)
CDF_eval <- function(time.eval, CDF.mx){
  if(length(time.eval) != nrow(CDF.mx)){
    stop("The number of time points does not equal the number of subjects!")
  }
  
  jumps = as.numeric(colnames(CDF.mx))
  CDF.mx = cbind(0, CDF.mx)
  
  probs = rep(NA, length(time.eval))
  for(i in 1:length(time.eval)){
    id = findInterval(time.eval[i], c(0, jumps, Inf))
    probs[i] = CDF.mx[i,id]
  }
  
  return(probs)
}



# # Return a matrix of CDF(time.eval_j|Z_i) - (i,j)-th element
# CDF_eval.mx <- function(time.eval, CDF.mx){
#   jumps = as.numeric(colnames(CDF.mx))
#   CDF.mx = cbind(0, CDF.mx)
#   id = findInterval(time.eval, c(0, jumps, Inf))
#   probs = CDF.mx[,id]
#   
#   if(length(time.eval)>1 & is.null(dim(probs))){
#     names(probs) = time.eval
#   }else if(length(time.eval)>1 & (!is.null(dim(probs)))){
#     colnames(probs) = time.eval
#   }
#   
#   return(probs)
# }


