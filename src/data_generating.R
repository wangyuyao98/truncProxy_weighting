
# simu1
para_set_simu1 <- list(mu_Z = 1,
                 sigma_Z = 0.35,
                 mu_U = 1,
                 sigma_U = 0.35,
                 mu_W1 = c(0, 0.3, 0.8),
                 sigma_W1 = 0.25,
                 mu_W2 = c(0, 0.5, 0.9),
                 sigma_W2 = 0.25,
                 mu_Q = c(0.1, 0.35, 0.9),
                 mu_TT = c(0.25, 0.3, 0.7),
                 T.min = 1,
                 Q.max = 2
)

# simu2
para_set_simu2 <- list(mu_Z = 0.6,
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
                       Q.max = 2
)


# simu3
para_set_simu3 <- list(mu_Z = 0.6,
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



data_gen <- function(N, para_set, multi = 20) {
  
  # generate Z*,U*
  Z <- pmax(0, para_set$mu_Z + rnorm(multi*N, 0, para_set$sigma_Z))
  U <- pmax(0, para_set$mu_U + rnorm(multi*N, 0, para_set$sigma_U))

  # generate W1*
  W1 <- cbind(1, Z, U) %*% para_set$mu_W1 + rnorm(multi*N, 0, para_set$sigma_W1)
  
  # generate W2*
  W2 <- cbind(1, Z, U) %*% para_set$mu_W2 + rnorm(multi*N, 0, para_set$sigma_W2)

  # generate T*
  TT <- para_set$T.min + rexp(multi*N, cbind(1, Z, U) %*% para_set$mu_TT)
  
  # generate Q*
  tau = para_set$Q.max
  Q2.max = tau 
  # U2 = runif(multi*N, min = 0, max = 1)
  Q2 = rexp(multi*N, cbind(1, Z, U) %*% para_set$mu_Q)
  Q2 = pmin(Q2, Q2.max)
  Q = tau - Q2
  
  if(is.null(para_set$scale_D)){
    D = rep(Inf, N)
  }else{
    if(is.null(para_set$shape_D)){
      para_set$shape_D = 1
    }
    
    D = rweibull(N, shape = para_set$shape_D, scale = para_set$scale_D)
    # D = rexp(N, rate = para_set$rate_D)
  }
  
  C = Q + D
  X = pmin(TT, C)
  delta = as.numeric(TT < C)
  
  dat.full = data.frame(Z, U, W1, W2, Q, TT, C, X, delta)
  dat.all = dat.full
  obs.id = which(dat.full$Q < dat.full$TT)
  dat.obs = dat.full[obs.id,]
  if(length(obs.id)<N){
    stop('Truncation rate is high. Need to increase `multi`.')
  }
  dat = dat.obs[1:N,]
  dat.full = dat.full[(1:obs.id[N]), ]
  
  
  return(list(dat = dat,
              dat.full = dat.full,
              dat.all = dat.all))
}



# The function for computing the true nuisance parameters
# return a vector of G(time.eval[i]|Z_i,U_i)
G_truth <- function(Z, U, time.eval){
  
  if(length(Z) != length(time.eval) | length(U) != length(time.eval)){
    stop("The length of Z, U, time.eval do not match. ")
  }
  
  lambda =  cbind(1, Z, U) %*% para_set$mu_Q
  Gtz.vec = pexp(para_set$Q.max - time.eval, lambda, lower = F)
  
  return(Gtz.vec)
}

Sd_truth <- function(time.eval){
  
  Sd.vec = pweibull(time.eval, shape = para_set$shape_D, scale = para_set$scale_D, lower = F)
  # Sd.vec = pexp(time.eval, rate = para_set$rate_D, lower = F)
  
  return(Sd.vec)
}





