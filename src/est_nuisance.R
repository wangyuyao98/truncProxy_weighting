# function for fitting a model for G on the reversed time scale

#' @title Estimate the Conditional CDF for the Left Truncation Time given Covariates
#' @description Estimate the conditional cumulative distribution function (CDF) of the left truncation time given covariates evaluated at given time points. The options implemented in this function are: Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen'.
#' @param dat.fit data frame that is used to fit the model for the full data conditional distribution of the event time given the covariates.
#' @param dat.est data frame that contains the subjects for which the estimated conditional CDF is computed.
#' @param time.eval vector of time points at which the conditional CDF is evaluated.
#' @param model method used to estimate the conditional CDF. The options available are "Cox" and "survPen", corresponding to Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen', respectively.
#' @param time.name name of the event time variable.
#' @param Q.name name of the left truncation time variable.
#' @param event.name name of the event indicator.
#' @param cov.names vector of the names of covariates.
#' @param trim constant for bounding the estimated conditional CDF from 0.
#' @param weights vector of case weights.
#' @param trunc a logical indicator indicating whether there is truncation or not. If \code{trunc = TRUE}, then the reversed time technique is applied.
#' @param formula.survPen the formula when applying the hazard model with penalized splines implemented in \code{survPen::survPen}.
#' @return \code{G_est()} returns a matrix of the estimated conditional CDF for subjects in \code{data.est} evaluated at the time points in the vector \code{time.eval}. Each row corresponds to a subject and each column corresponds to a time point. The column names of the matrix are the times in `\code{time.eval}'.
#' @export
#' @import survival stats survPen
#' @seealso  \code{\link{F_est}}
#' @examples
#' data("simu")
#' v = c(0.5, 1, 1.5, 2, 2.5, 3)
#' Gvz.mx = est_G(simu, simu[1:10,], v, "Cox", "time", "Q", "delta", c("Z1","Z2"))
est_G <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                  time.name, Q.name, event.name, cov.names, trim = 0, weights = NULL,
                  formula.survPen = NULL, tau = NULL, trunc = TRUE,
                  x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                  nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, df = 5){
    
    if(mean(dat.fit[,event.name])<1){
        stop("The truncation time is always observed, so dat.fit[,event.name] should be a vector of one's.")
    }
    
    names = c(time.name, Q.name, event.name, cov.names)
    if(sum(names == "Q2")){
        stop("The names of the variables cannot be 'Q2'. It is used in the middle of the computation.")
    }
    if(sum(names == "T2")){
        stop("The names of the variables cannot be 'T2'. It is used in the middle of the computation.")
    }
    
    v = time.eval
    
    if(is.null(tau)){
        tau = max(c(dat.fit[,time.name], dat.est[,time.name])) + 1
    }
    
    dat.fit$Q2 = tau - dat.fit[ ,Q.name]
    dat.fit$T2 = tau - dat.fit[ ,time.name]
    dat.fit$delta.1 = rep(1, nrow(dat.fit))
    
    dat.est$Q2 = tau - dat.est[ ,Q.name]
    dat.est$T2 = tau - dat.est[ ,time.name]
    dat.est$delta.1 = rep(1, nrow(dat.est))
    
    # if(model == "RF"){
    # 
    #     G_hat = G_hat.ltrcrrf(dat.fit, time.name, Q.name, event.name, cov.names,
    #                           mtry, ntree, tau)
    #     Gvz.mx = pmax(G_hat(newdata = dat.est, time.eval = v), trim)   #***
    # 
    # }else
    if(model == "Cox"){
        
        cov.names = c(cov.names, cov.names.binary)
        id.cox = ((dat.fit[, time.name] - dat.fit[, Q.name]) >= 10^{-7})
        dat.cox = dat.fit[id.cox, ]
        Z.Q = as.matrix(dat.est[,cov.names])
        
        if(trunc){
            formula.Q2 = formula(paste("Surv(tau-", time.name, ", tau-", Q.name,",", event.name, ") ~ ",
                                       paste(cov.names, collapse = " + "), collapse = ""))
        }else{
            formula.Q2 = formula(paste("Surv(Q2,", event.name, ") ~ ",
                                       paste(cov.names, collapse = " + "), collapse = ""))
        }
        
        
        # fit Cox-PH model for (tau-Q)|Z
        if(is.null(weights)){
            fit.Q2 = coxph(formula.Q2, data = dat.cox)
        }else{
            fit.Q2 = coxph(formula.Q2, data = dat.cox, weights = weights[id.cox])
        }
        
        basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
        beta.Q2 = coef(fit.Q2)
        
        Gvz.mx = Gz(v, Z.Q, basehaz.Q2, beta.Q2, tau)
    
        
    }else if(model == "addhaz"){
        
        if(trunc){
            formula.Q2 = formula(paste0(
                "Surv(T2, Q2,", event.name, ") ~ ",
                paste(paste0("const(", cov.names, ")"), collapse = " + ")
            ))
            
        }else{
            stop("This option has not been implemented.")
        }
        
        if(is.null(weights)){
            weights = rep(1, nrow(dat.fit))
        }
        
        fit.Q2 = aalen(formula.Q2, data = dat.fit, weights = weights)
        
        
    }else{
        stop("This Q model is not implemented in this function!")
    }
    
    colnames(Gvz.mx) = v
    Gvz.mx = pmax(Gvz.mx, trim)
    
    return(Gvz.mx)
    # return(list(Gvz.mx = Gvz.mx, b.Q = beta.Q2))
}










### Functions needed when using Cox models ---------------------------------------
# function for computing the baseline CDF or survival function
baseCDF.single <- function(t, basehaz){
    if(t<min(basehaz[,'time'])){
        return(0)   # check 1 or 0
    }else{
        id = which((t - basehaz[,'time'])>=0)
        index = which.min((t - basehaz[,'time'])[id])
        
        return(1-exp(-basehaz[index,'hazard']))
    }
}

baseCDF <- function(t, basehaz){
    cdf = sapply(t, baseCDF.single, basehaz = basehaz)
    return(cdf)
}

baseS.single <- function(t, basehaz){
    if(t<min(basehaz[,'time'])){
        return(1)   # check 1 or 0
    }else{
        id = which((t - basehaz[,'time'])>=0)
        index = which.min((t - basehaz[,'time'])[id])
        
        return(exp(-basehaz[index,'hazard']))
    }
}

baseS <- function(t, basehaz){
    S = sapply(t, baseS.single, basehaz = basehaz)
    return(S)
}


## function for computing G(t|z)
Gz <- function(t, Z, basehaz.Q2, beta.Q2, tau){
    if(is.matrix(Z)){
        result = matrix(nrow = nrow(Z), ncol = length(t))
        for(j in 1:length(t)){
            if(tau-t[j]<min(basehaz.Q2$time)){
                result[,j] = 1
            }else{
                result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(Z %*% beta.Q2)
            }
        }
    }else{
        result = matrix(nrow = 1, ncol = length(t))
        for(j in 1:length(t)){
            if(tau-t[j]<min(basehaz.Q2$time)){
                result[,j] = 1
            }else{
                result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(sum(Z*beta.Q2))
            }
        }
    }
    colnames(result) = t
    
    return(result)
}

## function for computing F(t|z)
Fz <- function(t, Z, basehaz.T, beta.T){
    if(is.matrix(Z)){
        result = matrix(nrow = nrow(Z), ncol = length(t))
        for(j in 1:length(t)){
            result[,j] = 1-(baseS(t[j], basehaz.T))^exp(Z %*% beta.T)
        }
    }else{
        result = matrix(nrow = 1, ncol = length(t))
        for(j in 1:length(t)){
            result[,j] = 1-(baseS(t[j], basehaz.T))^exp(sum(Z*beta.T))
        }
    }
    colnames(result) = t
    
    return(result)
}


