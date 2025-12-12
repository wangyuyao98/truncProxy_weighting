# Estimate the bridge function

# Z = c(1,W1,Z)
# W = c(1,W2,Z)
# qtime 

bbridge <- function(TT, Q, qtime, Z, W, weights, del = 0.001) {
    
    N <- length(TT)
    
    if(N != length(Q) | N != nrow(Z) | N != nrow(W)){
        stop("The length of TT, Q, Z, W should be the same, all equal to the sample size.")
    }
    
    K_q <- length(qtime)
    dim_Z <- ncol(Z)
    dim_W <- ncol(W)
    
    B_numer <- numeric(dim_W)
    B <- matrix(0, K_q, dim_Z)
    dB <- matrix(0, K_q, dim_Z)
    B_int <- matrix(0, N, K_q)
    dB_int <- matrix(0, N, K_q)
    dN <- matrix(0, N, K_q)
    risk_q <- matrix(0, N, K_q)
    
    for (j in (K_q-1):1) {  # check if here should be K_q:1, and see if the following code can be simplified.
        B_numer[] <- 0
        B_denom <- matrix(0, dim_W, dim_Z)
        
        for (i in 1:N) {
            dN[i, j] <- if (Q[i] == qtime[j]) -1 else 0
            risk_q[i, j] <- if (Q[i] <= qtime[j] && qtime[j] < TT[i]) 1 else 0
            
            for (dim_w in 1:dim_W) {  #!! check if the following j should be j+1
                B_numer[dim_w] <- B_numer[dim_w] + weights[i] * W[i, dim_w] * exp(B_int[i, j]) * dN[i, j]
                for (dim_z in 1:dim_Z) {
                    B_denom[dim_w, dim_z] <- B_denom[dim_w, dim_z] + 
                        weights[i] * W[i, dim_w] * Z[i, dim_z] * risk_q[i, j] * exp(B_int[i, j])
                }
            }
        }
        
        # Call the C++ pinv function here
        pinv_result <- pinv(B_denom, del)
        if ("error" %in% names(pinv_result)) {
            warning(paste("Pseudo-inverse failed at time point", j))
            next
        }
        tmp_pinv <- pinv_result$pinv
        
        # Compute dB
        dB[j, ] <- tmp_pinv %*% B_numer
        
        # Update B_int and B
        for (i in 1:N) {
            dB_int[i, j] <- Z[i, ] %*% dB[j, ]
            B_int[i, j] <- B_int[i, j+1] - dB_int[i, j]
        }
        
        if (j == K_q) {
            B[j, ] <- -dB[j, ]
        } else {
            B[j, ] <- B[j+1, ] - dB[j, ]
        }
    }
    
    return(list(
        dN = dN,
        risk_q = risk_q,
        dB = dB,
        B = B,
        B_int = B_int,
        dB_int = dB_int,
        B_numer = B_numer,
        B_denom = B_denom,
        qtime = qtime
    ))
}
