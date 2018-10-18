#' @export
DSTM_EM <- function(Z, Cov0, muinit, Ceta, M, sigma2_eps,
                    H, itermax = 10, tol = 0.1) {
    Cov0 <- as(Cov0, "dsyMatrix")
    Ceta <- as(Ceta, "dsyMatrix")
    Mest <- as(Mest, "dgeMatrix")


    cat("Initialising...\n")
    n <- ncol(H)
    nT <- ncol(Z)
    m <- nrow(Z)
    alpha_filter <- alpha_smooth <- alpha_pred <- matrix(0, n, nT)
    alpha_filter[, 1] <- muinit
    Cov_pred <- Q_pred <- Cov_filter <-  Lfilter <-
        Cov_smooth <- Cov_cross <- list()
    Cov_filter[[1]] <- Cov0
    Lfilter[[1]] <- t(chol(as(Cov0, "dsyMatrix")))
    alpha_resid <- obs_resid <- c()
    I <- Diagonal(n)
    Ceps <- sigma2_eps * Diagonal(m)
    Qeps <- solve(Ceps)
    Leta <- t(chol(Ceta))
    Leps <- t(chol(Ceps))
    L0 <- t(chol(Cov0))
    twonegloglik <- c()
    endloop <- FALSE

    for(iter in 1:itermax) {
        cat(paste0("Iteration ", iter, "\n"))
        cat("...Filtering...\n")

        for(i in 2:nT) {

            ## Forecasting
            alpha_pred[, i] <- Mest %*% alpha_filter[, i - 1] %>% as.numeric()
            Cov_pred[[i]] <- Ceta  + tcrossprod(Mest %*% Lfilter[[i-1]])
            Q_pred[[i]] <- chol2inv(chol(Cov_pred[[i]]))

            ## Filtering
            Rinv <- solve(chol(Q_pred[[i]] + crossprod(forwardsolve(Leps,H))))
            Temp1 <- Qeps - tcrossprod(Qeps %*% H %*% Rinv)
            K <- Cov_pred[[i]] %*% t(H) %*% Temp1
            KH <- K %*% H
            alpha_filter[, i] <- (alpha_pred[,i] +
                                      K %*% Z[,i] - KH %*% alpha_pred[,i]) %>%
                as.numeric()

            Cov_filter[[i]] <- Cov_pred[[i]] -
                crossprod(t(chol(Qeps)) %*% H %*% Cov_pred[[i]]) +
                tcrossprod(Cov_pred[[i]] %*% t(H) %*% Qeps  %*% H %*% Rinv)
            Lfilter[[i]] <- t(chol(Cov_filter[[i]]))
        }

        J <- dHPH <-list()
        alpha_smooth[, nT] <- alpha_filter[, nT]
        Cov_smooth[[nT]] <- Cov_filter[[nT]]
        dHPH[[nT]] <- rowSums((H %*% Cov_smooth[[nT]]) * H)
        alpha_resid[nT] <- crossprod(forwardsolve(Leta,
                                                  alpha_smooth[,nT] - M %*% alpha_smooth[,i]))
        obs_resid[nT] <-  crossprod(forwardsolve(Leps,
                                                 Z[,nT] - H %*% alpha_smooth[,nT]))

        cat("...Smoothing...\n")
        for(i in (nT - 1) : 1) {
            J[[i]] <- Cov_filter[[i]] %*% t(Mest) %*% Q_pred[[i + 1]]
            alpha_smooth[, i] <- (alpha_filter[, i] + J[[i]] %*% (alpha_smooth[, i + 1] - alpha_pred[, i + 1])) %>% as.numeric()
            Cov_smooth[[i]] <- Cov_filter[[i]] -
                tcrossprod( J[[i]] %*% t(chol((Cov_pred[[i+1]] - Cov_smooth[[i+1]]))))
            dHPH[[i]] <- rowSums((H %*% Cov_smooth[[i]]) * H)
            alpha_resid[i] <- crossprod(forwardsolve(Leta,
                                                     alpha_smooth[,i] - M %*% alpha_smooth[,i]))
            obs_resid[i] <-  crossprod(forwardsolve(Leps,
                                                    Z[,i] - H %*% alpha_smooth[,i]))
        }

        cat("...Computing cross-covariances...\n")
        Cov_cross[[nT]] <- (I - KH) %*% Mest %*% Cov_filter[[nT - 1]]
        for(i in nT : 3) {
            Cov_cross[[i - 1]] <- Cov_filter[[i - 1]] %*% J[[i - 2]] +
                J[[i - 1]] %*% (Cov_cross[[i]] - Mest %*% Cov_filter[[i-1]]) %*% J[[i - 2]]
        }

        twonegloglik[iter] <- (2*sum(log(diag(L0))) +
                                   2*nT*sum(log(diag(Leta))) +
                                   2*nT*sum(log(diag(Leps))) +
                                   crossprod(forwardsolve(L0,
                                                          alpha_smooth[, 1] - muinit)) +
                                   sum(alpha_resid) + sum(obs_resid)) %>% as.numeric()

        if(iter == itermax) {
            cat("Maximum iterations reached\n")
            endloop <- TRUE
        }
        if(iter > 1)
            if(abs(twonegloglik[iter] - twonegloglik[iter - 1]) < tol) {
                cat("Tolerance reached\n")
                endloop <- TRUE
            }

        if(endloop) {
            return(list(muinit = muinit,
                        Cov0 = as(Cov0, "matrix"),
                        alpha_smooth = alpha_smooth,
                        Cov_smooth = Cov_smooth,
                        Mest = as(Mest, "matrix"),
                        Ceta = as(Ceta, "matrix"),
                        sigma2_eps = sigma2_eps,
                        twonegloglik = twonegloglik))
        }

        S00 <- Reduce("+", Cov_smooth[-nT]) + tcrossprod(alpha_smooth[,-nT])
        R00 <- chol(S00)
        S11 <- Reduce("+", Cov_smooth[-1]) + tcrossprod(alpha_smooth[,-1])
        S10 <- Reduce("+", Cov_cross[-1]) +
            tcrossprod(alpha_smooth[,-1], alpha_smooth[,-nT])

        ## Mstep
        cat("...M-step...\n")
        muinit <- alpha_smooth[, 1]
        Cov0 <- Cov_smooth[[1]]
        Mest <- S10 %*% chol2inv(R00)
        Ceta <- (S11 - crossprod(forwardsolve(t(R00), t(S10)))) / nT
        sigma2_eps <- (sum((Z - H %*% alpha_smooth)^2) +
                           sum(sapply(dHPH,sum))) / (nT * m)
        Ceps <- sigma2_eps * Diagonal(m)
        Qeps <- solve(Ceps)
        Leta <- t(chol(Ceta))
        Leps <- t(chol(Ceps))
        L0 <- t(chol(Cov0))
        # Probably should double-check equations for negloglik and sigma2eps at some point

    }
}
