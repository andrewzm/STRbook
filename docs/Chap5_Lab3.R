## Wikle, C. K., Zammit-Mangion, A., and Cressie, N. (2019), 
## Spatio-Temporal Statistics with R, Boca Raton, FL: Chapman & Hall/CRC
## Copyright (c) 2019 Wikle, Zammit-Mangion, Cressie
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

library("ggplot2")
library("STRbook")
library("expm")
library("Matrix")

## ------------------------------------------------------------------------
data("SSTlandmask", package = "STRbook")
data("SSTlonlat", package = "STRbook")
data("SSTdata", package = "STRbook")
delete_rows <- which(SSTlandmask == 1)   # remove land values
SST_Oct97 <- SSTdata[-delete_rows, 334]  # save Oct 1997 SSTs
SSTdata <- SSTdata[-delete_rows, 1:328]  # until April 1997
SSTlonlat$mask <- SSTlandmask            # assign mask to df

## ------------------------------------------------------------------------
Z <- t(SSTdata)                         # data matrix
spat_mean <- apply(SSTdata, 1, mean)    # spatial mean
nT <- ncol(SSTdata)                     # no. of time points
Zspat_detrend <- Z - outer(rep(1, nT),  # detrend data
                           spat_mean)
Zt <- 1/sqrt(nT-1)*Zspat_detrend        # normalize
E <- svd(Zt)                            # SVD

## ------------------------------------------------------------------------
n <- 10

## ------------------------------------------------------------------------
options(digits = 3)

## ------------------------------------------------------------------------
TS <- Zspat_detrend %*% E$v[, 1:n]
summary(colMeans(TS))

## ------------------------------------------------------------------------
tau <- 6
nT <- nrow(TS)
TStplustau <- TS[-(1:tau), ] # TS with first tau time pts removed
TSt <- TS[-((nT-5):nT), ]    # TS with last tau time pts removed

## ------------------------------------------------------------------------
Cov0 <- crossprod(TS)/nT
Covtau <- crossprod(TStplustau,TSt)/(nT - tau)

## ------------------------------------------------------------------------
C0inv <- solve(Cov0)
Mest <- Covtau %*% C0inv
Ceta <- Cov0 - Covtau %*% C0inv %*% t(Covtau)

## ------------------------------------------------------------------------
image(Mest)
image(Ceta)

## ------------------------------------------------------------------------
SSTlonlat$pred <- NA
alpha_forecast <- Mest %*% TS[328, ]

## ------------------------------------------------------------------------
idx <- which(SSTlonlat$mask == 0)
SSTlonlat$curr[idx]  <- as.numeric(E$v[, 1:n] %*% TS[328, ] +
                                       spat_mean)
SSTlonlat$pred[idx]  <- as.numeric(E$v[, 1:n] %*% alpha_forecast +
                                       spat_mean)

## ------------------------------------------------------------------------
SSTlonlat$obs1[idx]  <- SSTdata[, 328]
SSTlonlat$obs2[idx]  <- SST_Oct97

## ------------------------------------------------------------------------
C <- Mest %*% Cov0 %*% t(Mest) + Ceta

## ------------------------------------------------------------------------
SSTlonlat$predse[idx] <-
    sqrt(diag(E$v[, 1:n] %*% C %*% t(E$v[, 1:n])))

## ------------------------------------------------------------------------
g1 <- ggplot(SSTlonlat) +
    geom_tile(aes(lon,lat,fill=curr)) +
    scale_fill_distiller(palette = "Spectral", limits = c(-2,2), name = "degC") +
    theme_bw() + coord_fixed()

g2 <- ggplot(SSTlonlat) +
    geom_tile(aes(lon,lat,fill=pred)) +
    scale_fill_distiller(palette = "Spectral", limits = c(-1.2,1.2), name = "degC") +
    theme_bw() + coord_fixed()


g3 <- ggplot(SSTlonlat) +
        geom_tile(aes(lon,lat,fill=obs1)) +
    scale_fill_distiller(palette = "Spectral", limits = c(-2,2), name = "degC") +
    theme_bw() + coord_fixed()

g4 <- ggplot(SSTlonlat) +
    geom_tile(aes(lon,lat,fill=obs2)) +
    scale_fill_distiller(palette = "Spectral",name = "degC") +
    theme_bw() + coord_fixed()


g5 <- ggplot(SSTlonlat) +
    geom_tile(aes(lon,lat,fill=predse)) +
    scale_fill_distiller(palette = "Spectral",name = "degC") +
    theme_bw() + coord_fixed()

## ------------------------------------------------------------------------
## DSTM_Results <- DSTM_EM(Z = SSTdata,
##                         Cov0 = Cov0,
##                         muinit = matrix(0, n, 1),
##                         M = Mest,
##                         Ceta = Ceta,
##                         sigma2_eps = 0.1,
##                         H = H <- E$v[, 1:n],
##                         itermax = 10,
##                         tol = 1)

## ------------------------------------------------------------------------
data("DSTM_EM_results", package = "STRbook")

## ------------------------------------------------------------------------
par(mfrow = c(1,3))
for(i in 1:3) {
  plot(DSTM_Results$alpha_smooth[i, ], type = 'l',
       xlab = "t", ylab = bquote(alpha[.(i)]))
  lines(TS[, i], lty = 'dashed', col = 'red')
}
par(mfrow = c(1,3))
for(i in 1:3) {
  plot(DSTM_Results$alpha_smooth[i,], type = 'l',
       xlab = "t", ylab = bquote(alpha[.(i)]))
  lines(TS[,i], lty = 'dashed', col = 'red')
}

## ------------------------------------------------------------------------
image(as(DSTM_Results$Mest, "Matrix"))
image(as(DSTM_Results$Mest %^% 6, "Matrix"))
image(as(Mest, "Matrix"))

## ------------------------------------------------------------------------
alpha <- DSTM_Results$alpha_smooth[, nT]
P <- DSTM_Results$Cov0
for(t in 1:6) {
   alpha <- DSTM_Results$Mest %*% alpha
   P <- DSTM_Results$Mest %*% P %*% t(DSTM_Results$Mest) +
       DSTM_Results$Ceta
}

