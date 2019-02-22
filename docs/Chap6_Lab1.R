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

library("FRK")
library("IDE")
library("scoringRules")
library("verification")
library("dplyr")
library("ggplot2")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

## ------------------------------------------------------------------------
data("radar_STIDF", package = "STRbook")
mtot <- length(radar_STIDF)
radar_STIDF$timeHM <- format(time(radar_STIDF), "%H:%M")

## ------------------------------------------------------------------------
valblock_idx <- which(radar_STIDF$timeHM %in% c("09:35",
                                                "09:45"))
obs_idx <- setdiff(1:mtot, valblock_idx)

## ------------------------------------------------------------------------
set.seed(1)
valrandom_idx <- sample(obs_idx,
                        0.1 * length(obs_idx),
                        replace = FALSE) %>% sort()
obs_idx <- setdiff(obs_idx, valrandom_idx)

## ------------------------------------------------------------------------
radar_obs <- radar_STIDF[obs_idx, ]
radar_valblock <- radar_STIDF[valblock_idx, ]
radar_valrandom <- radar_STIDF[valrandom_idx, ]

## ------------------------------------------------------------------------
## IDEmodel <- IDE(f = z ~ 1,
##                 data = radar_obs,
##                 dt = as.difftime(10, units = "mins"),
##                 grid_size = 41)
## 
## fit_results_radar2 <- fit.IDE(IDEmodel,
##                               parallelType = 1)

## ------------------------------------------------------------------------
data("IDE_Radar_results2", package = "STRbook")

## ------------------------------------------------------------------------
options(digits = 3)

## ------------------------------------------------------------------------
## load results with full data set
data("IDE_Radar_results", package = "STRbook")
with(fit_results_radar$IDEmodel, c(get("betahat")[1,1],
                                   unlist(get("k")),
                                   get("sigma2_eps"),
                                   get("sigma2_eta")))

with(fit_results_radar2$IDEmodel, c(get("betahat")[1,1],
                                   unlist(get("k")),
                                   get("sigma2_eps"),
                                   get("sigma2_eta")))

## ------------------------------------------------------------------------
pred_IDE_block <- predict(fit_results_radar2$IDEmodel,
                          newdata = radar_valblock)
pred_IDE_random <- predict(fit_results_radar2$IDEmodel,
                          newdata = radar_valrandom)

## ------------------------------------------------------------------------
G_spatial <- auto_basis(manifold = plane(),     # fns on plane
                        data = radar_obs,       # project
                        nres = 2,               # 2 res.
                        type = "bisquare",      # bisquare.
                        regular = 1)            # irregular

## ------------------------------------------------------------------------
t_grid <- matrix(seq(0, 12, length = 5))
G_temporal <- local_basis(manifold = real_line(), # fns on R1
                          type = "bisquare",      # bisquare
                          loc = t_grid,           # centroids
                          scale = rep(3.5, 5))    # aperture par.

## ------------------------------------------------------------------------
G <- TensorP(G_spatial, G_temporal)  # take the tensor product

## ------------------------------------------------------------------------
BAUs <- auto_BAUs(manifold = STplane(),   # ST field on plane
                  type = "grid",          # gridded (not "hex")
                  data = radar_obs,       # data
                  cellsize = c(1.65, 2.38, 10), # BAU cell size
                  nonconvex_hull = FALSE, # convex boundary
                  convex = 0,             # no hull extension
                  tunit = "mins")         # time unit is "mins"
BAUs$fs = 1       # fs variation prop. to 1

## ------------------------------------------------------------------------
sigma2_eps <- fit_results_radar2$IDEmodel$get("sigma2_eps")
radar_obs$std <- sqrt(sigma2_eps)

## ------------------------------------------------------------------------
S <- FRK(f = z ~ 1,
         BAUs = BAUs,
         data = list(radar_obs), # (list of) data
         basis = G,           # basis functions
         n_EM = 2,            # max. no. of EM iterations
         tol = 0.01)          # tol. on log-likelihood

## ------------------------------------------------------------------------
FRK_pred <- predict(S)

## ------------------------------------------------------------------------
df_block_over <- over(radar_valblock, FRK_pred)
df_random_over <- over(radar_valrandom, FRK_pred)

## ------------------------------------------------------------------------
radar_valblock_df <- radar_valblock %>%
             data.frame() %>%
             mutate(FRK_pred = df_block_over$mu,
                    FRK_predse = df_block_over$sd,
                    FRK_predZse = sqrt(FRK_predse^2 +
                                       sigma2_eps),
                    IDE_pred = pred_IDE_block$Ypred,
                    IDE_predse = pred_IDE_block$Ypredse,
                    IDE_predZse = sqrt(IDE_predse^2 +
                                       sigma2_eps))

## ------------------------------------------------------------------------
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,0,0))
conditional.quantile(radar_valblock_df$IDE_pred, 
                     radar_valblock_df$z, main = NULL,  xaxs="i", yaxs="i")

## ------------------------------------------------------------------------
radar_valblock_df_long <- radar_valblock_df %>%
                          dplyr::select(s1, s2, timeHM, z,
                                 FRK_pred, IDE_pred) %>%
                          gather(type, num, FRK_pred, IDE_pred)

## ------------------------------------------------------------------------
radar_valrandom_df <- radar_valrandom %>%
                      data.frame() %>%
                      mutate(FRK_pred = df_random_over$mu,
                         FRK_predse = df_random_over$sd,
                         FRK_predZse = sqrt(FRK_predse^2 + S@Ve[1,1] + S@sigma2fshat),
                         IDE_pred = pred_IDE_random$Ypred,
                         IDE_predse = pred_IDE_random$Ypredse,
                         IDE_predZse = sqrt(IDE_predse^2 + sigma2_eps))

radar_valrandom_df_long <- radar_valrandom_df %>%
                           dplyr::select(s1, s2, timeHM, z, FRK_pred, IDE_pred) %>%
                           gather(type, num, FRK_pred, IDE_pred)

## ------------------------------------------------------------------------
ggplot(radar_valblock_df_long) +
   geom_histogram(aes(z - num, fill = type),
                  binwidth = 2, position = "identity",
                  alpha = 0.4, colour = 'black') + theme_bw()

## ------------------------------------------------------------------------
g1 <- ggplot(radar_valblock_df_long) +
    geom_histogram(aes(z - num, fill = type),
                   binwidth = 2, position = "identity",
                   alpha = 0.4, colour = 'black') + theme_bw() +
    xlab(expression(Z - hat(Z)[p])) + ggtitle("Missing time points")
g2 <- ggplot(radar_valrandom_df_long) +
         geom_histogram(aes(z - num, fill = type),
                        binwidth = 2, position = "identity",
                        alpha = 0.4, colour = 'black') + theme_bw() +
    xlab(expression(Z - hat(Z)[p])) + ggtitle("Missing at random")

## ------------------------------------------------------------------------
 ggplot(radar_valblock_df) + geom_point(aes(z, FRK_pred))
 ggplot(radar_valblock_df) + geom_point(aes(z, IDE_pred))

## ------------------------------------------------------------------------
MSPE_time <- rbind(radar_valrandom_df_long,
                   radar_valblock_df_long) %>%
         group_by(timeHM, type) %>%
         dplyr::summarise(MSPE = mean((z - num)^2))

## ------------------------------------------------------------------------
ggplot(MSPE_time) +
 geom_line(aes(timeHM, MSPE, colour = type, group = type))


## ------------------------------------------------------------------------
ggplot(rbind(radar_valrandom_df_long, radar_valblock_df_long)) +
 geom_tile(aes(s1, s2, fill= z - num)) +
 facet_grid(type ~ timeHM) + coord_fixed() +
 fill_scale(name = "dBZ") +
 theme_bw()


## ------------------------------------------------------------------------
Bias <- function(x,y) mean(x - y)          # x: val. obs.
PCV <- function(x,y) mean((x - y)^2)       # y: predictions
SCV <- function(x,y,v) mean((x - y)^2 / v) # v: pred. variances

## ------------------------------------------------------------------------
## Compute CRPS. s is the pred. standard error
CRPS <- function(x, y, s) verification::crps(x, cbind(y, s))$CRPS

## ------------------------------------------------------------------------
Diagblock <- radar_valblock_df %>% summarise(
    Bias_FRK = Bias(FRK_pred, z),
    Bias_IDE = Bias(IDE_pred, z),
    PCV_FRK =  PCV(FRK_pred, z),
    PCV_IDE =  PCV(IDE_pred, z),
    SCV_FRK =  SCV(FRK_pred, z, FRK_predZse^2),
    SCV_IDE =  SCV(IDE_pred, z, IDE_predZse^2),
    CRPS_FRK = CRPS(z, FRK_pred, FRK_predZse),
    CRPS_IDE = CRPS(z, IDE_pred, IDE_predZse)
)

## ------------------------------------------------------------------------
Diagrandom <- radar_valrandom_df %>% summarise(
    Bias_FRK = Bias(FRK_pred, z),
    Bias_IDE = Bias(IDE_pred, z),
    PCV_FRK =  PCV(FRK_pred, z),
    PCV_IDE =  PCV(IDE_pred, z),
    SCV_FRK =  SCV(FRK_pred, z, FRK_predZse^2),
    SCV_IDE =  SCV(IDE_pred, z, IDE_predZse^2),
    CRPS_FRK = CRPS(z, FRK_pred,FRK_predZse),
    CRPS_IDE = CRPS(z, IDE_pred,IDE_predZse)
)
tab <- xtable::xtable(rbind(Diagblock, Diagrandom))
rownames(tab) <- c("Missing time points", "Missing at random")
tab <- strsplit(print(tab),"\n")[[1]]
idx <- which(grepl("Bias", tab))
tab[idx[1]+3] <- strsplit(tab[idx[1]+3], "\\\\")[[1]][1]
fileConn <- file("Chap6_table_results.tex")
writeLines(tab[(idx[1]+2):(idx[1] + 3)], fileConn)
close(fileConn)

## ------------------------------------------------------------------------
radar_val0935 <- subset(radar_valblock,
                        radar_valblock$timeHM == "09:35")
n_0935 <- length(radar_val0935)  # number of validation data

## ------------------------------------------------------------------------
pred_IDE_block <- predict(fit_results_radar2$IDEmodel,
                          newdata = radar_val0935,
                          covariances = TRUE)

## ------------------------------------------------------------------------
FRK_pred_block <- predict(S,
                          newdata = radar_val0935,
                          covariances = TRUE)

## ------------------------------------------------------------------------
Veps <- diag(rep(sigma2_eps, n_0935))

## ------------------------------------------------------------------------
L_IDE <- t(chol(pred_IDE_block$Cov + Veps))
L_FRK <- t(chol(FRK_pred_block$Cov + Veps))

## ------------------------------------------------------------------------
IntIDE <- coef(fit_results_radar2$IDEmodel)
IntFRK <- coef(S)

## ------------------------------------------------------------------------
nsim <- 100
E <- matrix(rnorm(n_0935*nsim), n_0935, nsim)
Sims_IDE <- IntIDE + pred_IDE_block$newdata$Ypred + L_IDE %*% E
Sims_FRK <- IntFRK + FRK_pred_block$newdata$mu + L_FRK %*% E

## ------------------------------------------------------------------------
## Put into long format
radar_val0935_long <- cbind(data.frame(radar_val0935),
                            IDE = Sims_IDE[,1],
                            FRK = Sims_FRK[,1]) %>%
                      gather(type, val, z, FRK, IDE)

## Plot
gsims <- ggplot(radar_val0935_long) +
    geom_tile(aes(s1, s2, fill = val)) +
  facet_grid(~ type) + theme_bw() + coord_fixed() +
  fill_scale(name = "dBZ")


## ------------------------------------------------------------------------
es_sample(radar_val0935$z, dat = as.matrix(Sims_IDE))
es_sample(radar_val0935$z, dat = as.matrix(Sims_FRK))

## ------------------------------------------------------------------------
distances <- radar_val0935 %>%
           coordinates() %>%
           dist() %>%
           as.matrix()
weights <- 0.5^distances

## ------------------------------------------------------------------------
vs_sample(radar_val0935$z, dat = as.matrix(Sims_IDE),
          w = weights, p = 1)

vs_sample(radar_val0935$z, dat = as.matrix(Sims_FRK),
          w = weights, p = 1)

## ------------------------------------------------------------------------
loglikFRK <- FRK:::loglik(S)
loglikIDE <- -fit_results_radar2$IDEmodel$negloglik()

## ------------------------------------------------------------------------
pIDE <- 7
pFRK <- 9

## ------------------------------------------------------------------------
m <- length(radar_obs)

## ------------------------------------------------------------------------
## Initialize table
Criteria <- data.frame(AIC = c(0, 0), BIC = c(0, 0),
                         row.names = c("FRK", "IDE"))

## Compute criteria
Criteria["FRK", "AIC"] <- -2*loglikFRK + 2*pFRK
Criteria["IDE", "AIC"] <- -2*loglikIDE + 2*pIDE
Criteria["FRK", "BIC"] <- -2*loglikFRK + pFRK*log(m)
Criteria["IDE", "BIC"] <- -2*loglikIDE + pIDE*log(m)
Criteria
