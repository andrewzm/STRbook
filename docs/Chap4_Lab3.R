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

library("dplyr")
library("ggplot2")
library("gstat")
library("RColorBrewer")
library("sp")
library("spacetime")
library("SpatioTemporal")
library("STRbook")
library("tidyr")

## ------------------------------------------------------------------------
data("NOAA_df_1990", package = "STRbook")   # load NOAA data
NOAA_sub <- filter(NOAA_df_1990,       # filter data to only
                   year == 1993 &      # contain July 1993
                   month == 7 &
                   proc == "Tmax")     # and just max. temp.

NOAA_sub_for_STdata <- NOAA_sub %>%
                       transmute(ID = as.character(id),
                                 obs = z,
                                 date = date)

## ------------------------------------------------------------------------
covars <-  dplyr::select(NOAA_sub, id, lat, lon) %>%
           unique() %>%
           dplyr::rename(ID = id)      # createSTdata expects "ID"

## ------------------------------------------------------------------------
STdata <- createSTdata(NOAA_sub_for_STdata, covars = covars)

## ------------------------------------------------------------------------
plot(STdata, "acf", ID = "3812")

## ------------------------------------------------------------------------
STdata <- updateTrend(STdata, n.basis = 2)

## ------------------------------------------------------------------------
## See how many stations are no longer lag-1 correlated
STdata1 <- createSTdata(NOAA_sub_for_STdata, covars = covars)
STdata2 <- updateTrend(STdata, n.basis = 2)
DFcheck <- left_join(NOAA_sub_for_STdata, STdata2$trend)
Insig1 <- Insig2 <- NULL
for(i in unique(STdata1$obs$ID)) {
  ACF1 <- acf(residuals(lm(obs ~ 1, data = filter(DFcheck, ID == i))), plot = 0)
  ACF2 <- acf(residuals(lm(obs ~ V1 + V2, data = filter(DFcheck, ID == i))), plot = 0)
  Insig1[[i]] <- ACF1$acf[2] < qnorm((1 + 0.95)/2)/sqrt(ACF1$n.used) # Taken from plot.act
  Insig2[[i]] <- ACF2$acf[2] < qnorm((1 + 0.95)/2)/sqrt(ACF2$n.used)
}
Sig1 <- 100 - round(mean(unlist(Insig1)),2)*100
Sig2 <- 100 - round(mean(unlist(Insig2)),2)*100

## ------------------------------------------------------------------------
plot(STdata, "acf", ID = "3812")

## ------------------------------------------------------------------------
STdata <- createSTdata(NOAA_sub_for_STdata, covars = covars)
plot(STdata, "acf", ID = "3812")

## ------------------------------------------------------------------------
STdata <- updateTrend(STdata, n.basis = 2)
plot(STdata, "acf", ID = "3812")

## ------------------------------------------------------------------------
data_to_plot <- gather(STdata$trend,trend, z, -date)
data_to_plot <- rbind(data=data.frame(date = as.Date(c("1993-07-01", "1993-07-31")),
                                      trend = "V0",
                                      z = c(1, 1)),
                      data_to_plot)
LineTrends <- ggplot(data_to_plot) +
    geom_line(aes(x=date,y=z,group=trend,linetype=trend)) +
    scale_linetype_manual(labels=c(expression(varphi[1](t)),
                                   expression(varphi[2](t)),
                                   expression(varphi[3](t))),
                          values=c(1,3,4),name="basis") +
    ylab("") + theme_bw()


## ------------------------------------------------------------------------
beta.lm <- estimateBetaFields(STdata)

## ------------------------------------------------------------------------
head(row.names(beta.lm$beta))
head(covars$ID)

## ------------------------------------------------------------------------
beta.lm$beta <- data.frame(beta.lm$beta)
beta.lm$beta.sd <- data.frame(beta.lm$beta.sd)
beta.lm$beta$ID <- as.integer(row.names(beta.lm$beta))
BETA <- cbind(beta.lm$beta, beta.lm$beta.sd)
colnames(BETA) <- c("alpha1", "alpha2", "alpha3", "ID",
                    "alpha1_CI", "alpha2_CI", "alpha3_CI")
BETA <- left_join(BETA, covars, by = "ID")


## ------------------------------------------------------------------------
g1 <- ggplot(BETA) + geom_point(aes(x=lat,y=alpha1)) +
    geom_errorbar(aes(x = lat,
                      ymin = alpha1 - 1.96*alpha1_CI,
                      ymax = alpha1 + 1.96*alpha1_CI)) +
    ylab(expression(alpha[1](s))) +
    xlab("lat (deg)") +
    theme_bw()

g2 <- ggplot(BETA) + geom_point(aes(x=lat,y=alpha2)) +
    geom_errorbar(aes(x = lat,
                      ymin = alpha2 - 1.96*alpha2_CI,
                      ymax = alpha2 + 1.96*alpha2_CI)) +
    ylab(expression(alpha[2](s))) +
    xlab("lat (deg)") +
    theme_bw()


g3 <- ggplot(BETA) + geom_point(aes(x=lat,y=alpha3)) +
    geom_errorbar(aes(x = lat,
                      ymin = alpha3 - 1.96*alpha3_CI,
                      ymax = alpha3 + 1.96*alpha3_CI)) +
    ylab(expression(alpha[3](s))) +
    xlab("lat (deg)") +
    theme_bw()

## ------------------------------------------------------------------------
cov.beta <- list(covf = "exp", nugget = FALSE)

## ------------------------------------------------------------------------
cov.nu <- list(covf = "exp",
               nugget = ~1,
               random.effect = FALSE) # No random mean
                                      # for each nu

## ------------------------------------------------------------------------
locations <- list(coords = c("lon", "lat"))
LUR <- list(~lat, ~lat, ~1)  # lat trend for phi1 and phi2 only
STmodel <- createSTmodel(STdata,              # data
                         LUR = LUR,           # spatial covariates
                         cov.beta = cov.beta, # cov. of alphas
                         cov.nu = cov.nu,     # cov. of nu
                         locations = locations) # coord. names

## ------------------------------------------------------------------------
parnames <- loglikeSTnames(STmodel, all = FALSE)
print(parnames)

## ------------------------------------------------------------------------
x.init <- matrix(3, 9, 1)
rownames(x.init) <- loglikeSTnames(STmodel, all = FALSE)
SpatioTemporalfit1 <- estimate(STmodel, x.init)

## ------------------------------------------------------------------------
x.final <- coef(SpatioTemporalfit1, pars = "cov")$par

## ------------------------------------------------------------------------
## Define space-time grid
spat_pred_grid <- expand.grid(lon = seq(-100, -80, length = 20),
                      lat = seq(32, 46, length = 20))
spat_pred_grid$id <- 1:nrow(spat_pred_grid)
temp_pred_grid <- as.Date("1993-07-01") + seq(3, 28, length = 6)

## Initialize data matrix
obs_pred_wide <- matrix(0, nrow = 6, ncol = 400)

## Set row names and column names
rownames(obs_pred_wide) <- as.character(temp_pred_grid)
colnames(obs_pred_wide) <- spat_pred_grid$id

covars_pred <- spat_pred_grid                     # covariates
STdata_pred <- createSTdata(obs = obs_pred_wide,  # ST object
                            covars = covars_pred)

## ------------------------------------------------------------------------
E <- predict(STmodel, x.final, STdata = STdata_pred)

## ------------------------------------------------------------------------
all_dates <- NOAA_sub$date %>% unique()     # dates
lookup <- data.frame(date = all_dates,      # covariate (linear)
                     V1 = scale(as.numeric(all_dates)))

## ------------------------------------------------------------------------
## Function that returns the covariates in a data frame
## at the required dates
fnc <- function(dates) {
  left_join(data.frame(date = dates),
            lookup, by = "date") %>%
  select(-date)
}

## ------------------------------------------------------------------------
STdata <- updateTrend(STdata, fnc = fnc)
