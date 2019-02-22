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

library("plyr")
library("dplyr")
library("IDE")
library("FRK")
library("ggplot2")
library("sp")
library("spacetime")
library("STRbook")

## ------------------------------------------------------------------------
SIM1 <- simIDE(T = 10, nobs = 100, k_spat_invariant = 1)

## ------------------------------------------------------------------------
print(SIM1$g_truth)
print(SIM1$g_obs)

## ------------------------------------------------------------------------
IDEmodel <- IDE(f = z ~ s1 + s2,
                data = SIM1$z_STIDF,
                dt = as.difftime(1, units = "days"),
                grid_size = 41)

## ------------------------------------------------------------------------
## fit_results_sim1 <- fit.IDE(IDEmodel,
##                            parallelType = 1)

## ------------------------------------------------------------------------
data("IDE_Sim1_results", package = "STRbook")

## ------------------------------------------------------------------------
show_kernel(fit_results_sim1$IDEmodel)

## ------------------------------------------------------------------------
fit_results_sim1$IDEmodel$get("k") %>% unlist()

## ------------------------------------------------------------------------
coef(fit_results_sim1$IDEmodel)

## ------------------------------------------------------------------------
abs_ev <- eigen(fit_results_sim1$IDEmodel$get("M"))$values %>%
          abs()
summary(abs_ev)

## ------------------------------------------------------------------------
ST_grid_df <- predict(fit_results_sim1$IDEmodel)

## ------------------------------------------------------------------------
gpred <- ggplot(ST_grid_df) +       # Plot the predictions
  geom_tile(aes(s1, s2, fill=Ypred)) +
  facet_wrap(~t) +
  fill_scale(name = "Ypred", limits = c(-0.1, 1.4)) +
  coord_fixed(xlim=c(0, 1), ylim = c(0, 1))

gpredse <- ggplot(ST_grid_df) +     # Plot the prediction s.es
  geom_tile(aes(s1, s2, fill = Ypredse)) +
  facet_wrap(~t) +
  fill_scale(name = "Ypredse") +
  coord_fixed(xlim=c(0, 1), ylim = c(0, 1))


## ------------------------------------------------------------------------
SIM2 <- simIDE(T = 15, nobs = 1000, k_spat_invariant = 0)

## ----results = 'hide', fig.keep = 'none'---------------------------------
print(SIM2$g_truth)
print(SIM2$g_obs)

## ------------------------------------------------------------------------
show_kernel(SIM2$IDEmodel, scale = 0.2)

## ------------------------------------------------------------------------
mbasis_1 <- auto_basis(manifold = plane(),   # fns on the plane
                       data = SIM2$z_STIDF,  # data
                       nres = 1,             # 1 resolution
                       type = 'bisquare')    # type of functions

## ------------------------------------------------------------------------
kernel_basis <- list(thetam1 = constant_basis(),
                     thetam2 = constant_basis(),
                     thetam3 = mbasis_1,
                     thetam4 = mbasis_1)

## ------------------------------------------------------------------------
IDEmodel <- IDE(f = z ~ s1 + s2 + 1,
                data = SIM2$z_STIDF,
                dt = as.difftime(1, units = "days"),
                grid_size = 41,
                kernel_basis = kernel_basis)

## ------------------------------------------------------------------------
## fit_results_sim2 <- fit.IDE(IDEmodel,
##                            parallelType = 1,
##                            itermax = 400)

## ------------------------------------------------------------------------
data("IDE_Sim2_results", package = "STRbook")

## ------------------------------------------------------------------------
show_kernel(fit_results_sim2$IDEmodel)

## ------------------------------------------------------------------------
K1 <- show_kernel(fit_results_sim2$IDEmodel, scale = 0.2) + ggplot2::coord_fixed()
K2 <- show_kernel(SIM2$IDEmodel, scale = 0.2) + ggplot2::coord_fixed()
g <- grid.arrange(K1, K2, ncol = 2)

## ------------------------------------------------------------------------
data("radar_STIDF", package = "STRbook")

## ------------------------------------------------------------------------
IDEmodel <- IDE(f = z ~ 1,
                data = radar_STIDF,
                dt = as.difftime(10, units = "mins"),
                grid_size = 41,
                forecast = 2,
                hindcast = 2)

## ------------------------------------------------------------------------
## fit_results_radar <- fit.IDE(IDEmodel,
##                              parallelType = 1)

## ------------------------------------------------------------------------
data("IDE_Radar_results", package = "STRbook")

## ------------------------------------------------------------------------
show_kernel(fit_results_radar$IDEmodel)

## ------------------------------------------------------------------------
options(digits = 2)

## ------------------------------------------------------------------------
shift_pars <- (fit_results_radar$IDEmodel$get("k") %>%
                   unlist())[3:4]
print(shift_pars)

## ------------------------------------------------------------------------
abs_ev <- eigen(fit_results_radar$IDEmodel$get("M"))$values %>%
          abs()
summary(abs_ev)

## ------------------------------------------------------------------------
ST_grid_df <- predict(fit_results_radar$IDEmodel)

## ------------------------------------------------------------------------
radar_df$time <- format(radar_df$t, "%H:%M")
ST_grid_df$time <- format(ST_grid_df$t, "%H:%M")

## ------------------------------------------------------------------------
## Add time records with missing data
radar_df <- rbind.fill(radar_df,
                       data.frame(time = c("08:05", "08:15",
                                           "10:25", "10:35")))

## Plot of data, with color scale capped to (-20, 60)
gobs <- ggplot(radar_df) +
  geom_tile(aes(s1, s2, fill = pmin(pmax(z, -20), 60))) +
  fill_scale(limits = c(-20, 60), name = "Z") +
  facet_wrap(~time) + coord_fixed() + theme_bw()

## Plot of predictions with color scale forced to (-20, 60)
gpred <- ggplot(ST_grid_df) +
  geom_tile(aes(s1, s2, fill = Ypred)) +
  facet_wrap(~time) + coord_fixed() + theme_bw() +
  fill_scale(limits = c(-20, 60), name = "Ypred")


