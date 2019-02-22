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

library("INLA")
library("dplyr")
library("tidyr")
library("ggplot2")
library("STRbook")
data("MOcarolinawren_long", package = "STRbook")

## ------------------------------------------------------------------------
coords <- unique(MOcarolinawren_long[c("loc.ID", "lon", "lat")])
boundary <- inla.nonconvex.hull(as.matrix(coords[, 2:3]))

## ------------------------------------------------------------------------
MOmesh <- inla.mesh.2d(boundary = boundary,
                       max.edge = c(0.8, 1.2), # max. edge length
                       cutoff = 0.1)           # min. edge length

plot(MOmesh,asp=1,main="");
lines(coords[c("lon","lat")],col="red",type="p")

## ------------------------------------------------------------------------
plot(MOmesh, asp = 1, main = "")
lines(coords[c("lon", "lat")], col = "red", type = "p")

## ------------------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = MOmesh,
                            alpha = 2,
                            prior.range = c(1, 0.01),
                            prior.sigma = c(4, 0.01))

## ------------------------------------------------------------------------
n_years <- length(unique(MOcarolinawren_long$t))
n_spatial <- MOmesh$n
s_index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = n_spatial,
                                n.group = n_years)


## ------------------------------------------------------------------------
coords.allyear <- MOcarolinawren_long[c("lon", "lat")] %>%
                  as.matrix()
PHI <- inla.spde.make.A(mesh = MOmesh,
                        loc = coords.allyear,
                        group = MOcarolinawren_long$t,
                        n.group = n_years)

## ------------------------------------------------------------------------
dim(PHI)

## ------------------------------------------------------------------------
nrow(MOcarolinawren_long)
length(s_index$spatial.field)

## ------------------------------------------------------------------------
## First stack: Estimation
n_data <- nrow(MOcarolinawren_long)
stack_est <- inla.stack(
                 data = list(cnt = MOcarolinawren_long$cnt),
                 A = list(PHI, 1),
                 effects = list(s_index,
                                list(Intercept = rep(1, n_data))),
                 tag = "est")

## ------------------------------------------------------------------------
df_pred <- data.frame(lon = rep(MOmesh$loc[,1], n_years),
                      lat = rep(MOmesh$loc[,2], n_years),
                      t = rep(1:n_years, each = MOmesh$n))
n_pred <- nrow(df_pred)
PHI_pred <- Diagonal(n = n_pred)

## ------------------------------------------------------------------------
## Second stack: Prediction
stack_pred <- inla.stack(
                 data = list(cnt = NA),
                 A = list(PHI_pred, 1),
                 effects = list(s_index,
                                list(Intercept = rep(1, n_pred))),
                 tag = "pred")

## ------------------------------------------------------------------------
stack <- inla.stack(stack_est, stack_pred)

## ------------------------------------------------------------------------
dim(stack_est$A)

## ------------------------------------------------------------------------
## PC prior on rho
rho_hyper <- list(theta=list(prior = 'pccor1',
                             param = c(0, 0.9)))

## Formula
formula <- cnt ~ -1 + Intercept +
                 f(spatial.field,
                   model = spde,
                   group = spatial.field.group,
                   control.group = list(model = "ar1",
                                        hyper = rho_hyper))

## ------------------------------------------------------------------------
## output <- inla(formula,
##                data = inla.stack.data(stack, spde = spde),
##                family = "nbinomial",
##                control.predictor = list(A = inla.stack.A(stack),
##                                         compute = TRUE))
##
## output2 <- list(summary.fixed = output$summary.fixed,
##                 summary.hyperpar = output$summary.hyperpar,
##                 summary.fitted.values = output$summary.fitted.values,
##                 marginals.fixed = output$marginals.fixed,
##                 marginals.hyperpar = output$marginals.hyperpar,
##                 internal.marginals.hyperpar = output$internal.marginals.hyperpar,
##                 marginals.spde2.blc = output$marginals.spde2.blc,
##                 marginals.spde3.blc = output$marginals.spde3.blc,
##                 internal.marginals.hyperpar = output$internal.marginals.hyperpar)
## class(output2) <- "inla"
## output <- output2

## ------------------------------------------------------------------------
data("INLA_output", package = "STRbook")

## ------------------------------------------------------------------------
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,4,2,2))
par(mfrow=c(2,2))
output.field <- inla.spde2.result(inla = output,
                                  name = "spatial.field",
                                  spde = spde,
                                  do.transf = TRUE)
plot(output$marginals.fix[[1]],type ='l',xlab=expression(beta[0]),ylab=expression(p(beta[0]*"|"*Z)))
plot(output$marginals.hyperpar[[4]],type = 'l',xlab = expression(rho),ylab = expression(p(rho*"|"*Z)))
plot(output.field$marginals.variance.nominal[[1]],type = 'l',xlab = expression(sigma^2),ylab = expression(p(sigma^2*"|"*Z)))
plot(output.field$marginals.range.nominal[[1]],type = 'l',xlab = expression(l),ylab = expression(p(l*"|"*Z)))

## Inference results on beta, rho, variance, and range parameters
round(output$summary.fixed,3)
round(output$summary.hyperpar,3)

## ------------------------------------------------------------------------
output.field <- inla.spde2.result(inla = output,
                                  name = "spatial.field",
                                  spde = spde,
                                  do.transf = TRUE)
## plot p(beta0 | Z)
plot(output$marginals.fix$Intercept,
     type = 'l',
     xlab = expression(beta[0]),
     ylab = expression(p(beta[0]*"|"*Z)))

## plot p(rho | Z)
plot(output$marginals.hyperpar$`GroupRho for spatial.field`,
     type = 'l',
     xlab = expression(rho),
     ylab = expression(p(rho*"|"*Z)))

## plot p(sigma^2 | Z)
plot(output.field$marginals.variance.nominal[[1]],
     type = 'l',
     xlab = expression(sigma^2),
     ylab = expression(p(sigma^2*"|"*Z)))

## plot p(range | Z)
plot(output.field$marginals.range.nominal[[1]],
     type = 'l',
     xlab = expression(l),
     ylab = expression(p(l*"|"*Z)))

## ------------------------------------------------------------------------
index_pred <- inla.stack.index(stack, "pred")$data
lp_mean <- output$summary.fitted.values$mean[index_pred]
lp_sd <- output$summary.fitted.values$sd[index_pred]

## ------------------------------------------------------------------------
grid_locs <- expand.grid(
                 lon = seq(min(MOcarolinawren_long$lon) - 0.2,
                           max(MOcarolinawren_long$lon) + 0.2,
                           length.out = 80),
                 lat = seq(min(MOcarolinawren_long$lat) - 0.2,
                           max(MOcarolinawren_long$lat) + 0.2,
                           length.out = 80))

## ------------------------------------------------------------------------
proj.grid <- inla.mesh.projector(MOmesh,
                    xlim = c(min(MOcarolinawren_long$lon) - 0.2,
                             max(MOcarolinawren_long$lon) + 0.2),
                    ylim = c(min(MOcarolinawren_long$lat) - 0.2,
                             max(MOcarolinawren_long$lat) + 0.2),
                    dims = c(80, 80))

## ------------------------------------------------------------------------
pred <- sd <-  NULL
for(i in 1:n_years) {
    ii <- (i-1)*MOmesh$n + 1
    jj <- i*MOmesh$n
    pred[[i]] <- cbind(grid_locs,
                       z = c(inla.mesh.project(proj.grid,
                                             lp_mean[ii:jj])),
                       t = i)
    sd[[i]] <- cbind(grid_locs,
                     z = c(inla.mesh.project(proj.grid,
                                             lp_sd[ii:jj])),
                     t = i)
}

## ------------------------------------------------------------------------
pred <- do.call("rbind", pred) %>% filter(!is.na(z))
sd <- do.call("rbind", sd) %>% filter(!is.na(z))

## ------------------------------------------------------------------------
df_tot <- group_by(MOcarolinawren_long,lon,lat) %>%
          dplyr::summarise(tot_cnt = sum(cnt,na.rm=T))

g1 <- ggplot() + geom_raster(data=pred,aes(lon,lat,fill=pmin(z,5))) + facet_wrap(~t,nrow=3,ncol=7) +
    geom_point(data=filter(MOcarolinawren_long,!is.na(cnt)),aes(lon,lat),colour="black",size=3) +
    geom_point(data=filter(MOcarolinawren_long,!is.na(cnt)),aes(lon,lat,colour=log(cnt)),size=2) +
    scale_colour_distiller(palette="Spectral",limits=c(-1,5)) +
    scale_fill_distiller(palette="Spectral",limits=c(-1,5),name=expression(log(Y[t]))) + theme_bw()

g2 <- ggplot() + geom_raster(data=sd,aes(lon,lat,fill=z)) + facet_wrap(~t,nrow=3,ncol=7) +
    scale_fill_distiller(palette="BrBG",limits=c(0,2.5),name=expression(s.e.)) + theme_bw()
