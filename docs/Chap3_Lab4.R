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

library("ape")
library("dplyr")
library("FRK")
library("ggplot2")
library("gstat")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

## ------------------------------------------------------------------------
data("MOcarolinawren_long", package = "STRbook")
MOcarolinawren_long <- MOcarolinawren_long %>%
                       filter(!is.na(cnt))

## ------------------------------------------------------------------------
G <- auto_basis(data = MOcarolinawren_long[,c("lon","lat")] %>%
                       SpatialPoints(),           # To sp obj
                nres = 1,                         # One resolution
                type = "Gaussian")                # Gaussian BFs

S <- eval_basis(basis = G,                       # basis functions
                s = MOcarolinawren_long[,c("lon","lat")] %>%
                     as.matrix()) %>%            # conv. to matrix
     as.matrix()                                 # conv. to matrix
colnames(S) <- paste0("B", 1:ncol(S)) # assign column names

## ------------------------------------------------------------------------
Wren_df <- cbind(MOcarolinawren_long,S) %>%
  select(-loc.ID, -t)
Wren_df[1:3, 1:5]

## ------------------------------------------------------------------------
Wren_GLM <- glm(cnt ~ (lon + lat + year)^2 + ., # formula
                family = poisson("log"),     # Poisson + log link
                data = Wren_df)              # data set

## ------------------------------------------------------------------------
options(digits = 3)

## ------------------------------------------------------------------------
Wren_GLM$deviance / Wren_GLM$df.residual

## ------------------------------------------------------------------------
## Wren_GLM_QP <- glm(cnt ~ (lon + lat + year)^2 + ., # formula
##                family = quasipoisson("log"),        # Poisson + log link
##                data = Wren_df)                 # data set

## ------------------------------------------------------------------------
Wren_GLM$df.residual

## ------------------------------------------------------------------------
Wren_GLM$deviance

## ------------------------------------------------------------------------
1 - pchisq(q = Wren_GLM$deviance, df = Wren_GLM$df.residual)

## ------------------------------------------------------------------------
pred_grid <- expand.grid(lon = seq(
                             min(MOcarolinawren_long$lon) - 0.2,
                             max(MOcarolinawren_long$lon) + 0.2,
                             length.out = 80),
                         lat = seq(
                             min(MOcarolinawren_long$lat) - 0.2,
                             max(MOcarolinawren_long$lat) + 0.2,
                             length.out = 80),
                         year = 1994:2014)

## ------------------------------------------------------------------------
S_pred <- eval_basis(basis = G,                    # basis functs
                s = pred_grid[,c("lon","lat")] %>% # pred locs
                     as.matrix()) %>%            # conv. to matrix
     as.matrix()                                 # as matrix
colnames(S_pred) <- paste0("B", 1:ncol(S_pred))  # assign  names
pred_grid <- cbind(pred_grid,S_pred)             # attach to grid

## ------------------------------------------------------------------------
wren_preds <- predict(Wren_GLM,
                      newdata = pred_grid,
                      type = "link",
                      se.fit = TRUE)

## ------------------------------------------------------------------------
pred_grid <- pred_grid %>%
             mutate(log_cnt = wren_preds$fit,
                    se = wren_preds$se.fit)


g1 <- ggplot() + geom_raster(data=pred_grid,
                             aes(lon, lat, fill = pmax(pmin(log_cnt,4),0))) +
  facet_wrap(~year,nrow=3,ncol=7) +
    geom_point(data = filter(MOcarolinawren_long, !is.na(cnt)),
               aes(lon, lat),colour="black", size=3) +
    geom_point(data=filter(MOcarolinawren_long,!is.na(cnt)),aes(lon,lat,colour=log(cnt)),size=2) +
    scale_colour_distiller(palette="Spectral",limits=c(0,4)) +
    scale_fill_distiller(palette="Spectral",limits=c(0,4),name=expression(log(Y[t]))) + theme_bw()

g2 <- ggplot() + geom_raster(data=pred_grid,aes(lon,lat,fill=se)) +
  facet_wrap(~year,nrow=3,ncol=7) +
    scale_fill_distiller(palette="BrBG",limits=c(0,1),name=expression(s.e.)) + theme_bw()

## ------------------------------------------------------------------------
Wren_df$residuals <- residuals(Wren_GLM)

## ------------------------------------------------------------------------
g2 <- ggplot(Wren_df) +
    geom_point(aes(lon, lat, colour = residuals)) +
    col_scale(name = "residuals") +
    facet_wrap(~year, nrow = 3) + theme_bw()

## ------------------------------------------------------------------------
P <- list()                                 # init list
years <- 1994:2014
for(i in seq_along(years)) {                # for each day
  Wren_year <- filter(Wren_df,
                     year == years[i])      # filter by year
  obs_dists <- Wren_year %>%                # take the data
    select(lon,lat) %>%                     # extract coords.
    dist() %>%                              # comp. dists.
    as.matrix()                             # conv. to matrix
  obs_dists.inv <- 1/obs_dists              # weight matrix
  diag(obs_dists.inv) <- 0                  # 0 on diag
  P[[i]] <- Moran.I(Wren_year$residuals,    # run Moran's I
                    obs_dists.inv) %>%
            do.call("cbind", .)             # conv. to df
}
do.call("rbind",P) %>% summary(digits = 2)

## ------------------------------------------------------------------------
Wren_STIDF <- STIDF(sp = SpatialPoints(
                            Wren_df[,c("lon","lat")],
                            proj4string = CRS("+proj=longlat")),
                    time = as.Date(Wren_df[, "year"] %>%
                                       as.character(),
                                   format = "%Y"),
                    data = Wren_df)

## ------------------------------------------------------------------------
tlags <- seq(0.01, 52.1429*6 + 0.01, by = 52.1429)
vv <- variogram(object = residuals ~ 1, # fixed effect component
                data = Wren_STIDF,      # data set
                tlags = tlags,          # temp. bins
                width = 25,             # spatial bin (25 km)
                cutoff = 150,           # use pts < 150 km apart
                tunit = "weeks")        # time unit


