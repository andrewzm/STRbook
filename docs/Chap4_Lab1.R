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

library("sp")
library("spacetime")
library("ggplot2")
library("dplyr")
library("gstat")
library("RColorBrewer")
library("STRbook")
library("tidyr")
library("grid")
library("gridExtra")

## ------------------------------------------------------------------------
data("STObj3", package = "STRbook")
STObj4 <- STObj3[, "1993-07-01::1993-07-31"]
vv <- variogram(object = z ~ 1 + lat, # fixed effect component
                data = STObj4,      # July data
                width = 80,         # spatial bin (80 km)
                cutoff = 1000,      # consider pts < 1000 km apart
                tlags = 0.01:6.01)  # 0 days to 6 days


## ------------------------------------------------------------------------
vignette("spatio-temporal-kriging")

## ------------------------------------------------------------------------
sepVgm <- vgmST(stModel = "separable",
                space = vgm(10, "Exp", 400, nugget = 0.1),
                time = vgm(10, "Exp", 1, nugget = 0.1),
                sill = 20)
sepVgm <- fit.StVariogram(vv, sepVgm)

## ------------------------------------------------------------------------
metricVgm <- vgmST(stModel = "metric",
                   joint = vgm(100, "Exp", 400, nugget = 0.1),
                   sill = 10,
                   stAni = 100)
metricVgm <- fit.StVariogram(vv, metricVgm)

## ------------------------------------------------------------------------
metricMSE <- attr(metricVgm, "optim")$value
sepMSE <- attr(sepVgm, "optim")$value

## ------------------------------------------------------------------------
plot(vv, list(sepVgm, metricVgm), main = "Semi-variance")


## ------------------------------------------------------------------------
spat_pred_grid <- expand.grid(
                      lon = seq(-100, -80, length = 20),
                      lat = seq(32, 46, length = 20)) %>%
            SpatialPoints(proj4string = CRS(proj4string(STObj3)))
gridded(spat_pred_grid) <- TRUE

## ------------------------------------------------------------------------
temp_pred_grid <- as.Date("1993-07-01") + seq(3, 28, length = 6)

## ------------------------------------------------------------------------
DE_pred <- STF(sp = spat_pred_grid,    # spatial part
               time = temp_pred_grid)  # temporal part

## ------------------------------------------------------------------------
STObj5 <- as(STObj4[, -14], "STIDF")         # convert to STIDF
STObj5 <- subset(STObj5, !is.na(STObj5$z))   # remove missing data

## ------------------------------------------------------------------------
pred_kriged <- krigeST(z ~ 1 + lat,         # latitude trend
                       data = STObj5,       # data set w/o 14 July
                       newdata = DE_pred,   # prediction grid
                       modelList = sepVgm,  # semivariogram
                       computeVar = TRUE)   # compute variances

## ------------------------------------------------------------------------
color_pal <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(16))

## ------------------------------------------------------------------------
stplot(pred_kriged,
      main = "Predictions (degrees Fahrenheit)",
      layout = c(3, 2),
      col.regions = color_pal)


## ------------------------------------------------------------------------
pred_kriged$se <- sqrt(pred_kriged$var1.var)
stplot(pred_kriged[, , "se"],
      main = "Prediction std. errors (degrees Fahrenheit)",
      layout = c(3, 2),
      col.regions = color_pal)



