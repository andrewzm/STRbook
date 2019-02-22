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
library("FRK")
library("ggplot2")
library("gstat")
library("RColorBrewer")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")
library("grid")
library("gridExtra")
library("sp")
library("FRK")
sp::set_Polypath(FALSE)

## ------------------------------------------------------------------------
data("STObj3", package = "STRbook")          # load STObj3
STObj4 <- STObj3[, "1993-07-01::1993-07-31"] # subset time
STObj5 <- as(STObj4[, -14], "STIDF")         # omit t = 14
STObj5 <- subset(STObj5, !is.na(STObj5$z))   # remove NAs

## ------------------------------------------------------------------------
BAUs <- auto_BAUs(manifold = STplane(),   # ST field on the plane
                  type = "grid",          # gridded (not "hex")
                  data = STObj5,          # data
                  cellsize = c(1, 0.75, 1), # BAU cell size
                  convex = -0.12,           # hull extension
                  tunit = "days")           # time unit is "days"

## ------------------------------------------------------------------------
plot(as(BAUs[, 1], "SpatialPixels"))    # plot pixel BAUs
plot(SpatialPoints(STObj5),
    add = TRUE, col = "red")           # plot data points

## ------------------------------------------------------------------------
BAUs_hex <- auto_BAUs(manifold = STplane(), # model on the plane
                  type = "hex",             # hex (not "grid")
                  data = STObj5,            # data
                  cellsize = c(1, 0.75, 1), # BAU cell size
                  nonconvex_hull = FALSE,   # convex hull
                  tunit = "days")           # time unit is "days"

## ------------------------------------------------------------------------
plot(as(BAUs_hex[, 1], "SpatialPolygons"))

## ------------------------------------------------------------------------
G_spatial <- auto_basis(manifold = plane(),      # fns on plane
                        data = as(STObj5, "Spatial"), # project
                        nres = 2,                     # 2 res.
                        type = "bisquare",            # bisquare.
                        regular = 0)                  # irregular

## ------------------------------------------------------------------------
t_grid <- matrix(seq(1, 31, length = 20))

## ------------------------------------------------------------------------
G_temporal <- local_basis(manifold = real_line(),  # fns on R1
                          type = "bisquare",       # bisquare
                          loc = t_grid,            # centroids
                          scale = rep(2, 20))      # aperture par.

## ------------------------------------------------------------------------
G <- TensorP(G_spatial, G_temporal)      # take the tensor product

## ------------------------------------------------------------------------
g1 <- show_basis(G_spatial) + xlab("lon (deg)") + ylab("lat (deg)") + coord_fixed()
g2 <- show_basis(G_temporal) + xlab("t (days)") +ylab(expression(varphi(t)))

## ------------------------------------------------------------------------
BAUs$fs = 1

## ------------------------------------------------------------------------
STObj5$std <- sqrt(0.049)

## ------------------------------------------------------------------------
f <- z ~ lat + 1

## ------------------------------------------------------------------------
S <- FRK(f = f,               # formula
         data = list(STObj5), # (list of) data
         basis = G,           # basis functions
         BAUs = BAUs,         # BAUs
         n_EM = 3,            # max. no. of EM iterations
         tol = 0.01)          # tol. on change in log-likelihood

## ------------------------------------------------------------------------
grid_BAUs <- predict(S)

## ------------------------------------------------------------------------
grid_BAUs@sp  <- SpatialPoints(grid_BAUs@sp)   # convert to Spatial points object
gridded(grid_BAUs@sp) <- TRUE                          # and assert that it is gridded
library(RColorBrewer)
 colour_regions <- (brewer.pal(11, "Spectral") %>%     # construct spectral brewer palette
                        colorRampPalette())(16) %>%    # at 16 levels
                         rev()                         # reverse direction

 grid_BAUs$mu_trunc <- pmin(pmax(grid_BAUs$mu,68),105)
 stplot(grid_BAUs[,c(4,9,14,19,24,29),"mu_trunc"],             # plot the FRK predictor
      main="Predictions (degrees Fahrenheit)",          # title
     layout=c(3,2),                                     # trellis layout
     col.regions=colour_regions,                        # color scale
     xlim=c(-100,-80),ylim=c(32,46),                    # axes limits
     aspect=1)                                          # fixed aspect ratio
 
 ## ------------------------------------------------------------------------
grid_BAUs$se <- pmax(pmin(sqrt(grid_BAUs$var),4.4),1.8)
stplot(grid_BAUs[,c(4,9,14,19,24,29),"se"],             # plot the FRK predictor
      main="Prediction std. errors (degrees Fahrenheit)",          # title
     layout=c(3,2),                                     # trellis layout
     col.regions=colour_regions,                        # color scale
     xlim=c(-100,-80),ylim=c(32,46),                    # axes limits
     aspect=1)                                          # fixed aspect ratio

