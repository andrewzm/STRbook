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

library("leaps")
library("lmtest")
library("nlme")
library("ape")
library("broom")
library("FRK")
library("purrr")
library("lattice")
library("ggplot2")
library("RColorBrewer")
library("dplyr")
library("gstat")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

## ------------------------------------------------------------------------
data("NOAA_df_1990", package = "STRbook")
Tmax <- filter(NOAA_df_1990,       # subset the data
              proc == "Tmax" &     # only max temperature
              month == 7 &         # July
              year == 1993)        # year of 1993

## ------------------------------------------------------------------------
G <- auto_basis(data = Tmax[,c("lon","lat")] %>%  # Take Tmax
                       SpatialPoints(),           # To sp obj
                nres = 1,                         # One resolution
                type = "Gaussian")                # Gaussian BFs

## ------------------------------------------------------------------------
S <- eval_basis(basis = G,                     # basis functions
                s = Tmax[,c("lon","lat")] %>%  # spat locations
                     as.matrix()) %>%          # conv. to matrix
     as.matrix()                               # results as matrix
colnames(S) <- paste0("B", 1:ncol(S)) # assign column names

## ------------------------------------------------------------------------
Tmax2 <- cbind(Tmax, S) %>%             # append S to Tmax
         select(-year, -month, -proc,   # and remove vars we
                -julian, -date)         # will not be using in
                                        # the model

## ------------------------------------------------------------------------
Tmax_no_14 <- filter(Tmax2, !(day == 14))  # remove day 14

## ------------------------------------------------------------------------
Tmax_July_lm <- lm(z ~ (lon + lat + day)^2 + .,     # model
                   data = select(Tmax_no_14, -id))  # omit id

## ------------------------------------------------------------------------
options(width = 55)

## ------------------------------------------------------------------------
Tmax_July_lm %>% summary()

## ------------------------------------------------------------------------
Tmax_July_gls <- gls(z ~ (lon + lat + day)^2 + .,
                     data = select(Tmax_no_14, -id),
                     correlation = corGaus(value = 0.5,
                                        form = ~ lon + lat + day,
                                        fixed = TRUE))

## ------------------------------------------------------------------------
Tmax_July_lm4 <- list()   # initialize
for(i in 0:4) {           # for four steps (after intercept model)
   ## Carry out stepwise forward selection for i steps
   Tmax_July_lm4[[i+1]] <- step(lm(z ~ 1,
                               data = select(Tmax_no_14, -id)),
                               scope = z ~(lon + lat + day)^2 + .,
                               direction = 'forward',
                               steps = i)
}

## ------------------------------------------------------------------------
regfit.full = regsubsets(z ~ 1 + (lon + lat + day)^2 + .,  # model
                         data = select(Tmax_no_14, -id),
                         method = "forward",
                         nvmax = 4)                  # 4 steps

## ------------------------------------------------------------------------
regfit.summary <- summary(regfit.full)

## ------------------------------------------------------------------------
set.seed(1) # Fix seed for reproducibility
Tmax_no_14_2 <- Tmax_no_14 %>%
                mutate(B13 = B5 + 0.01*rnorm(nrow(Tmax_no_14)))

## ------------------------------------------------------------------------
Tmax_July_lm3 <- lm(z ~ (lon + lat + day)^2 + .,
                   data = Tmax_no_14_2 %>%
                          select(-id))

## ------------------------------------------------------------------------
summary(Tmax_July_lm3)

## ------------------------------------------------------------------------
vcov(Tmax_July_lm3)[c("B5", "B13"),c("B5", "B13")] %>%
    cov2cor()

## ------------------------------------------------------------------------
Tmax_no_14$residuals <- residuals(Tmax_July_lm)

## ------------------------------------------------------------------------
g <- ggplot(filter(Tmax_no_14, day %in% 24:31)) +
  geom_point(aes(lon, lat, colour = residuals)) +
  facet_wrap(~ day, ncol=4) +
  col_scale(name = "degF") +
  geom_point(data = filter(Tmax_no_14,day %in% 24:31 &
                                        id %in% c(3810, 3889)),
               aes(lon, lat), colour = "black",
               pch = 2, size = 2.5) +
  theme_bw()

## ------------------------------------------------------------------------
print(g)

## ------------------------------------------------------------------------
P <- list()                                 # init list
days <- c(1:13, 15:31)                      # set of days
for(i in seq_along(days)) {                 # for each day
  Tmax_day <- filter(Tmax_no_14,
                     day == days[i])        # filter by day
  station.dists <- Tmax_day %>%             # take the data
    select(lon, lat) %>%                    # extract coords.
    dist() %>%                              # comp. dists.
    as.matrix()                             # conv. to matrix
  station.dists.inv <- 1/station.dists      # weight matrix
  diag(station.dists.inv) <- 0              # 0 on diag
  P[[i]] <- Moran.I(Tmax_day$residuals,     # run Moran's I
                    station.dists.inv) %>%
            do.call("cbind", .)             # conv. to df
}

## ------------------------------------------------------------------------
do.call("rbind", P) %>% head()

## ------------------------------------------------------------------------
station.dists <- Tmax_no_14 %>%  # take the data
  select(lon, lat, day) %>%      # extract coordinates
  dist() %>%                     # compute distances
  as.matrix()                    # convert to matrix

## ------------------------------------------------------------------------
station.dists.inv <- 1/station.dists
diag(station.dists.inv) <- 0
Moran.I(Tmax_no_14$residuals, station.dists.inv)$p.value

## ------------------------------------------------------------------------
TS1 <- filter(Tmax_no_14, id == 3810)$residuals
TS2 <- filter(Tmax_no_14, id == 3889)$residuals

## ------------------------------------------------------------------------
par(mar=c(4, 4, 1, 1))
plot(TS1,                           # Station 3810 residuals
     xlab = "day of July 1993",
     ylab = "residuals (degF)",
     type = 'o', ylim = c(-8, 7))
lines(TS2,                          # Station 3889 residuals
      xlab = "day of July 1993",
      ylab = "residuals (degF)",
      type = 'o', col = "red")


## ------------------------------------------------------------------------
nested_Tmax_no_14 <- group_by(Tmax_no_14, lon, lat) %>% nest()
head(nested_Tmax_no_14, 3)

## ------------------------------------------------------------------------
dwtest_one_station <- function(data)
                        dwtest(residuals ~ 1, data = data)

## ------------------------------------------------------------------------
dwtest_one_station(nested_Tmax_no_14$data[[1]])

## ------------------------------------------------------------------------
map(nested_Tmax_no_14$data, dwtest_one_station) %>% head()

## ------------------------------------------------------------------------
dwtest_one_station_tidy <- nested_Tmax_no_14$data[[1]] %>%
                           dwtest_one_station() %>%
                           tidy()

## ------------------------------------------------------------------------
dwtest_one_station_tidy[, 1:3]

## ------------------------------------------------------------------------
Tmax_DW_no_14 <- nested_Tmax_no_14 %>%
    mutate(dwtest = map(data, dwtest_one_station)) %>%
    mutate(test_df = map(dwtest, tidy)) %>%
    unnest(test_df)

## ------------------------------------------------------------------------
options(digits = 3)

## ------------------------------------------------------------------------
Tmax_DW_no_14 %>% select(-method, -alternative) %>% head(3)

## ------------------------------------------------------------------------
mean(Tmax_DW_no_14$p.value < 0.05/nrow(Tmax_DW_no_14)) * 100

## ------------------------------------------------------------------------
options(digits=1)

## ------------------------------------------------------------------------
data("STObj3", package = "STRbook")
STObj4 <- STObj3[, "1993-07-01::1993-07-31"]

## ------------------------------------------------------------------------
STObj4@data <- left_join(STObj4@data, Tmax_no_14)

## ------------------------------------------------------------------------
vv <- variogram(object = residuals ~ 1, # fixed effect component
                data = STObj4,     # July data
                width = 80,        # spatial bin (80 km)
                cutoff = 1000,     # consider pts < 1000 km apart
                tlags = 0.01:6.01) # 0 days to 6 days


## ------------------------------------------------------------------------
pred_grid <- expand.grid(lon = seq(-100, -80, length = 20),
                         lat = seq(32, 46, length = 20),
                         day = seq(4, 29, length = 6))

## ------------------------------------------------------------------------
Spred <- eval_basis(basis = G,                      # basis functs
                s = pred_grid[,c("lon","lat")] %>%  # pred locs
                     as.matrix()) %>%         # conv. to matrix
     as.matrix()                              # results as matrix
colnames(Spred) <- paste0("B", 1:ncol(Spred)) # assign col names
pred_grid <- cbind(pred_grid, Spred)          # attach to grid

## ------------------------------------------------------------------------
linreg_pred <- predict(Tmax_July_lm,
                       newdata = pred_grid,
                       interval = "prediction")

## ------------------------------------------------------------------------
## Assign prediction and prediction s.e. to the prediction grid
pred_grid$z_pred <- linreg_pred[,1]
pred_grid$z_err <- (linreg_pred[,3] - linreg_pred[,2]) / (2*1.96)
