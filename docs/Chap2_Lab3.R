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

library("CCA")
library("dplyr")
library("ggplot2")
library("gstat")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")
library("grid")
library("gridExtra")

## ------------------------------------------------------------------------
set.seed(1)

## -----------------------------------------------------------
data("NOAA_df_1990", package = "STRbook")
Tmax <- filter(NOAA_df_1990,     # subset the data
              proc == "Tmax" &   # only max temperature
              month %in% 5:9 &   # May to September
              year == 1993)      # year of 1993
Tmax$t <- Tmax$julian - 728049   # create a new time variable

## -----------------------------------------------------------
spat_av <- group_by(Tmax, lat, lon) %>%    # group by lon-lat
           summarise(mu_emp = mean(z))     # mean for each lon-lat

## ------------------------------------------------------------------------
lat_means <- ggplot(spat_av) +
             geom_point(aes(lat, mu_emp)) +
             xlab("Latitude (deg)") +
             ylab("Maximum temperature (degF)") + theme_bw()

lon_means <- ggplot(spat_av) +
             geom_point(aes(lon, mu_emp)) +
             xlab("Longitude (deg)") +
             ylab("Maximum temperature (degF)") + theme_bw()

## ------------------------------------------------------------------------
Tmax_av <- group_by(Tmax, date) %>%
           summarise(meanTmax = mean(z))

## ------------------------------------------------------------------------
gTmaxav <-
    ggplot() +
    geom_line(data = Tmax,aes(x = date, y = z, group = id),
              colour = "blue", alpha = 0.04) +
    geom_line(data = Tmax_av, aes(x = date, y = meanTmax)) +
    xlab("Month") + ylab("Maximum temperature (degF)") +
    theme_bw()

## ------------------------------------------------------------------------
lm1 <- lm(z ~ lat + t + I(t^2), data = Tmax) # fit a linear model
Tmax$residuals <- residuals(lm1)             # store the residuals

## ------------------------------------------------------------------------
spat_df <- filter(Tmax, t == 1) %>% # lon/lat coords of stations
            select(lon, lat)  %>%   # select lon/lat only
            arrange(lon, lat)       # sort ascending by lon/lat
m <- nrow(spat_av)                  # number of stations

## ------------------------------------------------------------------------
X <- select(Tmax, lon, lat, residuals, t) %>% # select columns
     spread(t, residuals) %>%                 # make time-wide
     select(-lon, -lat) %>%                   # drop coord info
     t()                                      # make space-wide

## ------------------------------------------------------------------------
Lag0_cov <- cov(X, use = 'complete.obs')
Lag1_cov <- cov(X[-1, ], X[-nrow(X),], use = 'complete.obs')

## ------------------------------------------------------------------------
spat_df$n <- 1:nrow(spat_df)    # assign an index to each station
lim_lon <- range(spat_df$lon)   # range of lon coordinates
lon_strips <- seq(lim_lon[1],   # create 4 long. strip boundaries
                  lim_lon[2],
                  length = 5)
spat_df$lon_strip <- cut(spat_df$lon,     # bin the lon into
                         lon_strips,      # their respective bins
                         labels = FALSE,  # don't assign labels
                         include.lowest = TRUE) # include edges

## ------------------------------------------------------------------------
head(spat_df)   # print the first 6 records of spat_df

plot_cov_strips(Lag0_cov, spat_df)  # plot the lag-0 matrices
plot_cov_strips(Lag1_cov, spat_df)  # plot the lag-1 matrices


## ------------------------------------------------------------------------
data("STObj3", package = "STRbook")
STObj4 <- STObj3[, "1993-07-01::1993-07-31"]

## ------------------------------------------------------------------------
vv <- variogram(object = z~1 + lat, # fixed effect component
                data = STObj4,      # July data
                width = 80,         # spatial bin (80 km)
                cutoff = 1000,      # consider pts < 1000 km apart
                tlags = 0.01:6.01)  # 0 days to 6 days

## ------------------------------------------------------------------------
data("SSTlandmask", package = "STRbook")
data("SSTlonlat", package = "STRbook")
data("SSTdata", package = "STRbook")

## ------------------------------------------------------------------------
delete_rows <- which(SSTlandmask == 1)
SSTdata <- SSTdata[-delete_rows, 1:396]

## ------------------------------------------------------------------------
## Put data into space-wide form
Z <- t(SSTdata)
dim(Z)

## ------------------------------------------------------------------------
## First find the matrix we need to subtract:
spat_mean <- apply(SSTdata, 1, mean)
nT <- ncol(SSTdata)

## Then subtract and standardize:
Zspat_detrend <- Z - outer(rep(1, nT), spat_mean)
Zt <- 1/sqrt(nT - 1)*Zspat_detrend

## ------------------------------------------------------------------------
E <- svd(Zt)

## ------------------------------------------------------------------------
V <- E$v
colnames(E$v) <- paste0("EOF", 1:ncol(SSTdata)) # label columns
EOFs <- cbind(SSTlonlat[-delete_rows, ], E$v)
head(EOFs[, 1:6])

## ------------------------------------------------------------------------
TS <- data.frame(E$u) %>%            # convert U to data frame
      mutate(t = 1:nrow(E$u)) %>%    # add a time field
      gather(EOF, PC, -t)            # put columns (except time)
                                     # into long-table format with
                                     # EOF-PC as key-value pair

## ------------------------------------------------------------------------
TS$nPC <- TS$PC * sqrt(nT-1)

## ------------------------------------------------------------------------
ggplot(EOFs) + geom_tile(aes(x = lon, y = lat, fill = EOF1)) +
   fill_scale(name = "degC") + theme_bw() +
   xlab("Longitude (deg)") + ylab("Latitude (deg)")

## ------------------------------------------------------------------------
nEOF <- 10
EOFset1 <- E$u[1:(nT-7), 1:nEOF] * sqrt(nT - 1)
EOFset2 <- E$u[8:nT, 1:nEOF] * sqrt(nT - 1)

## ------------------------------------------------------------------------
cc <- cancor(EOFset1, EOFset2)  # compute CCA
options(digits = 3)             # print to three d.p.
print(cc$cor[1:5])              # print
print(cc$cor[6:10])

## ------------------------------------------------------------------------
CCA_df <- data.frame(t = 1:(nT - 7),
                    CCAvar1 = (EOFset1 %*% cc$xcoef[,1])[,1],
                    CCAvar2 = (EOFset2 %*% cc$ycoef[,1])[,1])

## ------------------------------------------------------------------------
t_breaks <- seq(1, nT, by = 60)     # breaks for x-labels
year_breaks <- seq(1970,2002,by=5)  # labels for x-axis
g <- ggplot(CCA_df) +
     geom_line(aes(t, CCAvar1), col = "dark blue") +
     geom_line(aes(t, CCAvar2), col = "dark red") +
     scale_x_continuous(breaks = t_breaks, labels = year_breaks) +
     ylab("CCA variables") + xlab("Year") + theme_bw()

## ------------------------------------------------------------------------
EOFs_CCA <- EOFs[,1:4] # first two columns are lon-lat
EOFs_CCA[,3] <- c(as.matrix(EOFs[,3:12]) %*% cc$xcoef[,1])
EOFs_CCA[,4] <- c(as.matrix(EOFs[,3:12]) %*% cc$ycoef[,1])
