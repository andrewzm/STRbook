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

library("animation")
library("dplyr")
library("ggplot2")
library("gstat")
library("maps")
library("STRbook")
library("grid")
library("gridExtra")

## ------------------------------------------------------------------------
set.seed(1)

## ----message=FALSE, warning = FALSE--------------------------------------
data("NOAA_df_1990", package = "STRbook")
Tmax <- filter(NOAA_df_1990,     # subset the data
              proc == "Tmax" &   # only max temperature
              month %in% 5:9 &   # May to September
              year == 1993)      # year of 1993

## -----------------------------------------------------------
Tmax %>% select(lon, lat, date, julian, z) %>% head()

## -----------------------------------------------------------
Tmax$t <- Tmax$julian - 728049     # create a new time variable

## -----------------------------------------------------------
Tmax_1 <- subset(Tmax, t %in% c(1, 15, 30))  # extract data

## -----------------------------------------------------------
NOAA_plot <- ggplot(Tmax_1) +             # plot points
    geom_point(aes(x = lon,y = lat,       # lon and lat
                   colour = z),           # attribute color
               size = 2) +                # make all points larger
    col_scale(name = "degF") +            # attach color scale
    xlab("Longitude (deg)") +             # x-axis label
    ylab("Latitude (deg)") +              # y-axis label
    geom_path(data = map_data("state"),   # add US states map
          aes(x = long, y = lat, group = group)) +
    facet_grid(~date) +                   # facet by time
    coord_fixed(xlim = c(-105, -75),
                ylim = c(25, 50))  +      # zoom in
    theme_bw()                            # B&W theme

## ------------------------------------------------------------------------
data("BEA", package = "STRbook")
head(BEA %>% select(-Description), 3)

## ------------------------------------------------------------------------
data("MOcounties", package = "STRbook")
head(MOcounties %>% select(long, lat, NAME10), 3)

## -----------------------------------------------------------
County1 <- filter(MOcounties, NAME10 == "Clark, MO")
plot(County1$long, County1$lat)

## ------------------------------------------------------------------------
MOcounties <- left_join(MOcounties, BEA, by = "NAME10")

## -----------------------------------------------------------
g1 <- ggplot(MOcounties) +
    geom_polygon(aes(x = long, y = lat,     # county boundary
                     group = NAME10,        # county group
                     fill = log(X1970))) +  # log of income
    geom_path(aes(x = long, y = lat,        # county boundary
                  group = NAME10)) +        # county group
    fill_scale(limits = c(7.5,10.2),
               name = "log($)")  +
    coord_fixed() + ggtitle("1970") +       # annotations
    xlab("x (m)") + ylab("y (m)") + theme_bw()

## -----------------------------------------------------------
g2 <- ggplot(MOcounties) +
    geom_polygon(aes(x=long,y=lat,group=id,fill=log(X1980))) +
    geom_path(aes(x=long,y=lat,group=id)) +
    scale_fill_distiller(palette = "Spectral",limits=c(7.5,10.2),name="log($)")  +
    coord_fixed()  + ggtitle("1980") + xlab("x (m)")  + theme_bw() +
    ylab("y (m)")

g3 <- ggplot(MOcounties) +
    geom_polygon(aes(x=long,y=lat,group=id,fill=log(X1990))) +
    geom_path(aes(x=long,y=lat,group=id)) +
    scale_fill_distiller(palette = "Spectral",limits=c(7.5,10.2),name="log($)")  +
    coord_fixed() + ggtitle("1990")   + xlab("x (m)")  + theme_bw() +
    ylab("y (m)")

## -----------------------------------------------------------
UIDs <- unique(Tmax$id)                     # extract IDs
UIDs_sub <- sample(UIDs, 10)                # sample 10 IDs
Tmax_sub <- filter(Tmax, id %in% UIDs_sub)  # subset data

## ------------------------------------------------------------------------
TmaxTS <- ggplot(Tmax_sub) +
    geom_line(aes(x = t, y = z)) + # line plot of z against t
    facet_wrap(~id, ncol = 5) +    # facet by station
    xlab("Day number (days)") +    # x label
    ylab("Tmax (degF)") +          # y label
    theme_bw() +                   # BW theme
    theme(panel.spacing = unit(1, "lines")) # facet spacing

## ------------------------------------------------------------------------
lim_lat <- range(Tmax$lat)        # latitude range
lim_t <- range(Tmax$t)            # time range
lat_axis <- seq(lim_lat[1],       # latitude axis
                lim_lat[2],
                length=25)
t_axis <- seq(lim_t[1],           # time axis
              lim_t[2],
              length=100)
lat_t_grid <- expand.grid(lat = lat_axis,
                          t = t_axis)

## -----------------------------------------------------------
Tmax_grid <- Tmax
dists <- abs(outer(Tmax$lat, lat_axis, "-"))
Tmax_grid$lat <- lat_axis[apply(dists, 1, which.min)]

## ------------------------------------------------------------------------
Tmax_lat_Hov <- group_by(Tmax_grid, lat, t) %>%
                summarise(z = mean(z))

## ------------------------------------------------------------------------
Hovmoller_lat <- ggplot(Tmax_lat_Hov) +            # take data
        geom_tile(aes(x = lat, y = t, fill = z)) + # plot
        fill_scale(name = "degF") +     # add color scale
        scale_y_reverse() +             # rev y scale
        ylab("Day number (days)") +     # add y label
        xlab("Latitude (degrees)") +    # add x label
        theme_bw()                      # change theme

## -----------------------------------------------------------
station <- 13966
station_lon <- filter(Tmax,id == 13966)$lon[1]
station_lat <- filter(Tmax,id == 13966)$lat[1]
Hovmoller_lat <- Hovmoller_lat + geom_vline(xintercept = station_lat,linetype='dashed')

## -----------------------------------------------------------
Tmax_t <- function(tau) {
    Tmax_sub <- filter(Tmax, t == tau)        # subset data
    ggplot(Tmax_sub) +
        geom_point(aes(x = lon,y = lat, colour = z),   # plot
                   size = 4) +                         # pt. size
        col_scale(name = "z", limits = c(40, 110)) +
        theme_bw() # B&W theme
}

## -----------------------------------------------------------
help(saveHTML)

## -----------------------------------------------------------
gen_anim <- function() {
   for(t in lim_t[1]:lim_t[2]){  # for each time point
      plot(Tmax_t(t))            # plot data at this time point
   }
}
 
ani.options(interval = 0.2)     # 0.2s interval between frames
saveHTML(gen_anim(),            # run the main function
        autoplay = FALSE,      # do not play on load
        loop = FALSE,          # do not loop
        verbose = FALSE,       # no verbose
        outdir = ".",          # save to current dir
        single.opts = "'controls': ['first', 'previous',
                                    'play', 'next', 'last',
                                     'loop', 'speed'],
                                     'delayMin': 0",
        htmlfile = "NOAA_anim.html")  # save filename
