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
library("tidyr")
library("STRbook")

## ------------------------------------------------------------------------
locs <- read.table(system.file("extdata", "Stationinfo.dat",
                               package = "STRbook"),
                  col.names = c("id", "lat", "lon"))
Times <- read.table(system.file("extdata", "Times_1990.dat",
                                package = "STRbook"),
                  col.names = c("julian", "year", "month", "day"))
Tmax <- read.table(system.file("extdata", "Tmax_1990.dat",
                               package = "STRbook"))

## ------------------------------------------------------------------------
names(Tmax) <- locs$id

## ------------------------------------------------------------------------
Tmax <- cbind(Times, Tmax)
head(names(Tmax), 10)

## ------------------------------------------------------------------------
Tmax_long <- gather(Tmax, id, z, -julian, -year, -month, -day)
head(Tmax_long)

## ------------------------------------------------------------------------
Tmax_long$id <- as.integer(Tmax_long$id)

## -----------------------------------------------------------
nrow(Tmax_long)
Tmax_long <- filter(Tmax_long, !(z <= -9998))
nrow(Tmax_long)

## ------------------------------------------------------------------------
Tmax_long <- mutate(Tmax_long, proc = "Tmax")
head(Tmax_long)

## ------------------------------------------------------------------------
data(Tmin_long, package = "STRbook")
data(TDP_long, package = "STRbook")
data(Precip_long, package = "STRbook")

## ------------------------------------------------------------------------
NOAA_df_1990 <- rbind(Tmax_long, Tmin_long, TDP_long, Precip_long)

## ------------------------------------------------------------------------
summ <- group_by(NOAA_df_1990, year, proc) %>%  # groupings
        summarise(mean_proc = mean(z))          # operation

## ------------------------------------------------------------------------
NOAA_precip <- filter(NOAA_df_1990, proc == "Precip" & month == 6)
summ <- group_by(NOAA_precip, year, id) %>%
        summarise(days_no_precip = sum(z == 0))
head(summ)

## ------------------------------------------------------------------------
median(summ$days_no_precip)

## -------------------------------------------------------------
grps <- group_by(NOAA_precip, year, id)
summ <- summarise(grps, days_no_precip = sum(z == 0))

## ------------------------------------------------------------------------
NOAA_df_sorted <- arrange(NOAA_df_1990, julian, id)

## ------------------------------------------------------------------------
df1 <- select(NOAA_df_1990, julian, z)
df2 <- select(NOAA_df_1990, -julian)

## ------------------------------------------------------------------------
NOAA_df_1990 <- left_join(NOAA_df_1990, locs, by = "id")

## ------------------------------------------------------------------------
Tmax_long_sel <- select(Tmax_long, julian, id, z)
Tmax_wide <- spread(Tmax_long_sel, id, z)
dim(Tmax_wide)

## ------------------------------------------------------------------------
M <- select(Tmax_wide, -julian) %>% as.matrix()

## -----------------------------------------------------------
library("sp")
library("spacetime")

## ------------------------------------------------------------------------
NOAA_df_1990$date <- with(NOAA_df_1990,
                       paste(year, month, day, sep = "-"))
head(NOAA_df_1990$date, 4)   # show first four elements

## ------------------------------------------------------------------------
NOAA_df_1990$date <- as.Date(NOAA_df_1990$date)
class(NOAA_df_1990$date)

## ------------------------------------------------------------------------
Tmax_long2 <- filter(NOAA_df_1990, proc == "Tmax")
STObj <- stConstruct(x = Tmax_long2,           # data set
                    space = c("lon", "lat"),  # spatial fields
                    time = "date")            # time field
class(STObj)

## ------------------------------------------------------------------------
spat_part <- SpatialPoints(coords = Tmax_long2[, c("lon", "lat")])
temp_part <- Tmax_long2$date
STObj2 <- STIDF(sp = spat_part,
                time = temp_part,
                data = select(Tmax_long2, -date, -lon, -lat))
class(STObj2)

## ------------------------------------------------------------------------
spat_part <- SpatialPoints(coords = locs[, c("lon", "lat")])
temp_part <- with(Times,
                   paste(year, month, day, sep = "-"))
temp_part <- as.Date(temp_part)

## ------------------------------------------------------------------------
Tmax_long3 <- gather(Tmax, id, z, -julian, -year, -month, -day)

## ------------------------------------------------------------------------
Tmax_long3$id <- as.integer(Tmax_long3$id)
Tmax_long3 <- arrange(Tmax_long3,julian,id)

## ------------------------------------------------------------------------
all(unique(Tmax_long3$id) == locs$id)

## ------------------------------------------------------------------------
STObj3 <- STFDF(sp = spat_part,
                time = temp_part,
                data = Tmax_long3)
class(STObj3)

## ------------------------------------------------------------------------
proj4string(STObj3) <- CRS("+proj=longlat +ellps=WGS84")

## ------------------------------------------------------------------------
STObj3$z[STObj3$z == -9999] <- NA

