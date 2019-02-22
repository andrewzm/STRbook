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

library("broom")
library("dplyr")
library("ggplot2")
library("STRbook")
library("purrr")
library("tidyr")

## ------------------------------------------------------------------------
data("SST_df", package = "STRbook")

## ------------------------------------------------------------------------
data("SSTlandmask", package = "STRbook")
data("SSTdata", package = "STRbook")
data("SSTlonlat", package = "STRbook")

## ------------------------------------------------------------------------
lonlatmask_df <- data.frame(cbind(SSTlonlat, SSTlandmask))
names(lonlatmask_df) <- c("lon", "lat", "mask")

## ------------------------------------------------------------------------
SSTdata <- cbind(lonlatmask_df, SSTdata)

## ------------------------------------------------------------------------
SST_df <- gather(SSTdata, date, sst, -lon, -lat, -mask)

## ------------------------------------------------------------------------
SST_df %>% head(3)

## ------------------------------------------------------------------------
date_grid <- expand.grid(Month = c("Jan", "Feb", "Mar", "Apr",
                                   "May", "Jun", "Jul", "Aug",
                                   "Sep", "Oct", "Nov", "Dec"),
                         Year = 1970:2002,
                         stringsAsFactors =  FALSE)
date_grid$date <- paste0("V", 1:396)
SST_df <- left_join(SST_df, date_grid) %>%
          select(-date)

## ------------------------------------------------------------------------
SST_df$date <- paste(SST_df$Month, SST_df$Year)
SST_df %>% head(3)

## ------------------------------------------------------------------------
SST_df$sst<- ifelse(SST_df$mask == 0, SST_df$sst, NA)


## ------------------------------------------------------------------------
g <- ggplot(filter(SST_df, Year == 1997 &  # subset by month/year
                      Month %in% c("Apr","Aug","Jun","Oct"))) +
    geom_tile(aes(lon, lat,
                  fill = pmin(sst, 4))) +  # clamp SST at 4deg
    facet_wrap(~date, dir = "v") +         # facet by date
    fill_scale(limits = c(-4, 4),          # color limits
               name = "degC") +            # legend title
    theme_bw() + coord_fixed()             # fix scale and theme


## ------------------------------------------------------------------------
data("SOI", package = "STRbook")
SOI_df <- select(SOI, -Ann) %>%
          gather(Month, soi, -Year)

## ------------------------------------------------------------------------
SST_df <- left_join(SST_df, SOI_df,
                    by = c("Month", "Year"))


## ------------------------------------------------------------------------
SST_pre_May <- filter(SST_df, Year <= 1997) %>%
               filter(!(Year == 1997 &
                        Month %in% c("May", "Jun", "Jul",
                                     "Aug", "Sep", "Oct",
                                     "Nov", "Dec")))

## ------------------------------------------------------------------------
fit_one_pixel <- function(data)
                 mod <- lm(sst ~ 1 + soi, data = data)

pixel_lms <- SST_pre_May %>%
             filter(!is.na(sst)) %>%
             group_by(lon, lat) %>%
             nest() %>%
             mutate(model = map(data, fit_one_pixel)) %>%
             mutate(model_df = map(model, tidy))

## ------------------------------------------------------------------------
pixel_lms %>% head(3)

## ------------------------------------------------------------------------
lm_pars <- pixel_lms %>%
           unnest(model_df)

## ----echo = FALSE--------------------------------------------------------
options(digits = 3)
options(width = 60)

## ------------------------------------------------------------------------
head(lm_pars, 3)

## ------------------------------------------------------------------------
lm_pars <- left_join(lonlatmask_df, lm_pars)

## ------------------------------------------------------------------------
g2 <- ggplot(filter(lm_pars, term == "(Intercept)" | mask == 1)) +
    geom_tile(aes(lon, lat, fill = estimate)) +
    fill_scale() +
    theme_bw() + coord_fixed()

g3 <- ggplot(filter(lm_pars, term == "soi" | mask == 1)) +
    geom_tile(aes(lon, lat, fill = estimate)) +
    fill_scale() +
    theme_bw() + coord_fixed()

## ------------------------------------------------------------------------
soi_pred <- filter(SOI_df, Month == "Oct" & Year == "1997") %>%
            select(soi)

## ------------------------------------------------------------------------
predict_one_pixel <- function(lm, soi_pred) {
    predict(lm,                           # linear model
            newdata = soi_pred,           # pred. covariates
            interval = "prediction") %>%  # output intervals
    data.frame() %>%                      # convert to df
    mutate(se = (upr-lwr)/(2 * 1.96)) %>% # comp pred. se
    select(fit, se)                       # return fit & se
  }

## ------------------------------------------------------------------------
SST_Oct_1997 <- pixel_lms %>%
                mutate(preds = map(model,
                                   predict_one_pixel,
                                   soi_pred = soi_pred)) %>%
                unnest(preds)
