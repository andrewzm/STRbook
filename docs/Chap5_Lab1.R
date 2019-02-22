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
library("ggplot2")
library("STRbook")
set.seed(1)

## ------------------------------------------------------------------------
ds <- 0.01
s_grid <- seq(0, 1, by = ds)
N <- length(s_grid)

## ------------------------------------------------------------------------
nT <- 201
t_grid <- 0:(nT-1)
st_grid <- expand.grid(s = s_grid, t = t_grid)

## ------------------------------------------------------------------------
m <- function(s, x, thetap) {
  gamma <- thetap[1]                 # amplitude
  l <- thetap[2]                     # length scale
  offset <- thetap[3]                # offset
  D <- outer(s + offset, x, '-')     # displacements
  gamma * exp(-D^2/l)                # kernel eval.
}

## ------------------------------------------------------------------------
thetap <- list()
thetap[[1]] <- c(40, 0.0002, 0)
thetap[[2]] <- c(5.75, 0.01, 0)
thetap[[3]] <- c(8, 0.005, 0.1)
thetap[[4]] <- c(8, 0.005, -0.1)

## ------------------------------------------------------------------------
m_x_0.5 <- m(s = 0.5, x = s_grid,            # construct kernel
             thetap = thetap[[1]]) %>%       # at s = 0.5
           as.numeric()                      # convert to numeric
df <- data.frame(x = s_grid, m = m_x_0.5)      # allocate to df
ggplot(df) + geom_line(aes(x, m)) + theme_bw() # plot

## ------------------------------------------------------------------------
mplot <- list()
label <- c("(a)", "(b)", "(c)", "(d)")
for(i in 1:4)  {
  m_x_0.5 <- m(s = 0.5, x = s_grid, thetap = thetap[[i]]) %>% as.numeric()
  df <- data.frame(x = s_grid, m = m_x_0.5)
  mplot[[i]] <- ggplot(df) + geom_line(aes(x, m)) + theme_bw() +
  ggtitle(label[i])
}

## ------------------------------------------------------------------------
Sigma_eta <- 0.1 * exp( -abs(outer(s_grid, s_grid, '-') / 0.1))

## ------------------------------------------------------------------------
L <- t(chol(Sigma_eta))  # chol() returns upper Cholesky factor
sim <- L %*% rnorm(nrow(Sigma_eta))  # simulate

## ------------------------------------------------------------------------
Y <- list()

## ------------------------------------------------------------------------
for(i in 1:4) {                         # for each kernel
  M <- m(s_grid, s_grid, thetap[[i]])   # construct kernel
  Y[[i]] <- data.frame(s = s_grid,      # init. data frame with s
                       t = 0,           # init. time point 0, and
                       Y = 0)           # init. proc. value = 0
  for(j in t_grid[-1]) {                # for each time point
    prev_Y <- filter(Y[[i]],            # get Y at t - 1
                     t == j - 1)$Y
    eta <- L %*% rnorm(N)               # simulate eta
    new_Y <- (M %*% prev_Y * ds + eta) %>%
             as.numeric()               # Euler approximation

    Y[[i]] <- rbind(Y[[i]],             # update data frame
                    data.frame(s = s_grid,
                               t = j,
                               Y =  new_Y))
  }
}

## ------------------------------------------------------------------------
ggplot(Y[[1]]) + geom_tile(aes(s, t, fill = Y)) +
   scale_y_reverse() + theme_bw() +
   fill_scale(name = "Y")

## ------------------------------------------------------------------------
nobs <- 50
sobs <- sample(s_grid, nobs)

## ------------------------------------------------------------------------
Ht <- matrix(0, nobs, N)           # construct empty matrix
for(i in 1:nobs) {                 # for each obs
  idx <- which(sobs[i] == s_grid)  # find the element to set to 1
  Ht[i, idx] <- 1                  # set to 1
}

## ------------------------------------------------------------------------
z_df <- data.frame()               # init data frame
for(j in 0:(nT-1)) {               # for each time point
  Yt <- filter(Y[[1]], t == j)$Y   # get the simulated process
  zt <- Ht %*% Yt + rnorm(nobs)    # map to obs and add noise
  z_df <- rbind(z_df,              # update data frame
                data.frame(s = sobs, t = j, z = zt))
}

## ------------------------------------------------------------------------
ggplot(z_df) + geom_point(aes(s, t, colour = z))  +
  col_scale(name = "z") + scale_y_reverse() + theme_bw()

