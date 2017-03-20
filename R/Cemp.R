# ## I'm doing it this way because it's really awkward with NAs to deal with lags...
# ## DEPRECATED
# ## Function C_emp: compute the lagged empirical spatial 
# ## covariance matrix for a lag tau
# ## 'df' must contain fields lon,lat,t and residuals
# ## 'spat_df' must contain fields lon and lat
# ## 'tau' must be an integer
# C_emp <- function(df, spat_df, tau=0) {
#     m <- nrow(spat_df)
#     lim_t <- range(df$t)
#     C_Z <- matrix(0,m,m)                   # running-sum matrix
#     count <- 0                             # number of matrices computed
#     for(tt in (tau + lim_t[1]):lim_t[2]) { # for each time point
#         
#         df1 <-  subset(df,t == tt) %>%       # subset at first time point
#             select(lat,lon,residuals)    # select fields
#         df2 <-  subset(df, t== tt - tau) %>% # subset at lagged time point
#             select(lat,lon,residuals)    # select fields
#         
#         if(nrow(df1) == m &  nrow(df2) == m) {               
#             resid_df <- spat_df %>%                    # take spat_df
#                 left_join(df1,by=c("lon","lat")) %>%   # merge with first subset
#                 mutate(resid1 = residuals) %>%         # create new column resid1
#                 select(-residuals) %>%                 # remove column residuals
#                 left_join(df2,by=c("lon","lat")) %>%   # merge with second subset
#                 mutate(resid2 = residuals)             # create new column resid2
#             C_Z <- C_Z + outer(resid_df$resid1,        # add outer product to C_Z
#                                resid_df$resid2)
#             count <- count + 1                         # increment count
#         }
#     }
#     DD <- C_Z / count       # compute average
# }

# plot subsets of the covariance matrix M
#' @export
plot_cov_strips <- function(C,spat_df) {  
    # for each longitudinal strip 
    require(fields)   # load fields for plotting
    for(i in seq_along(unique(spat_df$lon_strip))){                        
        spat_strip <- spat_df %>%         # take spat_df
            filter(lon_strip == i)  %>%   # extract the ith strip
            arrange(lat)                  # sort by latitude
        idx <- spat_strip$n               # extract indices of locations
        jitter <- seq(0,0.0001,           # add jitter for locations that
                      length=length(idx)) # have same latitude component
        image.plot(spat_strip$lat+jitter, # plot the matrix using fields
                   spat_strip$lat+jitter,
                   C[idx,idx],            # subset and permute C
                   xlab="latitude",
                   ylab="latitude",
                   zlim=c(-15,85),
                   col=tim.colors(10),
                   cex=200)
    }
}