#' @export
LinePlotTheme <- function (df)
{
    g <- ggplot(df) + theme(panel.background = element_rect(fill = "white",
                                                            colour = "black"),
                            text = element_text(size = 20),
                            panel.grid.major = element_line(colour = "light gray",size = 0.05),
                            panel.border = element_rect(fill = NA,colour = "black"),
                            plot.margin = unit(c(5, 5, 5, 0), "mm"))
    return(g)
}

# world_map <- geom_path(data=map_data("world"),
#           aes(x=long,y=lat,group=group))

#' @export
col_scale <- function(palette = "Spectral", name = "", limits = NULL) {
    scale_colour_distiller(palette = palette,   # spectral colour scale
                           guide = "colourbar", # continuous colour bar
                           name = name,
                           limits = limits)
}

#' @export
fill_scale <- function(palette = "Spectral", name = "", limits = NULL) {
    scale_fill_distiller(palette = palette,   # spectral colour scale
                         guide = "colourbar", # continuous colour bar
                         name = name,
                         limits = limits)
}

## #' @name dist-matrix
## #' @title Distance Matrix Computation from Two Matrices
## #' @description This function extends \code{dist} to accept two arguments.
## #' @param x1 matrix of size N1 x n
## #' @param x2 matrix of size N2 x n
## #' @details Computes the distances between the coordinates in \code{x1} and the coordinates in \code{x2}. The matrices \code{x1} and \code{x2} do not need to have the same number of rows, but need to have the same number of columns (dimensions).
## #' @return Matrix of size N1 x N2
## #' @export
## #' @examples
## #' A <- matrix(rnorm(50),5,10)
## #' D <- distR(A,A[-3,])
## distR <- function (x1, x2 = NULL)  {
##     ## Try to coerce to matrix
##     if (!is.matrix(x1)) {
##         x1 <- as.matrix(x1)
##     }

##     ## If x2 is not specified set it equatl to x1
##     if (is.null(x2)) {
##         x2 <- x1
##     }

##     ## If it is specified, coerce it to matrix
##     if (!is.matrix(x2)) {
##         x2 <- as.matrix(x2)
##     }

##     ## Basic check
##     if(!(ncol(x1) == ncol(x2)))
##         stop("x1 and x2 have to have same number of columns")

##     ## Compute the distance in C (distR_C is a wrapper)
##     distR_C(x1,x2)
## }
