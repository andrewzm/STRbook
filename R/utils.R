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


col_scale <- function(leg_title) {
    scale_colour_distiller(palette="Spectral",   # spectral colour scale
                           guide="colourbar",    # continuous colour bar
                           name=leg_title)
}

fill_scale <- function(leg_title) {
    scale_fill_distiller(palette="Spectral",   # spectral colour scale
                           guide="colourbar",    # continuous colour bar
                           name=leg_title)
}

