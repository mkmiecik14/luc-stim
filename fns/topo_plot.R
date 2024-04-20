# EEG Topography Tools
# Matt Kmiecik
# Started 13 April 2024

# HEAVILY inspired by Matt Craddock's code (i.e., credit should go to him):
# https://www.mattcraddock.com/blog/2017/02/25/erp-visualization-creating-topographical-scalp-maps-part-1/
# I also used a similar script in the CRAMPP SSVEP project:
# https://github.com/mkmiecik14/ssvep/blob/main/topo_tools.R

topo_plot <- function(
    orig_data, # original data
    interp_data, # interpolated data
    dv, # name of the column that has the values to plot
    d = 1.7, # diameter of the maskRing = circleFun(diameter = 1.7), # creates the mask
    contour_alpha = 1/3, # alpha level of contour lines
    contour_color = "black", # color of contour lines
    headshape_size = .25, # headshape size
    electrode_size = 1, # size of electrode points
    nose_size = .25, # size of nose shape
    bwidth = .5, # width of colorbar
    bheight = .1, # height of colorbar
    color_pal =  NULL, # color palette: brewer.pal(n = 9, "Purples")
    color_pal_limits = c(-30, -10),
    color_pal_breaks = c(-30, -20, -10),
    legend_name = "DV",
    elec_shape_col = NULL,
    elec_shapes = NULL,
    nose_adj = -.1,
    elec_loc_path = "../output/chan-locs.rda"
    ){
  
  library(dplyr); library(ggplot2); library(RColorBrewer); library(scales)
  
  # Electrode Locations
  # Loads in electrode positions and converts from polar to cartesian
  # note: use little x and little y for the coordinates, not X and Y from MATLAB
  load(elec_loc_path) # change this to rds so you can assign
  elec_locs <- 
    chan_locs %>%
    mutate(
      radianTheta = pi/180*theta, 
      x = radius*sin(radianTheta), 
      y = radius*cos(radianTheta)
    )
  
  # joins original data with electrode positions
  data_elecs <- left_join(orig_data, elec_locs, by = c("elec" = "labels"))
  
  # theme_topo()
  # Plotting tools for the topo map with head shape:
  theme_topo <- 
    function(base_size = 12){
      theme_bw(base_size = base_size) %+replace%
        theme(
          rect = element_blank(),
          line = element_blank(), 
          axis.text = element_blank(),
          axis.title = element_blank()
        )
    }
  
  # circleFun()
  # function for drawing head for topo plots
  circleFun <- 
    function(center = c(0,0), diameter = 1, npoints = 100){
      r = diameter / 2
      tt <- seq(0, 2*pi, length.out = npoints)
      xx <- center[1] + r * cos(tt)
      yy <- center[2] + r * sin(tt)
      return(data.frame(x = xx, y = yy))
    }
  maskRing <- circleFun(diameter = d)
  
  # Headshape coordinates are determined
  headShape <- 
    circleFun(
      c(0, 0), 
      diameter = max(elec_locs$x) - min(elec_locs$x), #round(max(elec_locs$x))
      npoints = 100
    )
  
  # Nose coordinates
  nose <- data.frame(x = c(-0.075, 0, .075), y = c(.495, .575, .495) + nose_adj) 
  
  # color palette
  if (is.null(color_pal)) {
    color_pal <- brewer.pal(n = 9, "Purples")
  } else {
    # do nothing and use the color pal provided
  }
  
  # Produces blank topography map with electrodes
  # ggplot(headShape, aes(x, y)) +
  #   geom_path() +
  #   geom_text(data = elec_locs, aes(x, y, label = labels)) +
  #   geom_line(data = nose, aes(x, y, z = NULL)) +
  #   theme_topo() +
  #   coord_equal()
  plot <- 
    ggplot(interp_data, aes(x = x, y = y, fill = {{dv}}, group = 1)) +
    coord_equal() + # equalizes coordinates
    geom_raster(interpolate = TRUE) + # basis of the topo
    stat_contour(aes(z = {{dv}}), colour = contour_color, alpha = contour_alpha) + # contour lines
    # creates a mask to remove points outside defined circle
    geom_path(
      data = maskRing,
      aes(x, y, z = NULL, fill = NULL),
      colour = "white",
      size = 6
    ) +
    theme_topo() + # topo theme is added (white background etc.)
    # plots headshape
    geom_path(data = headShape, aes(x, y, z = NULL, fill = NULL), size = headshape_size) +
    # plots nose
    geom_path(data = nose, aes(x, y, z = NULL, fill = NULL), size = nose_size) +
    # colors here
    # note: oob = squish forces everything outside the colour limits to equal
    # nearest colour boundary (i.e., below min colours = min colour)
    scale_fill_gradientn(
      colours = color_pal,
      limits = color_pal_limits, # these should be determined from the uninterpolated data
      breaks = color_pal_breaks, labels = color_pal_breaks,
      guide = "colourbar",
      oob = squish,
      name = legend_name
    ) + 
    geom_point(
      data = data_elecs, 
      aes(x, y, shape = {{elec_shape_col}}),
      size = electrode_size
      ) +
    scale_shape_manual(values = elec_shapes) +
    guides(
      shape = "none", 
      fill = guide_colourbar(
        title.position = "top", 
        title.hjust = .5, 
        frame.colour = "black", 
        ticks.colour = "black", 
        barwidth = unit(bwidth, "in"),
        barheight = unit(bheight, "in"),
        theme = theme(legend.title.align = .5)
      )
    ) +
    theme(legend.position = "bottom")
      
 return(plot)
  
}

