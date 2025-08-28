# EEG Topography Tools
# Matt Kmiecik
# Started 13 April 2024

# HEAVILY inspired by Matt Craddock's code (i.e., credit should go to him):
# https://www.mattcraddock.com/blog/2017/02/25/erp-visualization-creating-topographical-scalp-maps-part-1/
# I also used a similar script in the CRAMPP SSVEP project:
# https://github.com/mkmiecik14/ssvep/blob/main/topo_tools.R

# Function for interpolating
# data is a data frame containing electrode coordinates at x and y
topo_interp <- 
  function(
    data = data, 
    meas = "value", 
    gridRes = 67, 
    size = .85,
    elec_loc_path = here::here("..", "output", "chan-locs.rda")
    ){
    
    # libraries ----
    library(akima); library(tidyverse)
  
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
    
  data_elecs <- left_join(data, elec_locs, by = c("elec" = "labels"))
  #return(data_elecs)
  
  
  
  # procedure ----
    tmp_topo <- 
      with(
        data_elecs,
        akima::interp(
          x = x, 
          y = y, 
          z = data_elecs[[meas]], 
          xo = seq(min(x)*2, max(x)*2, length = gridRes),
          yo = seq(min(y)*2, max(y)*2, length = gridRes),
          linear = FALSE,
          extrap = TRUE
        )
      )
    
    # Creating a matrix that is x by y filled with z
    interp_topo <- data.frame(x = tmp_topo$x, tmp_topo$z)
    names(interp_topo)[1:length(tmp_topo$y) + 1] <- tmp_topo$y
    interp_topo <- 
      interp_topo %>% gather(key = y, value = !!meas, -x, convert = TRUE)
    
    # mark grid elements that are outside of the plotting circle
    # changing the number after the less than sign will modulate the size of the
    # topo drawn
    interp_topo$incircle <- 
      sqrt(interp_topo$x^2 + interp_topo$y^2) < size # original value was .7
    
    # removes the elements outside the circle
    interp_topo <- interp_topo[interp_topo$incircle,] 
    
    # returns interpolated values
    return(interp_topo)
}
