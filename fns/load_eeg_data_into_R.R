# load_eeg_data_into_R.R
# Matt Kmiecik

# helper function to load data
load_eeg_data_into_R <- function(data_file = NULL, refresh = FALSE){
  if (refresh) {
    cat("Refreshing data...\n")
    source("scripts/prepro-to-r.r")
  } else{
    cat("Using previously saved data.\n")
  }
  res <- readr::read_rds(file = data_file)
  return(res)
}