# Function to make a nice HTML table for reports using datatable from DT
# Matt Kmiecik
# Started 19 May 2024

custom_datatable <- function(df, round = 2, sci_not_cols = NULL, ...){
  
  # libraries ----
  library(dplyr); library(DT)
  
  # rounds data
  if (is.null(sci_not_cols)) {
    df2 <- df %>% mutate(across(where(is.double), ~round(.x, round)))
  } else{
    df2 <- 
      df %>%
      mutate(
        across(
          .cols = where(is.double) & !matches(sci_not_cols), #not sci not
          ~round(.x, round)
        )
      ) %>%
      mutate(across(.cols = matches(sci_not_cols), ~format(.x, digits = 3)))
  }
  
  # styling
  df2_styled <- datatable(df2, ...)
  
  # returns
  return(df2_styled)
}