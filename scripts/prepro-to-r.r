# prepro-to-r.R
# Matt Kmiecik
# Purpose: organizes data from EEG spectral analysis to format suitable for R

# Packages ----
library(R.matlab)
library(tidyverse)
library(readxl)
library(purrr)

# Functions ----
source("fns/extract_spectral_data.R")

# Configuration ----
CONFIG <- list(
  stim_table_file = Sys.getenv("STIM_TABLE_FILE", "doc/ss-info.xlsx"),
  stim_table_sheet = "stimtable",
  spectral_dir = "output/spectral",
  spectral_pattern = "*spectral.mat",
  output_dir = "output/r-prepro",
  data_vars = c("psd_db", "psd_uv", "paf", "cog", "iaf")
)

# Create output directory
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# Helper Functions ----
load_stimulation_table <- function(file_path, sheet_name) {
  if (!file.exists(file_path)) {
    stop("Stimulation table file not found: ", file_path)
  }
  
  cat("Loading stimulation table from:", file_path, "\n")
  read_excel(file_path, sheet = sheet_name) %>%
    select(ss, session, stim_type, stim_version) %>%
    mutate(
      stim_type = case_when(
        stim_type == 1 ~ "tdcs", 
        stim_type == 3 ~ "tacs", 
        stim_type == 4 ~ "trns"
      ),
      stim = if_else(stim_version == "B", "sham", stim_type),
      stim = factor(stim), 
      stim = relevel(stim, ref = "sham")
    )
}

load_spectral_files <- function(spectral_dir, pattern) {
  if (!dir.exists(spectral_dir)) {
    stop("Spectral directory not found: ", spectral_dir)
  }
  
  files <- list.files(spectral_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No spectral files found in ", spectral_dir, " with pattern ", pattern)
  }
  
  cat("Found", length(files), "spectral files\n")
  return(files)
}

extract_all_spectral_data <- function(spec_files) {
  cat("Extracting spectral data from", length(spec_files), "files...\n")
  
  res <- spec_files %>%
    set_names(basename(.) %>% str_remove("-spectral\\.mat$")) %>%
    map(safely(extract_spectral_data)) %>%
    transpose()
  
  # Check for extraction failures
  failed <- map_lgl(res$error, ~ !is.null(.x))
  if (any(failed)) {
    failed_files <- names(res$error)[failed]
    warning("Failed to extract data from ", length(failed_files), " files: ", 
            paste(failed_files, collapse = ", "))
  }
  
  # Extract successful results
  successful_res <- res$result[!failed]
  names(successful_res) <- map_chr(successful_res, "subject")
  
  cat("Successfully extracted data from", length(successful_res), "files\n")
  return(successful_res)
}

organize_and_save_data <- function(res, data_vars, stim_table, output_dir) {
  cat("Organizing and saving data variables...\n")
  
  # Determine available data variables from the first result
  available_vars <- intersect(data_vars, names(res[[1]]))
  if (length(available_vars) < length(data_vars)) {
    missing_vars <- setdiff(data_vars, available_vars)
    warning("Some data variables not found: ", paste(missing_vars, collapse = ", "))
  }
  
  # Process each data variable
  results <- available_vars %>%
    set_names() %>%
    map(~ {
      cat(sprintf("Processing %s...\n", .x))
      
      org_data <- res %>%
        map(.x) %>%
        list_rbind(names_to = "ss") %>%
        left_join(stim_table, by = "ss")
      
      output_file <- file.path(output_dir, paste0(.x, ".rds"))
      write_rds(org_data, output_file)
      
      cat(sprintf("Saved %s to %s\n", .x, output_file))
      return(org_data)
    })
  
  cat("Data processing complete!\n")
  return(results)
}

# Main Processing Pipeline ----
stim_table <- load_stimulation_table(CONFIG$stim_table_file, CONFIG$stim_table_sheet)
spec_files <- load_spectral_files(CONFIG$spectral_dir, CONFIG$spectral_pattern)
res <- extract_all_spectral_data(spec_files)

# Process and save all data variables
results <- organize_and_save_data(res, CONFIG$data_vars, stim_table, CONFIG$output_dir)
