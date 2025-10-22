# regenerate_plots.R
# Matt Kmiecik
# Purpose: Regenerate plots from saved bandwise analysis results

#' Regenerate all plots from a saved bandwise analysis
#'
#' @param results Either (1) a file path to the saved .rds file from bandwise-analyses.R,
#'   or (2) the loaded results object itself (a list)
#' @param band_name Optional: specific band to plot (e.g., "alpha"). If NULL, plots all bands.
#' @param plot_types Character vector of plot types to generate. Options:
#'   "histogram", "qq", "estimates", "topos_eyes_closed", "topos_eyes_open", 
#'   "topos_interaction", "line_graphs"
#' @param p_fdr Logical: Use FDR-corrected p-values (TRUE) or uncorrected p-values (FALSE).
#'   Default is TRUE. Affects electrode shapes in interaction topographies.
#' @return A list of plot objects organized by band and plot type
#' @examples
#' # Method 1: Provide file path
#' plots <- regenerate_plots("output/r-analysis/bandwise-analyses-db.rds")
#' 
#' # Method 2: Provide the loaded object
#' results <- read_rds("output/r-analysis/bandwise-analyses-db.rds")
#' plots <- regenerate_plots(results)
#' 
#' # Use uncorrected p-values
#' plots <- regenerate_plots(results, p_fdr = FALSE)
#' 
#' # Regenerate all plots for alpha band only
#' plots <- regenerate_plots(results, band_name = "alpha")
#' 
#' # View specific plot
#' plots$alpha$histogram
#' plots$alpha$topos_eyes_closed$eyes_plot
regenerate_plots <- function(results = NULL, 
                             band_name = NULL,
                             plot_types = c("histogram", "qq", "estimates", 
                                            "topos_eyes_closed", "topos_eyes_open",
                                            "topos_interaction", "line_graphs"),
                             p_fdr = TRUE) {
  
  # Libraries ----
  library(tidyverse)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  
  # Functions ----
  source(here::here("fns", "topo_interp.R"))
  source(here::here("fns", "topo_plot.R"))
  
  # Load data ----
  if (is.null(results)) {
    stop("Please provide either a path to results file or the results object itself")
  }
  
  if (is.character(results)) {
    # If a file path is provided, load the file
    cat("Loading results from:", results, "\n")
    res <- read_rds(results)
  } else if (is.list(results)) {
    # If the object itself is provided, use it directly
    cat("Using provided results object\n")
    res <- results
  } else {
    stop("'results' must be either a file path (character) or a results object (list)")
  }
  
  # Determine which bands to process ----
  if (is.null(band_name)) {
    bands_to_process <- names(res)
    cat("Processing all bands:", paste(bands_to_process, collapse = ", "), "\n")
  } else {
    if (!band_name %in% names(res)) {
      stop("Band '", band_name, "' not found in results. Available bands: ", 
           paste(names(res), collapse = ", "))
    }
    bands_to_process <- band_name
    cat("Processing band:", band_name, "\n")
  }
  
  # Initialize output list ----
  all_plots <- vector("list", length(bands_to_process))
  names(all_plots) <- bands_to_process
  
  # Process each band ----
  for (band in bands_to_process) {
    
    cat("\n=== Regenerating plots for", toupper(band), "===\n")
    
    band_data <- res[[band]]
    plots <- list()
    
    # 1. HISTOGRAM ----
    if ("histogram" %in% plot_types) {
      cat("  - Creating power histogram...\n")
      
      plots$histogram <- 
        ggplot(band_data$power_histogram_data, aes(m)) + 
        geom_histogram() +
        geom_vline(
          xintercept = c(
            band_data$histogram_thresholds$lower, 
            band_data$histogram_thresholds$upper
          ), 
          linetype = 2, 
          color = band_data$CONFIG$rdgy[3]
        ) +
        labs(
          x = paste0("PSD (", band_data$unit, ")"), 
          y = "Count", 
          caption = paste0(
            "Lines depict ", 
            band_data$exclusion_threshold$thresh * 100, 
            "% cutoffs."
          )
        ) +
        ggtitle(paste(toupper(band), "- Electrodes x Blocks x Conditions")) +
        theme_bw()
    }
    
    # 2. QQ PLOTS ----
    if ("qq" %in% plot_types) {
      cat("  - Creating QQ plots...\n")
      
      elec_qqplot <- function(df, teyes, velec) {
        elecs <- unique(df$electrode)
        p <-
          ggplot(
            df %>% filter(eyes == teyes, electrode %in% elecs[velec]), 
            aes(sample = resid)
          ) + 
          stat_qq(size = .5) + 
          stat_qq_line(color = "red") +
          theme_bw() +
          labs(
            x = "Theoretical Quantiles", 
            y = "Sample Quantiles", 
            title = paste0(toupper(band), " - Eyes ", teyes)
          ) +
          facet_wrap(~electrode, ncol = 8)
        return(p)
      }
      
      plots$qq_plots <- list(
        eyes_open = elec_qqplot(
          df = band_data$resids, 
          teyes = "open", 
          velec = 1:64
        ),
        eyes_closed = elec_qqplot(
          df = band_data$resids, 
          teyes = "closed", 
          velec = 1:64
        )
      )
    }
    
    # 3. ESTIMATE PLOTS ----
    if ("estimates" %in% plot_types) {
      cat("  - Creating fixed-effect estimate plots...\n")
      
      # Helper function
      plot_fixed <- function(data, tterm, teyes) {
        this_est <- data %>% filter(term == tterm, eyes == teyes)
        tcol <- c("ns" = band_data$CONFIG$rdgy[8], "p<.05" = band_data$CONFIG$rdgy[3])
        
        # Determine which p-value column to use
        sig_col <- if ("p.fdr.sig" %in% names(this_est)) {
          sym("p.fdr.sig")
        } else {
          sym("p.value.sig")
        }
        
        this_plot <- 
          ggplot(
            this_est, 
            aes(estimate, reorder(electrode, estimate), color = !!sig_col)
          ) + 
          geom_vline(xintercept = 0, linetype = 2) +
          geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = .2) +
          geom_point() +
          labs(
            x = "Estimate", 
            y = "Electrode", 
            title = paste(toupper(band), "- Term:", tterm, "\nEyes:", teyes),
            caption = "95% CI error bars."
          ) +
          scale_color_manual(values = tcol) +
          theme_classic()
        return(this_plot)
      }
      
      # Generate all estimate plots
      plots$est_plots <- vector("list", length = nrow(band_data$est_plot_args))
      names(plots$est_plots) <- paste(
        band_data$est_plot_args$term, 
        band_data$est_plot_args$eyes, 
        sep = "."
      )
      
      for (i in 1:nrow(band_data$est_plot_args)) {
        plots$est_plots[[i]] <- 
          plot_fixed(
            band_data$ests, 
            band_data$est_plot_args[["term"]][i], 
            band_data$est_plot_args[["eyes"]][i]
          )
      }
    }
    
    # 4. TOPOGRAPHY PLOTS - EYES CLOSED ----
    if ("topos_eyes_closed" %in% plot_types) {
      cat("  - Creating eyes closed topography plots...\n")
      
      plots$topos_eyes_closed <- regenerate_topo_plots(
        data_orig = band_data$study_sum,
        data_interp = band_data$interp,
        data_sub_orig = band_data$study_sum_sub,
        data_sub_interp = band_data$interp_sub,
        cond = "closed",
        unit = band_data$unit,
        band_name = band,
        CONFIG = band_data$CONFIG
      )
    }
    
    # 5. TOPOGRAPHY PLOTS - EYES OPEN ----
    if ("topos_eyes_open" %in% plot_types) {
      cat("  - Creating eyes open topography plots...\n")
      
      plots$topos_eyes_open <- regenerate_topo_plots(
        data_orig = band_data$study_sum,
        data_interp = band_data$interp,
        data_sub_orig = band_data$study_sum_sub,
        data_sub_interp = band_data$interp_sub,
        cond = "open",
        unit = band_data$unit,
        band_name = band,
        CONFIG = band_data$CONFIG
      )
    }
    
    # 6. INTERACTION TOPOGRAPHIES ----
    if ("topos_interaction" %in% plot_types) {
      cat("  - Creating interaction topography plots...\n")
      
      plots$interaction_topos <- regenerate_interaction_topos(
        dub_sub_df = band_data$dub_sub_df,
        dub_sub_interp = band_data$dub_sub_interp,
        unit = band_data$unit,
        band_name = band,
        CONFIG = band_data$CONFIG,
        p_fdr = p_fdr  # Pass p_fdr argument
      )
    }
    
    # 7. LINE GRAPHS FOR SIGNIFICANT INTERACTIONS ----
    if ("line_graphs" %in% plot_types) {
      cat("  - Creating significant interaction line graphs...\n")
      
      # Choose which significance data to use based on p_fdr
      sig_data_to_use <- if (p_fdr) {
        band_data$sig_data_fdr
      } else {
        band_data$sig_data_uncorrected
      }
      
      if (!is.null(sig_data_to_use) && length(sig_data_to_use) > 0) {
        plots$interaction_line_graphs <- regenerate_line_graphs(
          sig_data_results = sig_data_to_use,
          band_name = band
        )
      }
    }
    
    all_plots[[band]] <- plots
  }
  
  cat("\n✓ Plot regeneration complete!\n")
  return(all_plots)
}

# Helper function to regenerate topography plots ----
regenerate_topo_plots <- function(data_orig, data_interp, data_sub_orig, 
                                  data_sub_interp, cond, unit, band_name, CONFIG) {
  
  # Constants
  nose_adj <- .02
  dia <- 1.45
  
  # Filter for condition
  this_orig <- data_orig %>% filter(eyes == cond)
  this_interp <- data_interp %>% filter(eyes == cond)
  this_sub_orig <- data_sub_orig %>% filter(eyes == cond)
  this_sub_interp <- data_sub_interp %>% filter(eyes == cond)
  
  # Determine color scale for main plot
  minmax <- 
    this_orig %>% 
    group_by(stim, eyes, task) %>%
    summarise(mean = mean(M), min = min(M), max = max(M), .groups = "drop")
  
  min_val <- min(minmax$min)
  max_val <- max(minmax$max)
  n <- if_else(abs(min_val) > abs(max_val), abs(min_val), abs(max_val))
  
  symmetric_vector <- function(magnitude, num_each_side = 3) {
    total_length <- 2 * num_each_side + 1
    return(seq(-magnitude, magnitude, length.out = total_length))
  }
  
  color_pal_breaks <- round(symmetric_vector(n))
  
  # Main eyes plot
  eyes_plot <- 
    topo_plot(
      orig_data = this_orig, 
      interp_data = this_interp, 
      dv = M,
      color_pal_limits = c(min(color_pal_breaks), max(color_pal_breaks)),
      color_pal_breaks = color_pal_breaks,
      elec_shape_col = NULL,
      elec_shapes = 19,
      bwidth = 1.5,
      bheight = .2,
      d = dia,
      contour_alpha = 1/3,
      contour_color = "black",
      headshape_size = .5,
      electrode_size = .25,
      nose_size = .5,
      nose_adj = nose_adj,
      legend_name = paste0("PSD (", unit, ")"),
      elec_loc_path = "doc/chan-locs.rda"
    ) + 
    facet_grid(stim ~ task) +
    labs(title = paste(toupper(band_name), "- Eyes", tools::toTitleCase(cond)))
  
  # Determine color scale for subtraction plot
  this_sub_minmax <- 
    this_sub_orig %>% 
    group_by(stim, eyes, task) %>%
    summarise(mean = mean(M), min = min(M), max = max(M), .groups = "drop")
  
  this_sub_min <- min(this_sub_minmax$min)
  this_sub_max <- max(this_sub_minmax$max)
  nn <- if_else(
    abs(this_sub_min) > abs(this_sub_max), 
    abs(this_sub_min), 
    abs(this_sub_max)
  )
  
  this_sub_color_pal_breaks <- round(symmetric_vector(nn))
  
  # Subtraction plot
  sub_plot <- 
    topo_plot(
      orig_data = this_sub_orig, 
      interp_data = this_sub_interp, 
      dv = M,
      color_pal_limits = c(
        min(this_sub_color_pal_breaks), 
        max(this_sub_color_pal_breaks)
      ),
      color_pal_breaks = this_sub_color_pal_breaks,
      color_pal = rev(brewer.pal(11, "RdBu")),
      elec_shape_col = NULL,
      elec_shapes = 19,
      bwidth = 1.5,
      bheight = .2,
      d = dia,
      contour_alpha = 1/3,
      contour_color = "black",
      headshape_size = .5,
      electrode_size = .25,
      nose_size = .5,
      nose_adj = nose_adj,
      legend_name = paste0("PSD Difference \n (", unit, ")"),
      elec_loc_path = "doc/chan-locs.rda"
    ) + 
    facet_grid(stim ~ task)
  
  # Combined plot
  comb_plot <- eyes_plot | sub_plot
  
  return(list(
    eyes_plot = eyes_plot,
    sub_plot = sub_plot,
    comb_plot = comb_plot
  ))
}

# Helper function to regenerate interaction topographies ----
regenerate_interaction_topos <- function(dub_sub_df, dub_sub_interp, unit, 
                                         band_name, CONFIG, p_fdr = TRUE) {
  
  nose_adj <- .02
  dia <- 1.45
  
  # Select appropriate p-value column based on p_fdr argument
  shape_col_name <- if (p_fdr) {
    "p.fdr.sig"
  } else {
    "p.value.sig"
  }
  
  # Check if the column exists
  if (!shape_col_name %in% names(dub_sub_df)) {
    warning(paste("Column", shape_col_name, "not found in data. Using default shapes."))
    shape_col_name <- NULL
  }
  
  # Dynamic caption
  caption_text <- if (p_fdr) {
    "Closed circles p<.05 (FDR-corrected)"
  } else {
    "Closed circles p<.05 (uncorrected)"
  }
  
  # P-value shapes: open circle for ns, closed circle for significant
  p_value_shapes <- c("ns" = 1, "p<.05" = 19)
  
  # Determine color scale
  minmax <- 
    dub_sub_df %>% 
    group_by(stim, eyes) %>%
    summarise(mean = mean(M_sub), min = min(M_sub), max = max(M_sub), .groups = "drop")
  
  min_val <- min(minmax$min)
  max_val <- max(minmax$max)
  n <- if_else(abs(min_val) > abs(max_val), abs(min_val), abs(max_val))
  
  symmetric_vector <- function(magnitude, num_each_side = 3) {
    total_length <- 2 * num_each_side + 1
    return(seq(-magnitude, magnitude, length.out = total_length))
  }
  
  color_pal_breaks <- round(symmetric_vector(n), 1)
  
  # Create plot
  # Note: We need to use the column name as a symbol for tidy evaluation
  interaction_plot <- 
    topo_plot(
      orig_data = dub_sub_df, 
      interp_data = dub_sub_interp, 
      dv = M_sub,
      color_pal_limits = c(min(color_pal_breaks), max(color_pal_breaks)),
      color_pal_breaks = color_pal_breaks,
      color_pal = rev(brewer.pal(11, "RdBu")),
      elec_shape_col = if (!is.null(shape_col_name)) !!sym(shape_col_name) else NULL,
      elec_shapes = if (!is.null(shape_col_name)) p_value_shapes else 19,
      bwidth = 1.5,
      bheight = .2,
      d = dia,
      contour_alpha = 1/3,
      contour_color = "black",
      headshape_size = .5,
      electrode_size = 1,  # Larger size to see shapes
      nose_size = .5,
      nose_adj = nose_adj,
      legend_name = paste0("PSD Difference \n (", unit, ")"),
      elec_loc_path = "doc/chan-locs.rda"
    ) + 
    facet_grid(stim ~ eyes) +
    labs(
      title = paste(toupper(band_name), "- Interaction (Stim × Task)"),
      caption = if (!is.null(shape_col_name)) caption_text else NULL
    )
  
  return(interaction_plot)
}

# Helper function to regenerate line graphs ----
regenerate_line_graphs <- function(sig_data_results, band_name) {
  
  cpal <- palette.colors(palette = "Okabe-Ito")
  pd <- position_dodge(width = .2)
  
  plots <- map(sig_data_results, function(data) {
    ggplot(data, aes(task, M, group = stim, color = stim)) +
      geom_point(position = pd) +
      geom_errorbar(
        aes(ymin = M - SEM, ymax = M + SEM), 
        width = .1, 
        position = pd
      ) +
      geom_path(position = pd) +
      scale_color_manual(values = cpal) +
      labs(
        x = "Task", 
        y = "Mean PSD",
        title = paste(toupper(band_name), "- Significant Interaction"),
        caption = "Error bars are SEM."
      ) +
      theme_bw() +
      theme(legend.position = "bottom") +
      facet_wrap(~electrode)
  })
  
  return(plots)
}

# Convenience function to save plots ----
save_plots <- function(plots, output_dir = "output/r-analysis/plots", 
                       formats = c("png", "pdf"), width = 10, height = 8) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (band in names(plots)) {
    band_dir <- file.path(output_dir, band)
    if (!dir.exists(band_dir)) {
      dir.create(band_dir, recursive = TRUE)
    }
    
    band_plots <- plots[[band]]
    
    for (plot_type in names(band_plots)) {
      
      # Handle nested lists (like qq_plots, est_plots, etc.)
      if (is.list(band_plots[[plot_type]]) && 
          !inherits(band_plots[[plot_type]], "ggplot")) {
        
        nested_dir <- file.path(band_dir, plot_type)
        if (!dir.exists(nested_dir)) {
          dir.create(nested_dir, recursive = TRUE)
        }
        
        for (sub_name in names(band_plots[[plot_type]])) {
          for (fmt in formats) {
            filename <- file.path(nested_dir, paste0(sub_name, ".", fmt))
            ggsave(
              filename, 
              plot = band_plots[[plot_type]][[sub_name]], 
              width = width, 
              height = height
            )
            cat("Saved:", filename, "\n")
          }
        }
        
      } else if (inherits(band_plots[[plot_type]], "ggplot")) {
        # Save individual plot
        for (fmt in formats) {
          filename <- file.path(band_dir, paste0(plot_type, ".", fmt))
          ggsave(filename, plot = band_plots[[plot_type]], width = width, height = height)
          cat("Saved:", filename, "\n")
        }
      }
    }
  }
  
  cat("\n✓ All plots saved to:", output_dir, "\n")
}