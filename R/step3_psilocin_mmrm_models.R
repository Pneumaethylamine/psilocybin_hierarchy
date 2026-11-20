# ---- paths: set these environment variables (see README) ----
project_dir <- Sys.getenv("UPF_PROJECT_DIR")  # project root (data/, outputs/)
np2_dir     <- Sys.getenv("NP2_DATA_DIR")     # NeuroPharm2 derived data (restricted access)
roi_dir     <- Sys.getenv("ROI_DIR")          # Schaefer atlas files
# Participant(s) omitted from the lines of the exploratory hierarchy-timeline plot (not in the paper);
# comma-separated IDs, kept out of the code
timeline_plot_exclude <- strsplit(Sys.getenv("TIMELINE_PLOT_EXCLUDE"), ",")[[1]]

# ============================================================================
# STREAMLINED R ANALYSIS PIPELINE FOR PSILOCIN STUDY
# Updated with merged hierarchy analyses and improved visualizations
# ============================================================================

# -----------------------------------------------------------------------------
# 1. SETUP AND DATA LOADING
# -----------------------------------------------------------------------------

# Load required libraries
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(ggplot2)
library(scales)
library(MuMIn)  # for R2 calculation
require(LMMstar)
require(mmrm)
library(readr)

# Set paths - find most recent output folder
base_out_dir <- file.path(project_dir, "outputs")

# Get all date folders and find the most recent one
date_folders <- list.dirs(base_out_dir, full.names = FALSE, recursive = FALSE)
date_folders <- date_folders[grepl("^\\d{4}-\\d{2}-\\d{2}$", date_folders)]  # Only YYYY-MM-DD format

if(length(date_folders) == 0) {
  stop("No date-specific output folders found in ", base_out_dir)
}

# Sort and get most recent
most_recent_date <- sort(date_folders, decreasing = TRUE)[2] ## change this
out_dir <- file.path(base_out_dir, most_recent_date)
fig_dir <- file.path(out_dir, "figures")

cat(sprintf("Using most recent output folder: %s\n", most_recent_date))
cat(sprintf("Output directory: %s\n\n", out_dir))

# Ensure directories exist
if(!dir.exists(out_dir)) {
  stop(sprintf("Output directory does not exist: %s", out_dir))
}
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Reading data from: %s\n", out_dir))
cat(sprintf("Saving figures to: %s\n\n", fig_dir))

# Define color schemes and themes
network_order <- c("Default", "Cont", "DorsAttn", "SalVentAttn", "Vis", "SomMot", "Limbic", "SubC")

network_colors <- c(
  "Default" = "#1f77b4",    # DMN - classic dark blue
  "Cont" = "#d62728",       # Control - bright red
  "DorsAttn" = "#ff7f0e",   # Dorsal Attention - bright orange
  "SalVentAttn" = "#e377c2", # Salience/Ventral Attention - pink/magenta
  "Vis" = "#8e9e8b",        # Visual - muted grey-green
  "SomMot" = "#b5b8c9",     # Somatomotor - pale grey-blue
  "Limbic" = "#bcbd22",     # Limbic - muted yellow-green
  "SubC" = "#a0a0a0"        # Subcortical - medium grey
)

custom_theme <- theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold")
  )

theme_set(theme_gray(base_size = 16))

# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPROCESSING FUNCTIONS
# -----------------------------------------------------------------------------

load_and_prepare_data <- function() {
  # Load combined scanner metrics from MATLAB
  combined_data <- read_xlsx(file.path(out_dir, "combined_scanner_metrics.xlsx"))
  
  # Remove old modularity columns if they exist
  combined_data <- combined_data %>%
    select(-any_of(c("modularity_fc", "modularity_forward_MI", 
                     "modularity_fc_abs", "modularity_forward_MI_abs")))
  
  # Load new modularity results
  cat("Loading new modularity results...\n")
  mod_fc <- read_csv(file.path(out_dir, "modularity_summary_fc.csv"))
  mod_fmi <- read_csv(file.path(out_dir, "modularity_summary_fmi.csv"))
  
  # Convert CIMBI_ID to character for joining
  mod_fc <- mod_fc %>% mutate(CIMBI_ID = as.character(CIMBI_ID))
  mod_fmi <- mod_fmi %>% mutate(CIMBI_ID = as.character(CIMBI_ID))
  
  # Merge modularity data with main data
  combined_data <- combined_data %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
    left_join(
      mod_fc %>% select(CIMBI_ID, time_since_admin, modularity_abs, modularity_norm) %>%
        rename(modularity_fc_abs = modularity_abs, modularity_fc = modularity_norm),
      by = c("CIMBI_ID", "time_since_admin")
    ) %>%
    left_join(
      mod_fmi %>% select(CIMBI_ID, time_since_admin, modularity_abs, modularity_norm) %>%
        rename(modularity_forward_MI_abs = modularity_abs, modularity_forward_MI = modularity_norm),
      by = c("CIMBI_ID", "time_since_admin")
    )
  
  # Convert to tibble and handle factors
  combined_data <- combined_data %>%
    mutate(
      CIMBI_ID = as.factor(CIMBI_ID),
      Sex = as.factor(Sex),
      scanner = as.factor(scanner)
    )
  
  # Calculate time_since_baseline
  combined_data <- combined_data %>%
    group_by(CIMBI_ID) %>%
    mutate(time_since_baseline = time_since_admin - min(time_since_admin)) %>%
    ungroup()
  
  # Create time windows
  combined_data <- combined_data %>%
    mutate(time_window = case_when(
      time_since_admin < 3 ~ "predrug",
      time_since_admin >= 47 & time_since_admin < 95 ~ "upslope",
      time_since_admin >= 95 & time_since_admin < 140 ~ "peak",
      time_since_admin >= 140 & time_since_admin < 250 ~ "downslope",
      time_since_admin >= 282 ~ "return",
      TRUE ~ "other"
    )) %>%
    mutate(time_window = factor(time_window, 
                                levels = c("predrug", "upslope", "peak", "downslope", "return", "other")))
  
  # Remove "other" time window scans
  combined_data <- combined_data %>%
    filter(time_window != "other")
  
  # Keep only first scan per subject per time window
  combined_data <- combined_data %>%
    group_by(CIMBI_ID, time_window) %>%
    arrange(time_since_admin) %>%
    slice(1) %>%
    ungroup()
  
  # Assign scan numbers
  combined_data <- combined_data %>%
    mutate(scan = case_when(
      time_window == "predrug" ~ 0,
      time_window == "upslope" ~ 1,
      time_window == "peak" ~ 2,
      time_window == "downslope" ~ 3,
      time_window == "return" ~ 4
    ))
  
  # Scale variables
  combined_data <- combined_data %>%
    mutate(
      psi_conc_scaled = scale(psi_conc),
      Age_scaled = scale(Age)
    )
  
  # Calculate occupancy
  combined_data$occupancy <- calculate_occupancy(combined_data$psi_conc)
  combined_data$occupancy_ofmax <- (combined_data$occupancy/0.766)
  
  # Load FD data
  copbet_path <- file.path(np2_dir, 'CopBET/')
  copbet_files <- list.files(path = copbet_path, pattern = "^CopBETtbl.*\\.csv$", full.names = TRUE)
  
  if(length(copbet_files) > 0) {
    cat("Loading FD data from CopBET table...\n")
    fd_data <- read.csv(copbet_files[1])
    
    fd_subset <- fd_data %>% 
      select(CIMBI, time_since_admin, FD) %>%
      rename(CIMBI_ID = CIMBI) %>%
      mutate(CIMBI_ID = as.character(CIMBI_ID))
    
    combined_data <- combined_data %>%
      mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
      left_join(fd_subset, by = c("CIMBI_ID", "time_since_admin"))
    
    cat(sprintf("FD data merged: %d scans with FD values\n", sum(!is.na(combined_data$FD))))
  } else {
    warning("No CopBET files found - FD data not available")
    combined_data$FD <- NA
  }
  
  # Print diagnostic info
  cat("\n=== SCAN STRUCTURE CHECK ===\n")
  
  duplicates <- combined_data %>%
    group_by(CIMBI_ID, time_window) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if(nrow(duplicates) > 0) {
    cat("WARNING: Found duplicates within time windows:\n")
    print(duplicates)
  } else {
    cat("✓ No duplicates found within time windows\n")
  }
  
  scan_summary <- combined_data %>%
    group_by(CIMBI_ID) %>%
    summarise(
      n_scans = n(), 
      scan_numbers = paste(sort(scan), collapse = ","),
      time_windows = paste(time_window, collapse = ",")
    )
  
  cat("\nScans per subject:\n")
  print(scan_summary)
  
  # Print available variables
  cat("\n=== AVAILABLE VARIABLES ===\n")
  cat("Modularity variables:\n")
  mod_vars <- grep("modularity", names(combined_data), value = TRUE)
  print(mod_vars)
  cat("\nOther global measures:\n")
  other_vars <- c("insideout", "global_GCOR")
  print(intersect(other_vars, names(combined_data)))
  
  return(combined_data)
}

calculate_occupancy <- function(psi_conc, occmax = 76.6, ec50 = 1.95) {
  occmax * (psi_conc / (ec50 + psi_conc))
}

# -----------------------------------------------------------------------------
# 3. ANALYSIS FUNCTIONS
# -----------------------------------------------------------------------------

analyze_timeline_variable <- function(data, var_name, nice_name) {
  
  data <- data %>%
    mutate(time_window = factor(time_window, 
                                levels = c("predrug", "upslope", "peak", "downslope", "return")))
  
  summary_stats <- data %>%
    filter(!is.na(time_window)) %>%
    group_by(time_window) %>%
    summarise(
      mean_value = mean(.data[[var_name]], na.rm = TRUE),
      se_value = sd(.data[[var_name]], na.rm = TRUE) / sqrt(n()),
      n_obs = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      lower_ci = mean_value - 1.96 * se_value,
      upper_ci = mean_value + 1.96 * se_value
    )
  
  plot <- ggplot(data %>% filter(!is.na(time_window)), aes(x = time_window, y = .data[[var_name]])) +
    geom_point(aes(colour = time_since_admin), alpha = 0.4, size = 2) +
    geom_line(aes(group = CIMBI_ID), alpha = 0.15, color = "grey60", linewidth = 0.3) +
    geom_point(data = summary_stats, aes(x = time_window, y = mean_value), 
               size = 6, shape = 21, fill = "black", color = "black") +
    scale_colour_gradientn(
      colors = c("black", "purple", "darkorange3"),
      values = scales::rescale(c(0, 160, max(data$time_since_admin, na.rm = TRUE))),
      guide = "none"
    ) +
    {if(var_name == "insideout") scale_y_continuous(labels = function(x) {
      ifelse(x == 0, "0", {
        exponent <- floor(log10(abs(x)))
        coefficient <- x / 10^exponent
        parse(text = paste0(round(coefficient, 1), " %*% 10^", exponent))
      })
    })} +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 18),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "Time Window", y = nice_name)
  
  # Statistical comparison
  predrug_data <- data %>%
    filter(time_window == "predrug") %>%
    select(CIMBI_ID, all_of(var_name)) %>%
    rename_with(~ "predrug_value", all_of(var_name))
  
  time_windows <- c("upslope", "peak", "downslope", "return")
  results_list <- list()
  
  cohens_d_paired <- function(x, y) {
    diff = x - y
    d = mean(diff) / sd(diff)
    return(d)
  }
  
  for(tw in time_windows) {
    test_data <- data %>%
      filter(time_window == tw) %>%
      inner_join(predrug_data, by = "CIMBI_ID")
    
    if(nrow(test_data) > 0) {
      ttest <- t.test(test_data[[var_name]], test_data$predrug_value, paired = TRUE)
      d <- cohens_d_paired(test_data[[var_name]], test_data$predrug_value)
      
      results_list[[tw]] <- list(ttest = ttest, d = d)
    }
  }
  
  results_df <- data.frame(
    variable = nice_name,
    comparison = paste("predrug vs", names(results_list)),
    t_statistic = sapply(results_list, function(x) x$ttest$statistic),
    p_value = sapply(results_list, function(x) x$ttest$p.value),
    cohens_d = sapply(results_list, function(x) x$d),
    df = sapply(results_list, function(x) x$ttest$parameter)
  ) %>%
    mutate(across(where(is.numeric), \(x) round(x, 4)))
  
  return(list(plot = plot, results = results_df))
}

run_mixed_model <- function(data, dependent_var, predictor_var, var_names) {
  
  cat(sprintf("    Preparing data for %s...\n", dependent_var))
  
  model_data <- data %>%
    select(all_of(c(dependent_var, predictor_var, "FD", "Age", "Sex", "scanner", 
                    "CIMBI_ID", "scan", "time_since_admin"))) %>%
    na.omit()
  
  cat(sprintf("    Data prepared: %d observations, %d subjects\n", 
              nrow(model_data), length(unique(model_data$CIMBI_ID))))
  
  model_data$CIMBI_ID <- as.factor(model_data$CIMBI_ID)
  model_data$scan <- as.factor(model_data$scan)
  
  formula_strMMRM <- paste(dependent_var, "~", predictor_var, "+ FD + Age + Sex + scanner + us(scan | CIMBI_ID)")
  cat(sprintf("    Fitting main model...\n"))
  
  model <- mmrm(as.formula(formula_strMMRM), data = model_data)
  
  if(model$opt_details$convergence != 0) {
    stop(sprintf("Model did not converge for %s. Convergence code: %d", 
                 dependent_var, model$opt_details$convergence))
  }
  
  cat(sprintf("    Main model fitted successfully.\n"))
  
  coef_summary <- summary(model)$coefficients
  pred_row <- which(rownames(coef_summary) == predictor_var)
  
  if(length(pred_row) == 0) {
    stop(sprintf("Could not find predictor '%s' in model coefficients", predictor_var))
  }
  
  results <- data.frame(
    Measure = dependent_var,
    Predictor = predictor_var,
    Estimate = coef_summary[pred_row, "Estimate"],
    SE = coef_summary[pred_row, "Std. Error"],
    DF = coef_summary[pred_row, "df"],
    tStat = coef_summary[pred_row, "t value"],
    pValue = coef_summary[pred_row, "Pr(>|t|)"]
  )
  
  # Calculate partial correlation using Satterthwaite df
  # rho = sign(beta) * sqrt(t^2 / (t^2 + df))
  t_val <- results$tStat
  df_val <- results$DF
  partial_r <- sign(results$Estimate) * sqrt(t_val^2 / (t_val^2 + df_val))
  
  results$partial_r <- partial_r
  
  model_data$partial_residuals <- residuals(model) + coef(model)[predictor_var] * model_data[[predictor_var]]
  
  plot_data <- data.frame(
    predictor = model_data[[predictor_var]],
    partial_residuals = model_data$partial_residuals,
    time_since_admin = model_data$time_since_admin
  )
  
  p_val <- results$pValue
  if(p_val < 0.001) {
    p_text <- sprintf("p < 0.001")
  } else {
    p_text <- sprintf("p = %.3f", p_val)
  }
  
  x_label <- case_when(
    predictor_var == "psi_conc" ~ "Plasma Psilocin Level (mg/L)",
    predictor_var == "SDI_score" ~ "SDI Score",
    predictor_var == "occupancy" ~ "5-HT2AR Occupancy (%)",
    predictor_var == "occupancy_ofmax" ~ "5-HT2AR Occupancy (% of maximum)",
    TRUE ~ predictor_var
  )
  
  use_scientific <- (dependent_var == "insideout")
  
  p <- ggplot(plot_data, aes(x = predictor, y = partial_residuals)) +
    geom_point(aes(colour = time_since_admin), alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    xlab(x_label) +
    ylab(paste("Partial Residuals of", var_names[[dependent_var]], "\n(Adjusted for age, sex, scanner and motion)")) +
    annotate("text", x = Inf, y = Inf, label = p_text, hjust = 1.1, vjust = 1.1) +
    theme_bw() +
    scale_colour_gradientn(
      colors = c("black", "purple", "darkorange3"),
      values = scales::rescale(c(0, 160, max(model_data$time_since_admin, na.rm = TRUE))),
      breaks = c(0, 120, 240, 360),
      labels = c("Predrug", "+120 mins", "+240 mins", "+360 mins")
    ) +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 12),
      legend.text = element_text(size = 16),
      plot.title = element_blank(),
      legend.position = "none"
    )
  
  if(use_scientific) {
    p <- p + scale_y_continuous(labels = function(x) {
      ifelse(x == 0, "0", {
        exponent <- floor(log10(abs(x)))
        coefficient <- x / 10^exponent
        parse(text = paste0(round(coefficient, 1), " %*% 10^", exponent))
      })
    })
  }
  
  return(list(results = results, plot = p, model = model, model_data = model_data))
}

# -----------------------------------------------------------------------------
# 4. REGIONAL ANALYSIS FUNCTIONS
# -----------------------------------------------------------------------------

run_regional_analysis <- function() {
  cat("\n=== RUNNING REGIONAL ANALYSIS WITH UNSTRUCTURED COVARIANCE ===\n")
  
  roi_table <- read.csv(file.path(roi_dir, 'schaefer2018/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv'))
  
  network_labels <- sapply(1:200, function(i) {
    roi_name <- roi_table$ROI.Name[i]
    parts <- strsplit(roi_name, "_")[[1]]
    if(length(parts) >= 4) {
      return(parts[3])
    } else {
      return("Unknown")
    }
  })
  
  network_labels <- c(network_labels, rep("SubC", 32))
  
  combined_metadata <- combined_data %>%
    select(CIMBI_ID, time_since_admin, scan, Age, Sex, scanner, FD) %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID))
  
  measures <- list(
    directionality = list(
      file = "subject_directionality_scores.csv",
      name = "Directionality"
    ),
    regional_GCOR = list(
      file = "subject_regional_GCOR_scores.csv",
      name = "Regional_GCOR"
    ),
    forward_MI = list(
      file = "subject_forward_MI_scores.csv",
      name = "Forward_MI"
    ),
    reverse_MI = list(
      file = "subject_reverse_MI_scores.csv",
      name = "Reverse_MI"
    )
  )
  
  all_results <- list()
  
  for(measure_key in names(measures)) {
    measure_info <- measures[[measure_key]]
    cat(sprintf("\nProcessing %s...\n", measure_info$name))
    
    regional_data <- read.csv(file.path(out_dir, measure_info$file))
    
    regional_long <- regional_data %>%
      pivot_longer(cols = starts_with("Region_"),
                   names_to = "Region",
                   values_to = "Score") %>%
      mutate(
        Region = as.numeric(gsub("Region_", "", Region)),
        CIMBI_ID = as.character(CIMBI_ID)
      ) %>%
      left_join(combined_metadata, by = c("CIMBI_ID", "time_since_admin")) %>%
      filter(!is.na(Score) & !is.na(FD) & !is.na(scan))
    
    n_regions <- 232
    estimates <- numeric(n_regions)
    ses <- numeric(n_regions)
    pvalues <- numeric(n_regions)
    dfs <- numeric(n_regions)
    tstats <- numeric(n_regions)
    partial_rs <- numeric(n_regions)
    
    for(reg in 1:n_regions) {
      if(reg %% 20 == 0) {
        cat(sprintf("  Processing region %d/%d\n", reg, n_regions))
      }
      
      region_data <- regional_long %>%
        filter(Region == reg) %>%
        na.omit()
      
      region_data$CIMBI_ID <- as.factor(region_data$CIMBI_ID)
      region_data$scan <- as.factor(region_data$scan)
      
      tryCatch({
        formula_strMMRM <- "Score ~ psi_conc + FD + Age + Sex + scanner + us(scan | CIMBI_ID)"
        model <- mmrm(as.formula(formula_strMMRM), data = region_data)
        
        if(model$opt_details$convergence != 0) {
          warning(sprintf("Region %d did not converge. Convergence code: %d", 
                          reg, model$opt_details$convergence))
          estimates[reg] <- NA
          ses[reg] <- NA
          pvalues[reg] <- NA
          dfs[reg] <- NA
          tstats[reg] <- NA
          partial_rs[reg] <- NA
          next
        }
        
        coef_table <- summary(model)$coefficients
        psi_row <- which(rownames(coef_table) == "psi_conc")
        
        if(length(psi_row) > 0) {
          estimates[reg] <- coef_table[psi_row, "Estimate"]
          ses[reg] <- coef_table[psi_row, "Std. Error"]
          pvalues[reg] <- coef_table[psi_row, "Pr(>|t|)"]
          dfs[reg] <- coef_table[psi_row, "df"]
          tstats[reg] <- coef_table[psi_row, "t value"]
          
          # Calculate partial correlation using Satterthwaite df
          t_val  <- tstats[reg]
          df_val <- dfs[reg]
          partial_rs[reg] <- sign(estimates[reg]) * sqrt(t_val^2 / (t_val^2 + df_val))
        } else {
          warning(sprintf("Could not find psi_conc coefficient for region %d", reg))
          estimates[reg] <- NA
          ses[reg] <- NA
          pvalues[reg] <- NA
          dfs[reg] <- NA
          tstats[reg] <- NA
          partial_rs[reg] <- NA
        }
        
      }, error = function(e) {
        warning(sprintf("Failed to fit model for region %d: %s", reg, e$message))
        estimates[reg] <- NA
        ses[reg] <- NA
        pvalues[reg] <- NA
        dfs[reg] <- NA
        tstats[reg] <- NA
        partial_rs[reg] <- NA
      })
    }
    
    bonferroni_threshold <- 0.05 / n_regions
    bonferroni_significant <- ifelse(pvalues < bonferroni_threshold, 1, 0)
    
    results_df <- data.frame(
      Region = 1:n_regions,
      Network = network_labels,
      Measure = measure_info$name,
      Estimate = estimates,
      SE = ses,
      PValue = pvalues,
      DF = dfs,
      tStat = tstats,
      partial_r = partial_rs,
      Holm_Significant = bonferroni_significant
    ) %>%
      mutate(
        CI_lower = Estimate - 1.96 * SE,
        CI_upper = Estimate + 1.96 * SE
      )
    
    all_results[[measure_key]] <- results_df
    
    n_uncorr <- sum(pvalues < 0.05, na.rm = TRUE)
    n_bonf <- sum(bonferroni_significant, na.rm = TRUE)
    cat(sprintf("  %s: %d uncorrected, %d Bonferroni-corrected\n", 
                measure_info$name, n_uncorr, n_bonf))
  }
  
  combined_results <- bind_rows(all_results)
  
  write.csv(combined_results, 
            file.path(out_dir, "all_regional_results_for_R.csv"),
            row.names = FALSE)
  
  cat(sprintf("\nRegional analysis complete. Results saved to:\n%s\n", 
              file.path(out_dir, "all_regional_results_for_R.csv")))
  
  return(combined_results)
}

# -----------------------------------------------------------------------------
# 5. CORRELATION HEATMAP FUNCTION
# -----------------------------------------------------------------------------

create_correlation_heatmap <- function(data, fig_dir, out_dir) {
  measure_order <- c("global_GCOR", "insideout", "modularity_fc_abs", "modularity_forward_MI_abs")

  nice_labels <- c(
    "global_GCOR"               = "Global GCOR",
    "insideout"                 = "INSIDEOUT",
    "modularity_fc_abs"         = "Modularity (FC)",
    "modularity_forward_MI_abs" = "Modularity (Fwd MI)"
  )

  nice_tw <- c(
    "predrug"   = "Pre-drug",
    "upslope"   = "Upslope",
    "peak"      = "Peak",
    "downslope" = "Downslope",
    "return"    = "Return"
  )

  time_windows <- c("predrug", "upslope", "peak", "downslope", "return")

  # Helper: build heatmap from a correlation matrix
  # ci_df (optional): data.frame with columns measure_1, measure_2, ci_lower, ci_upper
  make_heatmap <- function(cor_matrix, title_text, ci_df = NULL) {
    cor_df <- as.data.frame(cor_matrix) %>%
      mutate(var1 = rownames(.)) %>%
      pivot_longer(cols = -var1, names_to = "var2", values_to = "correlation") %>%
      mutate(
        var1 = factor(var1, levels = rev(measure_order)),
        var2 = factor(var2, levels = measure_order)
      )

    if (!is.null(ci_df)) {
      # Allow CI lookup in both directions (upper-triangle pairs only stored once)
      ci_both <- bind_rows(
        ci_df %>% rename(var1 = measure_1, var2 = measure_2),
        ci_df %>% rename(var1 = measure_2, var2 = measure_1)
      )
      cor_df <- cor_df %>%
        left_join(ci_both %>% select(var1, var2, ci_lower, ci_upper),
                  by = c("var1", "var2")) %>%
        mutate(cell_label = ifelse(
          !is.na(ci_lower),
          sprintf("%.2f\n[%.2f, %.2f]", correlation, ci_lower, ci_upper),
          sprintf("%.2f", correlation)
        ))
      txt_size <- 6
    } else {
      cor_df    <- cor_df %>% mutate(cell_label = sprintf("%.2f", correlation))
      txt_size  <- 10
    }

    ggplot(cor_df, aes(x = var2, y = var1, fill = correlation)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = cell_label), color = "black", size = txt_size,
                fontface = "bold", lineheight = 0.85) +
      scale_fill_gradient2(
        low = "#2166ac", mid = "ivory", high = "#b2182b",
        midpoint = 0, limits = c(-1, 1), name = "Correlation"
      ) +
      scale_x_discrete(labels = nice_labels, position = "bottom") +
      scale_y_discrete(labels = nice_labels) +
      coord_fixed() +
      theme_minimal() +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 17),
        axis.text.y  = element_text(size = 17),
        axis.title   = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text  = element_text(size = 13),
        panel.grid   = element_blank(),
        plot.title   = element_text(hjust = 0.5, size = 16, face = "bold")
      ) +
      labs(title = title_text)
  }

  # --- Per-timepoint correlations ---
  cor_matrices    <- list()
  all_cor_results <- list()
  plots           <- list()

  for(tw in time_windows) {
    tw_data <- data %>%
      filter(time_window == tw) %>%
      select(all_of(measure_order)) %>%
      na.omit()

    cat(sprintf("  %s: %d observations\n", tw, nrow(tw_data)))
    if(nrow(tw_data) < 3) {
      cat(sprintf("  Skipping %s — too few observations\n", tw))
      next
    }

    cor_mat            <- cor(tw_data)
    cor_matrices[[tw]] <- cor_mat

    # Pairwise cor.test for this timepoint
    tw_results <- data.frame()
    for(i in 1:(length(measure_order) - 1)) {
      for(j in (i + 1):length(measure_order)) {
        var1 <- measure_order[i]
        var2 <- measure_order[j]
        test <- cor.test(tw_data[[var1]], tw_data[[var2]])
        tw_results <- rbind(tw_results, data.frame(
          time_window = tw,
          measure_1   = var1,
          measure_2   = var2,
          correlation = test$estimate,
          ci_lower    = test$conf.int[1],
          ci_upper    = test$conf.int[2],
          t_statistic = test$statistic,
          df          = test$parameter,
          p_value     = test$p.value,
          n           = nrow(tw_data)
        ))
      }
    }
    all_cor_results[[tw]] <- tw_results
  }

  # --- Mean correlation matrix via Fisher z-transform ---
  z_matrices      <- lapply(cor_matrices, atanh)
  mean_z_matrix   <- Reduce("+", z_matrices) / length(z_matrices)
  mean_cor_matrix <- tanh(mean_z_matrix)
  diag(mean_cor_matrix) <- 1  # numerical clean-up

  # --- Bootstrap CIs (5000 bootstraps, percentile method) ---
  # Resample subjects with replacement; keep all timepoints for each drawn subject
  cat("  Computing bootstrap CIs (5000 iterations)...\n")
  n_boot_cor   <- 5000
  boot_subjects_cor <- unique(data$CIMBI_ID)
  n_subj_cor   <- length(boot_subjects_cor)
  pairs        <- combn(measure_order, 2, simplify = FALSE)
  pair_names   <- sapply(pairs, function(p) paste(p, collapse = "."))

  # Slim data for speed — only columns we need
  data_slim <- data %>% select(CIMBI_ID, time_window, all_of(measure_order))

  # Storage
  boot_tw_cors   <- setNames(
    lapply(pair_names, function(x) matrix(NA_real_, nrow = n_boot_cor,
                                           ncol  = length(time_windows),
                                           dimnames = list(NULL, time_windows))),
    pair_names
  )
  boot_mean_cors <- matrix(NA_real_, nrow = n_boot_cor, ncol = length(pair_names),
                           dimnames = list(NULL, pair_names))

  set.seed(123)
  for(b in 1:n_boot_cor) {
    drawn <- sample(boot_subjects_cor, n_subj_cor, replace = TRUE)
    boot_b <- bind_rows(lapply(seq_along(drawn), function(i) {
      data_slim %>% filter(CIMBI_ID == drawn[i])
    }))

    for(pi in seq_along(pairs)) {
      v1    <- pairs[[pi]][1]
      v2    <- pairs[[pi]][2]
      pname <- pair_names[pi]
      z_acc <- numeric(0)
      for(tw in time_windows) {
        d_tw <- boot_b %>%
          filter(time_window == tw) %>%
          select(all_of(c(v1, v2))) %>%
          na.omit()
        if(nrow(d_tw) >= 3) {
          r <- cor(d_tw[[v1]], d_tw[[v2]])
          boot_tw_cors[[pname]][b, tw] <- r
          z_acc <- c(z_acc, atanh(r))
        }
      }
      if(length(z_acc) > 0) boot_mean_cors[b, pname] <- tanh(mean(z_acc))
    }
  }

  # Percentile CIs — per timepoint
  boot_ci_tw <- bind_rows(lapply(pair_names, function(pname) {
    parts <- strsplit(pname, "\\.")[[1]]
    bind_rows(lapply(time_windows, function(tw) {
      vals <- boot_tw_cors[[pname]][, tw]
      data.frame(time_window = tw, measure_1 = parts[1], measure_2 = parts[2],
                 ci_lower = quantile(vals, 0.025, na.rm = TRUE),
                 ci_upper = quantile(vals, 0.975, na.rm = TRUE))
    }))
  }))

  # Percentile CIs — mean across timepoints
  boot_ci_mean <- bind_rows(lapply(pair_names, function(pname) {
    parts <- strsplit(pname, "\\.")[[1]]
    vals  <- boot_mean_cors[, pname]
    data.frame(measure_1 = parts[1], measure_2 = parts[2],
               ci_lower = quantile(vals, 0.025, na.rm = TRUE),
               ci_upper = quantile(vals, 0.975, na.rm = TRUE))
  }))

  # --- Generate heatmaps with bootstrap CIs ---
  for(tw in time_windows) {
    if(is.null(cor_matrices[[tw]])) next
    ci_tw       <- boot_ci_tw %>%
      filter(time_window == tw) %>%
      select(measure_1, measure_2, ci_lower, ci_upper)
    p_tw        <- make_heatmap(cor_matrices[[tw]],
                                sprintf("Intercorrelations — %s", nice_tw[tw]),
                                ci_df = ci_tw)
    plots[[tw]] <- p_tw
    ggsave(file.path(fig_dir, sprintf("correlation_heatmap_%s.png", tw)),
           p_tw, width = 8, height = 7, dpi = 300)
  }

  # Save intermediate data so the mean heatmap can be regenerated without rerunning bootstraps
  saveRDS(mean_cor_matrix, file.path(out_dir, "correlation_mean_cor_matrix.rds"))
  write_csv(boot_ci_mean,  file.path(out_dir, "correlation_boot_ci_mean.csv"))

  p_mean          <- make_heatmap(mean_cor_matrix,
                                  "Intercorrelation \u2014 Mean Across Timepoints",
                                  ci_df = boot_ci_mean)
  plots[["mean"]] <- p_mean
  ggsave(file.path(fig_dir, "correlation_heatmap_mean.png"),
         p_mean, width = 8, height = 7, dpi = 300)

  # --- Save outputs ---
  all_cor_results_df <- bind_rows(all_cor_results) %>%
    select(-ci_lower, -ci_upper) %>%
    left_join(boot_ci_tw, by = c("time_window", "measure_1", "measure_2"))
  write_csv(all_cor_results_df,
            file.path(out_dir, "global_measures_correlations_by_timepoint.csv"))

  # Mean per pair (same column structure as before for supplementary table compatibility)
  mean_cor_results <- data.frame()
  for(i in 1:(length(measure_order) - 1)) {
    for(j in (i + 1):length(measure_order)) {
      var1   <- measure_order[i]
      var2   <- measure_order[j]
      z_vals <- sapply(z_matrices, function(m) m[var1, var2])
      ci_row <- boot_ci_mean %>% filter(measure_1 == var1, measure_2 == var2)
      mean_cor_results <- rbind(mean_cor_results, data.frame(
        measure_1    = var1,
        measure_2    = var2,
        correlation  = tanh(mean(z_vals)),
        ci_lower     = if(nrow(ci_row) == 1) ci_row$ci_lower else NA_real_,
        ci_upper     = if(nrow(ci_row) == 1) ci_row$ci_upper else NA_real_,
        t_statistic  = NA_real_,
        df           = NA_real_,
        p_value      = NA_real_,
        n_timepoints = length(z_vals)
      ))
    }
  }
  write_csv(mean_cor_results, file.path(out_dir, "global_measures_correlations.csv"))

  mean_cor_df <- as.data.frame(mean_cor_matrix) %>%
    mutate(measure = rownames(.)) %>%
    select(measure, everything())
  write_csv(mean_cor_df, file.path(out_dir, "global_measures_correlation_matrix.csv"))

  # --- Line plot: pairwise correlations across timepoints ---
  line_data <- all_cor_results_df %>%
    mutate(
      label       = paste(measure_1, "\u00d7", measure_2),
      time_window = factor(time_window, levels = time_windows)
    )

  p_lines <- ggplot(line_data, aes(x = time_window, y = correlation,
                                    color = label, fill = label, group = label)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_discrete(labels = nice_tw) +
    scale_color_brewer(palette = "Dark2", name = "Metric pair") +
    scale_fill_brewer(palette = "Dark2", guide = "none") +
    labs(x = "Timepoint", y = "Pearson r",
         title = "Pairwise Metric Correlations Across Timepoints") +
    theme_bw() +
    theme(
      axis.text        = element_text(size = 14),
      axis.title       = element_text(size = 16),
      legend.text      = element_text(size = 11),
      legend.title     = element_text(size = 12),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  plots[["lines"]] <- p_lines
  ggsave(file.path(fig_dir, "correlation_by_timepoint_lines.png"),
         p_lines, width = 9, height = 6, dpi = 300)

  cat("Correlation heatmaps saved (per timepoint + mean across timepoints)\n")

  return(list(
    plots                = plots,
    matrix               = mean_cor_matrix,
    results              = mean_cor_results,
    results_by_timepoint = all_cor_results_df
  ))
}

# -----------------------------------------------------------------------------
# 6. NULL DISTRIBUTION VS GCOR PLOTS
# -----------------------------------------------------------------------------

# This script contains the updated plot_null_distributions_vs_gcor function
# with modifications to the actual modularity plots as requested


plot_null_distributions_vs_gcor <- function(data, out_dir, fig_dir) {
  cat("\n=== PLOTTING NULL DISTRIBUTIONS VS GCOR ===\n")
  
  # Load null distribution files
  null_fc <- read_csv(file.path(out_dir, "null_distributions_fc.csv"))
  null_fmi <- read_csv(file.path(out_dir, "null_distributions_fmi.csv"))
  
  # Get GCOR values for each subject
  gcor_data <- data %>%
    select(CIMBI_ID, time_since_admin, global_GCOR, insideout) %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID))
  
  # Merge GCOR and INSIDEOUT with null distributions
  null_fc <- null_fc %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
    left_join(gcor_data, by = c("CIMBI_ID", "time_since_admin"))
  
  null_fmi <- null_fmi %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
    left_join(gcor_data, by = c("CIMBI_ID", "time_since_admin"))
  
  # Shared colour scale parameters (used across all null/absolute plots for consistency)
  time_max_fc  <- max(null_fc$time_since_admin,  na.rm = TRUE)
  time_max_fmi <- max(null_fmi$time_since_admin, na.rm = TRUE)

  time_colours  <- c("black", "purple", "darkorange3")
  time_values_fc  <- scales::rescale(c(0, 160, time_max_fc))
  time_values_fmi <- scales::rescale(c(0, 160, time_max_fmi))

  shared_theme <- theme_bw() +
    theme(
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 12),
      legend.position  = "right"
    )

  # Plot 1: FC - ALL null modularity values vs GCOR
  p_fc_null <- ggplot(null_fc %>% filter(distribution_type == "null"),
                      aes(x = global_GCOR, y = modularity)) +
    geom_point(aes(color = time_since_admin), alpha = 0.3, size = 0.8) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fc, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (FC)", subtitle = "Null") +
    shared_theme

  # Plot 2a: FC - All absolute modularity iterations vs GCOR
  p_fc_absolute_all <- ggplot(null_fc %>% filter(distribution_type == "actual"),
                              aes(x = global_GCOR, y = modularity)) +
    geom_point(aes(color = time_since_admin), alpha = 0.3, size = 0.8) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fc, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (FC)", subtitle = "Absolute") +
    shared_theme

  # Plot 2b: FC - Mean absolute modularity vs GCOR (as used in regression)
  fc_mean <- null_fc %>%
    filter(distribution_type == "actual") %>%
    group_by(CIMBI_ID, time_since_admin, global_GCOR) %>%
    summarise(mean_modularity = mean(modularity, na.rm = TRUE), .groups = 'drop')

  p_fc_absolute_mean <- ggplot(fc_mean,
                               aes(x = global_GCOR, y = mean_modularity)) +
    geom_point(aes(color = time_since_admin), alpha = 0.7, size = 2) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fc, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (FC)", subtitle = "Absolute (mean)") +
    shared_theme

  # Plot 2c: FC - Combined null (grey) + absolute mean (time colour) on one panel
  # Uses per-scan mean so each scan contributes one absolute point against the null cloud
  p_fc_combined <- ggplot() +
    geom_point(data = null_fc %>% filter(distribution_type == "null"),
               aes(x = global_GCOR, y = modularity),
               color = "grey70", alpha = 0.15, size = 0.5) +
    geom_point(data = fc_mean,
               aes(x = global_GCOR, y = mean_modularity, color = time_since_admin),
               alpha = 0.7, size = 1.5) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fc, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (FC)", subtitle = "Null (grey) vs Absolute") +
    shared_theme

  # Plot 3: Forward MI - ALL null modularity values vs GCOR
  p_fmi_null <- ggplot(null_fmi %>% filter(distribution_type == "null"),
                       aes(x = global_GCOR, y = modularity)) +
    geom_point(aes(color = time_since_admin), alpha = 0.3, size = 0.8) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fmi, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (Fwd MI)", subtitle = "Null") +
    shared_theme

  # Plot 4a: Forward MI - All absolute modularity iterations vs GCOR
  p_fmi_absolute_all <- ggplot(null_fmi %>% filter(distribution_type == "actual"),
                               aes(x = global_GCOR, y = modularity)) +
    geom_point(aes(color = time_since_admin), alpha = 0.3, size = 0.8) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fmi, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (Fwd MI)", subtitle = "Absolute") +
    shared_theme

  # Plot 4b: Forward MI - Mean absolute modularity vs GCOR (as used in regression)
  fmi_mean <- null_fmi %>%
    filter(distribution_type == "actual") %>%
    group_by(CIMBI_ID, time_since_admin, global_GCOR, insideout) %>%
    summarise(mean_modularity = mean(modularity, na.rm = TRUE), .groups = 'drop')

  p_fmi_absolute_mean <- ggplot(fmi_mean,
                                aes(x = global_GCOR, y = mean_modularity)) +
    geom_point(aes(color = time_since_admin), alpha = 0.7, size = 2) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fmi, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (Fwd MI)", subtitle = "Absolute (mean)") +
    shared_theme

  # Plot 4c: Forward MI - Combined null (grey) + absolute mean (time colour) on one panel
  p_fmi_combined <- ggplot() +
    geom_point(data = null_fmi %>% filter(distribution_type == "null"),
               aes(x = global_GCOR, y = modularity),
               color = "grey70", alpha = 0.15, size = 0.5) +
    geom_point(data = fmi_mean,
               aes(x = global_GCOR, y = mean_modularity, color = time_since_admin),
               alpha = 0.7, size = 1.5) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fmi, name = "Time (min)") +
    labs(x = "Global GCOR", y = "Modularity (Fwd MI)", subtitle = "Null (grey) vs Absolute") +
    shared_theme

  # Plot 5: Forward MI - Mean absolute modularity vs INSIDEOUT
  p_fmi_vs_insideout <- ggplot(fmi_mean,
                               aes(x = insideout, y = mean_modularity)) +
    geom_point(aes(color = time_since_admin), alpha = 0.7, size = 2) +
    scale_colour_gradientn(colors = time_colours, values = time_values_fmi, name = "Time (min)") +
    scale_x_continuous(labels = function(x) {
      ifelse(x == 0, "0", {
        exponent <- floor(log10(abs(x)))
        coefficient <- x / 10^exponent
        parse(text = paste0(round(coefficient, 1), " %*% 10^", exponent))
      })
    }) +
    labs(x = "INSIDEOUT", y = "Modularity (Fwd MI)") +
    shared_theme

  # Save plots — all individual panels at 4.5x4 for easy side-by-side assembly
  ggsave(file.path(fig_dir, "null_modularity_vs_gcor_fc_all_runs.png"),
         p_fc_null, width = 4.5, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "null_modularity_vs_gcor_fmi_all_runs.png"),
         p_fmi_null, width = 4.5, height = 4, dpi = 300)

  ggsave(file.path(fig_dir, "absolute_modularity_vs_gcor_fc_all_iterations.png"),
         p_fc_absolute_all, width = 4.5, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "absolute_modularity_vs_gcor_fmi_all_iterations.png"),
         p_fmi_absolute_all, width = 4.5, height = 4, dpi = 300)

  ggsave(file.path(fig_dir, "absolute_modularity_vs_gcor_fc_mean.png"),
         p_fc_absolute_mean, width = 4.5, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "absolute_modularity_vs_gcor_fmi_mean.png"),
         p_fmi_absolute_mean, width = 4.5, height = 4, dpi = 300)

  # Combined null + absolute on one panel
  ggsave(file.path(fig_dir, "combined_null_absolute_modularity_vs_gcor_fc.png"),
         p_fc_combined, width = 4.5, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "combined_null_absolute_modularity_vs_gcor_fmi.png"),
         p_fmi_combined, width = 4.5, height = 4, dpi = 300)

  # Modularity vs INSIDEOUT plot
  ggsave(file.path(fig_dir, "absolute_modularity_vs_insideout_fmi_mean.png"),
         p_fmi_vs_insideout, width = 4.5, height = 4, dpi = 300)

  # Print summary statistics
  cat(sprintf("\nFC absolute - All iterations: %d points\n", nrow(null_fc %>% filter(distribution_type == "actual"))))
  cat(sprintf("FC absolute - Mean values: %d scans\n", nrow(fc_mean)))
  cat(sprintf("FMI absolute - All iterations: %d points\n", nrow(null_fmi %>% filter(distribution_type == "actual"))))
  cat(sprintf("FMI absolute - Mean values: %d scans\n", nrow(fmi_mean)))
  
  # Calculate all correlations
  cat("\n=== CALCULATING CORRELATIONS ===\n")
  
  cor_fc_null <- cor.test(null_fc %>% filter(distribution_type == "null") %>% pull(modularity),
                          null_fc %>% filter(distribution_type == "null") %>% pull(global_GCOR))
  cor_fc_absolute_all <- cor.test(null_fc %>% filter(distribution_type == "actual") %>% pull(modularity),
                                  null_fc %>% filter(distribution_type == "actual") %>% pull(global_GCOR))
  cor_fc_absolute_mean <- cor.test(fc_mean$mean_modularity, fc_mean$global_GCOR)

  cor_fmi_null <- cor.test(null_fmi %>% filter(distribution_type == "null") %>% pull(modularity),
                           null_fmi %>% filter(distribution_type == "null") %>% pull(global_GCOR))
  cor_fmi_absolute_all <- cor.test(null_fmi %>% filter(distribution_type == "actual") %>% pull(modularity),
                                   null_fmi %>% filter(distribution_type == "actual") %>% pull(global_GCOR))
  cor_fmi_absolute_mean <- cor.test(fmi_mean$mean_modularity, fmi_mean$global_GCOR)
  
  # Correlation between Modularity (Fwd MI) and INSIDEOUT
  fmi_insideout_clean <- fmi_mean %>% filter(!is.na(insideout) & !is.na(mean_modularity))
  cor_fmi_vs_insideout <- cor.test(fmi_insideout_clean$mean_modularity, fmi_insideout_clean$insideout)
  
  cat(sprintf("\nModularity (Fwd MI) vs INSIDEOUT:\n"))
  cat(sprintf("  Correlation: r = %.4f, 95%% CI [%.4f, %.4f]\n", 
              cor_fmi_vs_insideout$estimate,
              cor_fmi_vs_insideout$conf.int[1],
              cor_fmi_vs_insideout$conf.int[2]))
  cat(sprintf("  t(%d) = %.2f, p = %.4e\n",
              cor_fmi_vs_insideout$parameter,
              cor_fmi_vs_insideout$statistic,
              cor_fmi_vs_insideout$p.value))
  cat(sprintf("  N = %d scans\n\n", nrow(fmi_insideout_clean)))
      
      # Save correlation results
      cor_results <- data.frame(
        network = c("FC", "FC", "FC", "Forward MI", "Forward MI", "Forward MI", "Forward MI"),
        distribution_type = c("null", "absolute_all_iterations", "absolute_mean",
                              "null", "absolute_all_iterations", "absolute_mean", "absolute_mean"),
        comparison = c("vs_GCOR", "vs_GCOR", "vs_GCOR", "vs_GCOR", "vs_GCOR", "vs_GCOR", "vs_INSIDEOUT"),
        correlation = c(cor_fc_null$estimate, cor_fc_absolute_all$estimate, cor_fc_absolute_mean$estimate,
                        cor_fmi_null$estimate, cor_fmi_absolute_all$estimate, cor_fmi_absolute_mean$estimate,
                        cor_fmi_vs_insideout$estimate),
        p_value = c(cor_fc_null$p.value, cor_fc_absolute_all$p.value, cor_fc_absolute_mean$p.value,
                    cor_fmi_null$p.value, cor_fmi_absolute_all$p.value, cor_fmi_absolute_mean$p.value,
                    cor_fmi_vs_insideout$p.value),
        ci_lower = c(cor_fc_null$conf.int[1], cor_fc_absolute_all$conf.int[1], cor_fc_absolute_mean$conf.int[1],
                     cor_fmi_null$conf.int[1], cor_fmi_absolute_all$conf.int[1], cor_fmi_absolute_mean$conf.int[1],
                     cor_fmi_vs_insideout$conf.int[1]),
        ci_upper = c(cor_fc_null$conf.int[2], cor_fc_absolute_all$conf.int[2], cor_fc_absolute_mean$conf.int[2],
                     cor_fmi_null$conf.int[2], cor_fmi_absolute_all$conf.int[2], cor_fmi_absolute_mean$conf.int[2],
                     cor_fmi_vs_insideout$conf.int[2]),
        n_points = c(nrow(null_fc %>% filter(distribution_type == "null")),
                     nrow(null_fc %>% filter(distribution_type == "actual")),
                     nrow(fc_mean),
                     nrow(null_fmi %>% filter(distribution_type == "null")),
                     nrow(null_fmi %>% filter(distribution_type == "actual")),
                     nrow(fmi_mean),
                     nrow(fmi_insideout_clean))
      )
      
      write_csv(cor_results, file.path(out_dir, "modularity_gcor_correlations.csv"))
      
      cat("Null distribution plots saved (showing all individual runs)\n")
      cat("\nCorrelation summary:\n")
      print(cor_results)
      
      return(list(
        plots = list(fc_null             = p_fc_null,
                     fc_absolute_all     = p_fc_absolute_all,
                     fc_absolute_mean    = p_fc_absolute_mean,
                     fc_combined         = p_fc_combined,
                     fmi_null            = p_fmi_null,
                     fmi_absolute_all    = p_fmi_absolute_all,
                     fmi_absolute_mean   = p_fmi_absolute_mean,
                     fmi_combined        = p_fmi_combined,
                     fmi_vs_insideout    = p_fmi_vs_insideout),
        correlations = cor_results
      ))
      }

# -----------------------------------------------------------------------------
# 7. REGIONAL FOREST PLOT FUNCTION (FIXED)
# -----------------------------------------------------------------------------

create_regional_plots <- function(all_results, out_dir) {
  cat("\n=== CREATING REGIONAL FOREST PLOTS ===\n")
  
  library(ggplot2)
  library(dplyr)
  
  results_dir <- file.path(out_dir, "figures")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  network_colors <- c(
    "Default" = "#1f77b4",
    "Cont" = "#d62728",
    "DorsAttn" = "#ff7f0e",
    "SalVentAttn" = "#e377c2",
    "Vis" = "#8e9e8b",
    "SomMot" = "#b5b8c9",
    "Limbic" = "#bcbd22",
    "SubC" = "#a0a0a0"
  )
  
  network_order <- c("Default", "Cont", "DorsAttn", "SalVentAttn", "Limbic", "SomMot", "Vis", "SubC")
  
  measures <- unique(all_results$Measure)
  
  for(measure in measures) {
    cat(sprintf("\nCreating forest plot for %s...\n", measure))
    
    plot_data <- all_results %>%
      filter(Measure == measure) %>%
      mutate(
        Network = factor(Network, levels = network_order),
        sig_bonf = PValue < (0.05 / 232),
        sig_uncorr = PValue < 0.05 & !sig_bonf,
        alpha_level = case_when(
          sig_bonf ~ 1.0,
          sig_uncorr ~ 0.6,
          TRUE ~ 0.2
        )
      ) %>%
      arrange(Network, Estimate) %>%
      mutate(
        # Create ordering variable within each network
        network_numeric = as.numeric(Network),
        row_order = row_number()
      )
    
    # Add small gaps between networks
    network_breaks <- plot_data %>%
      group_by(Network) %>%
      summarise(
        max_row = max(row_order),
        min_row = min(row_order),
        .groups = 'drop'
      ) %>%
      mutate(
        gap_at = max_row + 0.5,
        label_pos = (min_row + max_row) / 2
      )
    
    cat(sprintf("  Bonferroni significant (p < %.2e): %d regions\n", 
                0.05/232, sum(plot_data$sig_bonf)))
    cat(sprintf("  Uncorrected significant (p < 0.05): %d regions\n", 
                sum(plot_data$sig_uncorr)))
    cat(sprintf("  Non-significant: %d regions\n", 
                sum(!plot_data$sig_bonf & !plot_data$sig_uncorr)))
    
    p_forest <- ggplot(plot_data, aes(x = Estimate, y = row_order)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      # Add horizontal lines to separate networks
      geom_hline(data = network_breaks[-nrow(network_breaks), ], 
                 aes(yintercept = gap_at), 
                 linetype = "solid", color = "grey80", linewidth = 0.5) +
      geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Network, alpha = alpha_level),
                     height = 0, linewidth = 0.8) +
      geom_point(data = plot_data %>% filter(sig_bonf),
                 aes(color = Network, alpha = alpha_level), 
                 size = 2, shape = 16) +
      geom_point(data = plot_data %>% filter(!sig_bonf),
                 aes(color = Network, alpha = alpha_level), 
                 size = 0.5, shape = 124) +
      scale_color_manual(values = network_colors, name = "Network") +
      scale_alpha_identity() +
      scale_y_continuous(
        breaks = network_breaks$label_pos,
        labels = network_breaks$Network,
        expand = c(0.01, 0.01)
      ) +
      labs(
        x = "Effect Estimate (β) with 95% CI",
        y = "Network",
        title = sprintf("Regional Drug Effects on %s", measure)
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        legend.position = "none",
        plot.title = element_text(size = 14, face = "bold")
      )
    
    filename <- paste0("regional_forest_plot_", gsub(" ", "_", measure), ".png")
    ggsave(file.path(results_dir, filename), 
           p_forest, width = 10, height = 12, dpi = 300)
    cat(sprintf("  Saved: %s\n", filename))
  }
  
  cat("\n=== FOREST PLOTS COMPLETE ===\n")
}

# -----------------------------------------------------------------------------
# 8. HIERARCHY ANALYSIS FUNCTIONS
# -----------------------------------------------------------------------------

calculate_hierarchy_positions <- function() {
  cat("\n=== CALCULATING BASELINE HIERARCHY POSITIONS ===\n")
  
  # Load forward (OUT) and reverse (IN) MI scores
  forward_scores <- read_csv(file.path(out_dir, "subject_forward_MI_scores.csv"), 
                             show_col_types = FALSE)
  reverse_scores <- read_csv(file.path(out_dir, "subject_reverse_MI_scores.csv"), 
                             show_col_types = FALSE)
  
  # Remove scan with missing psi_conc
  reverse_scores <- reverse_scores %>% filter(!is.na(psi_conc) & !is.nan(psi_conc))
  forward_scores <- forward_scores %>% filter(!is.na(psi_conc) & !is.nan(psi_conc))
  
  # Filter to baseline only (psi_conc = 0)
  baseline_forward <- forward_scores %>% filter(psi_conc == 0)
  baseline_reverse <- reverse_scores %>% filter(psi_conc == 0)
  
  # Calculate mean OUT and IN across baseline scans
  region_cols <- grep("^Region_", names(baseline_forward), value = TRUE)
  
  mean_OUT <- colMeans(baseline_forward[, region_cols], na.rm = TRUE)
  mean_IN <- colMeans(baseline_reverse[, region_cols], na.rm = TRUE)
  
  # Calculate hierarchy (OUT - IN)
  hierarchy_positions <- data.frame(
    Region = 1:232,
    mean_OUT = mean_OUT,
    mean_IN = mean_IN,
    mean_hierarchy = mean_OUT - mean_IN
  )
  
  # Load network information
  regional_results <- read_csv(file.path(out_dir, "all_regional_results_for_R.csv"),
                               show_col_types = FALSE)
  
  # Extract unique region-network mapping
  network_mapping <- regional_results %>%
    select(Region, Network) %>%
    distinct()
  
  # Add network to hierarchy_positions
  hierarchy_positions <- hierarchy_positions %>%
    left_join(network_mapping, by = "Region")
  
  # Save
  write_csv(hierarchy_positions, file.path(out_dir, "hierarchy_positions.csv"))
  
  cat(sprintf("Calculated hierarchy positions for %d regions\n", nrow(hierarchy_positions)))
  
  return(list(
    hierarchy_positions = hierarchy_positions,
    forward_scores = forward_scores,
    reverse_scores = reverse_scores
  ))
}

run_hierarchy_analysis <- function(hierarchy_data) {
  cat("\n=== RUNNING HIERARCHY × PSILOCIN INTERACTION ANALYSIS ===\n")
  
  hierarchy_positions <- hierarchy_data$hierarchy_positions
  forward_scores <- hierarchy_data$forward_scores
  reverse_scores <- hierarchy_data$reverse_scores
  
  # Filter data (same as main analysis)
  forward_scores <- forward_scores %>%
    mutate(time_window = case_when(
      time_since_admin < 3 ~ "predrug",
      time_since_admin >= 47 & time_since_admin < 95 ~ "upslope",
      time_since_admin >= 95 & time_since_admin < 140 ~ "peak",
      time_since_admin >= 140 & time_since_admin < 250 ~ "downslope",
      time_since_admin >= 282 ~ "return",
      TRUE ~ "other"
    )) %>%
    filter(time_window != "other") %>%
    group_by(CIMBI_ID, time_window) %>%
    arrange(time_since_admin) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(scan = case_when(
      time_window == "predrug" ~ 0,
      time_window == "upslope" ~ 1,
      time_window == "peak" ~ 2,
      time_window == "downslope" ~ 3,
      time_window == "return" ~ 4
    ))
  
  reverse_scores <- reverse_scores %>%
    mutate(time_window = case_when(
      time_since_admin < 3 ~ "predrug",
      time_since_admin >= 47 & time_since_admin < 95 ~ "upslope",
      time_since_admin >= 95 & time_since_admin < 140 ~ "peak",
      time_since_admin >= 140 & time_since_admin < 250 ~ "downslope",
      time_since_admin >= 282 ~ "return",
      TRUE ~ "other"
    )) %>%
    filter(time_window != "other") %>%
    group_by(CIMBI_ID, time_window) %>%
    arrange(time_since_admin) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(scan = case_when(
      time_window == "predrug" ~ 0,
      time_window == "upslope" ~ 1,
      time_window == "peak" ~ 2,
      time_window == "downslope" ~ 3,
      time_window == "return" ~ 4
    ))
  
  # Load covariates
  copbet_path <- file.path(np2_dir, 'CopBET/')
  copbet_files <- list.files(path = copbet_path, pattern = "^CopBETtbl.*\\.csv$", 
                             full.names = TRUE)
  
  if(length(copbet_files) > 0) {
    fd_data <- read.csv(copbet_files[1]) %>%
      select(CIMBI, time_since_admin, FD) %>%
      rename(CIMBI_ID = CIMBI) %>%
      mutate(CIMBI_ID = as.character(CIMBI_ID))
  } else {
    stop("FD data not found!")
  }
  
  combined_metadata <- combined_data %>%
    select(CIMBI_ID, time_since_admin, Age, Sex, scanner) %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID))
  
  covariates <- combined_metadata %>%
    left_join(fd_data, by = c("CIMBI_ID", "time_since_admin"))
  
  # Join covariates
  reverse_with_covariates <- reverse_scores %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
    left_join(covariates, by = c("CIMBI_ID", "time_since_admin"))
  
  forward_with_covariates <- forward_scores %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
    left_join(covariates, by = c("CIMBI_ID", "time_since_admin"))
  
  # Run regional models for IN
  cat("\nFitting MMRM models for IN (Reverse MI)...\n")
  
  mmrm_data_in <- reverse_with_covariates %>%
    mutate(
      CIMBI_ID = as.factor(CIMBI_ID),
      scan = as.factor(scan),
      Sex = as.factor(Sex),
      scanner = as.factor(scanner)
    )
  
  regional_results_in <- data.frame(
    Region = 1:232,
    beta = numeric(232),
    se = numeric(232),
    t_value = numeric(232),
    df = numeric(232),
    p_value = numeric(232),
    partial_r = numeric(232),
    converged = logical(232)
  )
  
  for(reg in 1:232) {
    if(reg %% 20 == 0) {
      cat(sprintf("  Processing region %d/232...\n", reg))
    }
    
    col_name <- sprintf("Region_%d", reg)
    region_data <- mmrm_data_in %>%
      mutate(IN_score = .data[[col_name]])
    
    tryCatch({
      model <- mmrm(
        IN_score ~ psi_conc + FD + Age + Sex + scanner + us(scan | CIMBI_ID),
        data = region_data
      )
      
      if(model$opt_details$convergence == 0) {
        coef_table <- summary(model)$coefficients
        psi_row <- which(rownames(coef_table) == "psi_conc")
        
        if(length(psi_row) > 0) {
          regional_results_in$beta[reg] <- coef_table[psi_row, "Estimate"]
          regional_results_in$se[reg] <- coef_table[psi_row, "Std. Error"]
          regional_results_in$t_value[reg] <- coef_table[psi_row, "t value"]
          regional_results_in$df[reg] <- coef_table[psi_row, "df"]
          regional_results_in$p_value[reg] <- coef_table[psi_row, "Pr(>|t|)"]
          
          # Calculate partial correlation using Satterthwaite df
          t_val  <- regional_results_in$t_value[reg]
          df_val <- regional_results_in$df[reg]
          regional_results_in$partial_r[reg] <- sign(regional_results_in$beta[reg]) * sqrt(t_val^2 / (t_val^2 + df_val))
          
          regional_results_in$converged[reg] <- TRUE
        }
      }
    }, error = function(e) {
      regional_results_in$converged[reg] <<- FALSE
    })
  }
  
  # Add hierarchy information
  regional_results_in <- regional_results_in %>%
    left_join(hierarchy_positions %>% select(Region, mean_hierarchy, Network), 
              by = "Region")
  
  # Test correlation with hierarchy
  regional_results_in_converged <- regional_results_in %>% filter(converged)
  
  hierarchy_cor_in <- cor.test(
    regional_results_in_converged$mean_hierarchy,
    regional_results_in_converged$beta
  )
  
  cat("\n=== IN (REVERSE MI) RESULTS ===\n")
  cat(sprintf("Correlation: r = %.4f, 95%% CI [%.4f, %.4f]\n", 
              hierarchy_cor_in$estimate, 
              hierarchy_cor_in$conf.int[1], 
              hierarchy_cor_in$conf.int[2]))
  cat(sprintf("t(%d) = %.2f, p = %.4e\n\n", 
              hierarchy_cor_in$parameter, 
              hierarchy_cor_in$statistic, 
              hierarchy_cor_in$p.value))
  
  # Run regional models for OUT
  cat("Fitting MMRM models for OUT (Forward MI)...\n")
  
  mmrm_data_out <- forward_with_covariates %>%
    mutate(
      CIMBI_ID = as.factor(CIMBI_ID),
      scan = as.factor(scan),
      Sex = as.factor(Sex),
      scanner = as.factor(scanner)
    )
  
  regional_results_out <- data.frame(
    Region = 1:232,
    beta = numeric(232),
    se = numeric(232),
    t_value = numeric(232),
    df = numeric(232),
    p_value = numeric(232),
    partial_r = numeric(232),
    converged = logical(232)
  )
  
  for(reg in 1:232) {
    if(reg %% 20 == 0) {
      cat(sprintf("  Processing region %d/232...\n", reg))
    }
    
    col_name <- sprintf("Region_%d", reg)
    region_data <- mmrm_data_out %>%
      mutate(OUT_score = .data[[col_name]])
    
    tryCatch({
      model <- mmrm(
        OUT_score ~ psi_conc + FD + Age + Sex + scanner + us(scan | CIMBI_ID),
        data = region_data
      )
      
      if(model$opt_details$convergence == 0) {
        coef_table <- summary(model)$coefficients
        psi_row <- which(rownames(coef_table) == "psi_conc")
        
        if(length(psi_row) > 0) {
          regional_results_out$beta[reg] <- coef_table[psi_row, "Estimate"]
          regional_results_out$se[reg] <- coef_table[psi_row, "Std. Error"]
          regional_results_out$t_value[reg] <- coef_table[psi_row, "t value"]
          regional_results_out$df[reg] <- coef_table[psi_row, "df"]
          regional_results_out$p_value[reg] <- coef_table[psi_row, "Pr(>|t|)"]
          
          # Calculate partial correlation using Satterthwaite df
          t_val  <- regional_results_out$t_value[reg]
          df_val <- regional_results_out$df[reg]
          regional_results_out$partial_r[reg] <- sign(regional_results_out$beta[reg]) * sqrt(t_val^2 / (t_val^2 + df_val))
          
          regional_results_out$converged[reg] <- TRUE
        }
      }
    }, error = function(e) {
      regional_results_out$converged[reg] <<- FALSE
    })
  }
  
  regional_results_out <- regional_results_out %>%
    left_join(hierarchy_positions %>% select(Region, mean_hierarchy, Network), 
              by = "Region")
  
  regional_results_out_converged <- regional_results_out %>% filter(converged)
  
  hierarchy_cor_out <- cor.test(
    regional_results_out_converged$mean_hierarchy,
    regional_results_out_converged$beta
  )
  
  cat("=== OUT (FORWARD MI) RESULTS ===\n")
  cat(sprintf("Correlation: r = %.4f, 95%% CI [%.4f, %.4f]\n", 
              hierarchy_cor_out$estimate, 
              hierarchy_cor_out$conf.int[1], 
              hierarchy_cor_out$conf.int[2]))
  cat(sprintf("t(%d) = %.2f, p = %.4e\n\n", 
              hierarchy_cor_out$parameter, 
              hierarchy_cor_out$statistic, 
              hierarchy_cor_out$p.value))
  
  # Save results
  write_csv(regional_results_in, file.path(out_dir, "regional_drug_effects_IN_hierarchy.csv"))
  write_csv(regional_results_out, file.path(out_dir, "regional_drug_effects_OUT_hierarchy.csv"))
  
  # Save correlation summary
  hierarchy_correlation_summary <- data.frame(
    measure = c("IN", "OUT"),
    correlation_r = c(hierarchy_cor_in$estimate, hierarchy_cor_out$estimate),
    correlation_ci_lower = c(hierarchy_cor_in$conf.int[1], hierarchy_cor_out$conf.int[1]),
    correlation_ci_upper = c(hierarchy_cor_in$conf.int[2], hierarchy_cor_out$conf.int[2]),
    correlation_t = c(hierarchy_cor_in$statistic, hierarchy_cor_out$statistic),
    correlation_df = c(hierarchy_cor_in$parameter, hierarchy_cor_out$parameter),
    correlation_p = c(hierarchy_cor_in$p.value, hierarchy_cor_out$p.value)
  )
  
  write_csv(hierarchy_correlation_summary, file.path(out_dir, "hierarchy_correlation_summary.csv"))
  
  return(list(
    results_in = regional_results_in,
    results_out = regional_results_out,
    cor_in = hierarchy_cor_in,
    cor_out = hierarchy_cor_out,
    forward_with_covariates = forward_with_covariates,
    reverse_with_covariates = reverse_with_covariates
  ))
}

create_hierarchy_plots <- function(hierarchy_results, hierarchy_positions) {
  cat("\n=== CREATING HIERARCHY PLOTS ===\n")
  
  regional_results_in_converged <- hierarchy_results$results_in %>% filter(converged)
  regional_results_out_converged <- hierarchy_results$results_out %>% filter(converged)
  
  # Plot 1: IN beta vs baseline hierarchy
  p_in_scatter <- ggplot(regional_results_in_converged, 
                         aes(x = mean_hierarchy, y = beta)) +
    geom_point(aes(color = Network), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = network_colors, name = "Network") +
    annotate("segment", x = -0.001, xend = -0.003, 
             y = min(regional_results_in_converged$beta) * 1.1, 
             yend = min(regional_results_in_converged$beta) * 1.1,
             arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    annotate("text", x = -0.0018, y = min(regional_results_in_converged$beta) * 1.15,
             label = "Receivers", size = 5, fontface = "bold") +
    annotate("segment", x = 0.001, xend = 0.003, 
             y = min(regional_results_in_converged$beta) * 1.1, 
             yend = min(regional_results_in_converged$beta) * 1.1,
             arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    annotate("text", x = 0.002, y = min(regional_results_in_converged$beta) * 1.15,
             label = "Senders", size = 5, fontface = "bold") +
    annotate("segment", x = min(regional_results_in_converged$mean_hierarchy) * 1.1, 
             xend = min(regional_results_in_converged$mean_hierarchy) * 1.1,
             y = 0, yend = min(regional_results_in_converged$beta) * 0.8,
             arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    annotate("text", x = min(regional_results_in_converged$mean_hierarchy) * 1.05, 
             y = min(regional_results_in_converged$beta) * 0.4,
             label = "Receives less", size = 5, fontface = "bold", angle = 90) +
    labs(
      x = "Baseline Hierarchy Position (OUT - IN)",
      y = "Psilocin Effect on IN (β)"
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.position = "right"
    )
  
  ggsave(file.path(fig_dir, "hierarchy_drug_effect_IN.png"), 
         p_in_scatter, width = 5, height = 5, dpi = 600)
  
  # Plot 2: OUT beta vs baseline hierarchy
  p_out_scatter <- ggplot(regional_results_out_converged, 
                          aes(x = mean_hierarchy, y = beta)) +
    geom_point(aes(color = Network), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = network_colors, name = "Network") +
    annotate("segment", x = -0.001, xend = -0.003, 
             y = min(regional_results_out_converged$beta) * 1.1, 
             yend = min(regional_results_out_converged$beta) * 1.1,
             arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    annotate("text", x = -0.0018, y = min(regional_results_out_converged$beta) * 1.15,
             label = "Receivers", size = 5, fontface = "bold") +
    annotate("segment", x = 0.001, xend = 0.003, 
             y = min(regional_results_out_converged$beta) * 1.1, 
             yend = min(regional_results_out_converged$beta) * 1.1,
             arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    annotate("text", x = 0.002, y = min(regional_results_out_converged$beta) * 1.15,
             label = "Senders", size = 5, fontface = "bold") +
    annotate("segment", x = max(regional_results_out_converged$mean_hierarchy) * 1.1, 
             xend = max(regional_results_out_converged$mean_hierarchy) * 1.1,
             y = 0, yend = min(regional_results_out_converged$beta) * 0.8,
             arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    annotate("text", x = max(regional_results_out_converged$mean_hierarchy) * 1.05, 
             y = min(regional_results_out_converged$beta) * 0.4,
             label = "Sends less", size = 5, fontface = "bold", angle = 90) +
    labs(
      x = "Baseline Hierarchy Position (OUT - IN)",
      y = "Psilocin Effect on OUT (β)"
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.position = "none"
    )
  
  ggsave(file.path(fig_dir, "hierarchy_drug_effect_OUT.png"), 
         p_out_scatter, width = 5, height = 5, dpi = 600)
  
  # Plot 3: Subject trajectories for extreme regions
  cat("Creating subject trajectory plot...\n")
  
  forward_with_covariates <- hierarchy_results$forward_with_covariates
  reverse_with_covariates <- hierarchy_results$reverse_with_covariates
  
  # Identify top 10% senders and receivers
  n_regions_10pct <- ceiling(232 * 0.10)
  
  top_senders <- hierarchy_positions %>%
    arrange(desc(mean_hierarchy)) %>%
    slice(1:n_regions_10pct) %>%
    pull(Region)
  
  top_receivers <- hierarchy_positions %>%
    arrange(mean_hierarchy) %>%
    slice(1:n_regions_10pct) %>%
    pull(Region)
  
  # Convert to long format
  forward_long <- forward_with_covariates %>%
    select(CIMBI_ID, time_since_admin, time_window, starts_with("Region_")) %>%
    pivot_longer(cols = starts_with("Region_"),
                 names_to = "Region",
                 values_to = "OUT_score") %>%
    mutate(Region = as.numeric(gsub("Region_", "", Region)))
  
  reverse_long <- reverse_with_covariates %>%
    select(CIMBI_ID, time_since_admin, time_window, starts_with("Region_")) %>%
    pivot_longer(cols = starts_with("Region_"),
                 names_to = "Region",
                 values_to = "IN_score") %>%
    mutate(Region = as.numeric(gsub("Region_", "", Region)))
  
  # Calculate hierarchy score
  hierarchy_timeline <- forward_long %>%
    left_join(reverse_long, by = c("CIMBI_ID", "time_since_admin", "time_window", "Region")) %>%
    mutate(hierarchy_score = OUT_score - IN_score)
  
  # Get subject-level averages for top 10% groups
  subject_senders <- hierarchy_timeline %>%
    filter(Region %in% top_senders) %>%
    filter(!is.na(time_window), time_window != "other") %>%
    mutate(time_window = factor(time_window, 
                                levels = c("predrug", "upslope", "peak", "downslope", "return"))) %>%
    group_by(CIMBI_ID, time_window) %>%
    summarise(mean_hierarchy_score = mean(hierarchy_score, na.rm = TRUE), .groups = 'drop') %>%
    mutate(group = "Top 10% Senders")
  
  subject_receivers <- hierarchy_timeline %>%
    filter(Region %in% top_receivers) %>%
    filter(!is.na(time_window), time_window != "other") %>%
    mutate(time_window = factor(time_window, 
                                levels = c("predrug", "upslope", "peak", "downslope", "return"))) %>%
    group_by(CIMBI_ID, time_window) %>%
    summarise(mean_hierarchy_score = mean(hierarchy_score, na.rm = TRUE), .groups = 'drop') %>%
    mutate(group = "Top 10% Receivers")
  
  subject_extreme <- bind_rows(subject_senders, subject_receivers)
  
  # Calculate group means
  summary_stats_subjects <- subject_extreme %>%
    group_by(group, time_window) %>%
    summarise(
      mean_hierarchy = mean(mean_hierarchy_score, na.rm = TRUE),
      se_hierarchy = sd(mean_hierarchy_score, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  p_timeline_subjects <- ggplot() +
    geom_line(data = subject_extreme %>% filter(!CIMBI_ID %in% timeline_plot_exclude), 
              aes(x = time_window, y = mean_hierarchy_score, 
                  group = interaction(CIMBI_ID, group), color = group), 
              alpha = 0.2, linewidth = 0.4) +
    geom_point(data = subject_extreme, 
               aes(x = time_window, y = mean_hierarchy_score, 
                   group = interaction(CIMBI_ID, group), color = group), 
               alpha = 0.5) +
    geom_line(data = summary_stats_subjects, 
              aes(x = time_window, y = mean_hierarchy, color = group, group = group),
              linewidth = 1.5) +
    geom_point(data = summary_stats_subjects, 
               aes(x = time_window, y = mean_hierarchy, color = group),
               size = 5, shape = 21, fill = "white", stroke = 2) +
    scale_color_manual(
      values = c("Top 10% Senders" = "#aa0000", "Top 10% Receivers" = "#005500"),
      name = "Baseline Regional\nHierarchy Group"
    ) +
    labs(
      x = "Time Window",
      y = "Hierarchy Score (OUT - IN)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(fig_dir, "hierarchy_timeline_subjects.png"), 
         p_timeline_subjects, width = 6, height = 6, dpi = 600)
  
  # Export hierarchy data for surface visualization
  hierarchy_surface_export <- hierarchy_positions %>%
    select(Region, mean_OUT, mean_IN, mean_hierarchy) %>%
    mutate(hierarchy_rank = rank(mean_hierarchy)) %>%
    pivot_longer(
      cols = c(mean_OUT, mean_IN, mean_hierarchy, hierarchy_rank),
      names_to = "Measure",
      values_to = "Estimate"
    ) %>%
    mutate(Measure = case_when(
      Measure == "mean_OUT" ~ "Baseline_OUT",
      Measure == "mean_IN" ~ "Baseline_IN",
      Measure == "mean_hierarchy" ~ "Baseline_Hierarchy",
      Measure == "hierarchy_rank" ~ "Hierarchy_Rank"
    )) %>%
    select(Region, Measure, Estimate) %>%
    arrange(Measure, Region)
  
  hierarchy_surface_path <- file.path(out_dir, 'hierarchy_baseline_for_surfaces.csv')
  write_csv(hierarchy_surface_export, hierarchy_surface_path)
  
  cat(sprintf("Hierarchy surface export saved to: %s\n", hierarchy_surface_path))
  cat("Hierarchy plots complete!\n\n")
  
  return(list(
    p_in = p_in_scatter,
    p_out = p_out_scatter,
    p_timeline = p_timeline_subjects
  ))
}

# -----------------------------------------------------------------------------
# REGIONAL INSIDEOUT BOXPLOT BY NETWORK
# -----------------------------------------------------------------------------

plot_regional_insideout_by_network <- function() {
  cat("\n=== PLOTTING REGIONAL INSIDEOUT BY NETWORK ===\n")

  # Load per-region INSIDEOUT (directionality) scores
  dir_scores <- read_csv(file.path(out_dir, "subject_directionality_scores.csv"),
                         show_col_types = FALSE) %>%
    filter(!is.na(psi_conc), !is.nan(psi_conc), psi_conc == 0)

  region_cols_dir <- grep("^Region_", names(dir_scores), value = TRUE)
  mean_insideout  <- colMeans(dir_scores[, region_cols_dir], na.rm = TRUE)

  # Load network labels
  network_mapping <- read_csv(file.path(out_dir, "all_regional_results_for_R.csv"),
                              show_col_types = FALSE) %>%
    select(Region, Network) %>%
    distinct()

  plot_data <- data.frame(
    Region    = seq_along(mean_insideout),
    insideout = as.numeric(mean_insideout)
  ) %>%
    left_join(network_mapping, by = "Region") %>%
    filter(!is.na(Network))

  # Order networks by median INSIDEOUT score
  network_order <- plot_data %>%
    group_by(Network) %>%
    summarise(med = median(insideout, na.rm = TRUE), .groups = "drop") %>%
    arrange(med) %>%
    pull(Network)

  plot_data <- plot_data %>%
    mutate(Network = factor(Network, levels = network_order))

  # Network-level summary for annotation
  network_summary <- plot_data %>%
    group_by(Network) %>%
    summarise(
      n         = n(),
      median_io = median(insideout, na.rm = TRUE),
      mean_io   = mean(insideout,   na.rm = TRUE),
      sd_io     = sd(insideout,     na.rm = TRUE),
      .groups   = "drop"
    )

  p <- ggplot(plot_data, aes(x = Network, y = insideout, fill = Network)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.6,
                 width = 0.6, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.4, color = "grey30") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    labs(
      x     = "Network",
      y     = "Baseline INSIDEOUT score (OUT \u2212 IN)",
      title = "Regional Baseline INSIDEOUT Score by Network"
    ) +
    theme_bw() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 13),
      axis.text.y      = element_text(size = 13),
      axis.title       = element_text(size = 15),
      legend.position  = "none",
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(fig_dir, "regional_insideout_by_network.png"),
         p, width = 10, height = 6, dpi = 300)

  write_csv(network_summary,
            file.path(out_dir, "regional_insideout_network_summary.csv"))

  cat(sprintf("Boxplot saved: regional_insideout_by_network.png\n"))
  cat(sprintf("Network summary saved: regional_insideout_network_summary.csv\n"))
  print(network_summary)

  return(p)
}

# -----------------------------------------------------------------------------
# SCANNER SENSITIVITY: mr001-only GCOR regional analysis
# -----------------------------------------------------------------------------

run_gcor_scanner_sensitivity <- function(combined_data, regional_results_full,
                                         scanner_label = "mr001") {
  cat(sprintf("\n=== GCOR SCANNER SENSITIVITY ANALYSIS (%s only) ===\n", scanner_label))

  # Subjects from this scanner
  mr_subjects <- combined_data %>%
    filter(as.character(scanner) == scanner_label) %>%
    pull(CIMBI_ID) %>%
    as.character() %>%
    unique()
  cat(sprintf("  Subjects on %s: %d\n", scanner_label, length(mr_subjects)))

  # Load regional GCOR scores and filter to scanner subjects
  gcor_scores <- read_csv(file.path(out_dir, "subject_regional_GCOR_scores.csv"),
                          show_col_types = FALSE) %>%
    filter(!is.na(psi_conc), !is.nan(psi_conc)) %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
    filter(CIMBI_ID %in% mr_subjects)

  # Time-window / scan index (same logic as main script)
  gcor_scores <- gcor_scores %>%
    mutate(time_window = case_when(
      time_since_admin < 3                                  ~ "predrug",
      time_since_admin >= 47  & time_since_admin < 95      ~ "upslope",
      time_since_admin >= 95  & time_since_admin < 140     ~ "peak",
      time_since_admin >= 140 & time_since_admin < 250     ~ "downslope",
      time_since_admin >= 282                               ~ "return",
      TRUE ~ "other"
    )) %>%
    filter(time_window != "other") %>%
    group_by(CIMBI_ID, time_window) %>%
    arrange(time_since_admin) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(scan = case_when(
      time_window == "predrug"   ~ 0L,
      time_window == "upslope"   ~ 1L,
      time_window == "peak"      ~ 2L,
      time_window == "downslope" ~ 3L,
      time_window == "return"    ~ 4L
    ))

  # Join covariates (no scanner covariate — single scanner)
  covariates <- combined_data %>%
    mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
    select(CIMBI_ID, time_since_admin, FD, Age, Sex) %>%
    filter(CIMBI_ID %in% mr_subjects)

  gcor_scores <- gcor_scores %>%
    left_join(covariates, by = c("CIMBI_ID", "time_since_admin"))

  region_cols_g <- grep("^Region_", names(gcor_scores), value = TRUE)
  n_reg         <- length(region_cols_g)
  n_subj_mr     <- length(mr_subjects)

  # Network labels from full results
  network_mapping <- regional_results_full %>%
    filter(Measure == "Regional_GCOR") %>%
    select(Region, Network) %>%
    distinct()

  # Fit per-region MMRM (no scanner term)
  estimates_mr <- rep(NA_real_, n_reg)
  cat(sprintf("  Fitting %d regional MMRMs...\n", n_reg))

  for (reg in seq_len(n_reg)) {
    region_data <- gcor_scores %>%
      mutate(
        Score    = .data[[region_cols_g[reg]]],
        CIMBI_ID = as.factor(CIMBI_ID),
        scan     = as.factor(scan),
        Sex      = as.factor(Sex)
      ) %>%
      select(CIMBI_ID, scan, Score, psi_conc, FD, Age, Sex) %>%
      na.omit()

    tryCatch({
      model <- suppressWarnings(mmrm(
        Score ~ psi_conc + FD + Age + Sex + us(scan | CIMBI_ID),
        data = region_data
      ))
      if (model$opt_details$convergence == 0) {
        cf  <- summary(model)$coefficients
        row <- which(rownames(cf) == "psi_conc")
        if (length(row) > 0) estimates_mr[reg] <- cf[row, "Estimate"]
      }
    }, error = function(e) {})

    # toeph fallback
    if (is.na(estimates_mr[reg])) {
      tryCatch({
        model <- suppressWarnings(mmrm(
          Score ~ psi_conc + FD + Age + Sex + toeph(scan | CIMBI_ID),
          data = region_data
        ))
        if (model$opt_details$convergence == 0) {
          cf  <- summary(model)$coefficients
          row <- which(rownames(cf) == "psi_conc")
          if (length(row) > 0) estimates_mr[reg] <- cf[row, "Estimate"]
        }
      }, error = function(e) {})
    }
  }

  n_converged <- sum(!is.na(estimates_mr))
  cat(sprintf("  Converged: %d / %d regions\n", n_converged, n_reg))

  # --- Global GCOR x PPL: MMRM on mean-across-regions score ---
  cat("\n  --- Global GCOR (mean across regions) ---\n")
  global_gcor <- gcor_scores %>%
    mutate(Global_GCOR = rowMeans(across(all_of(region_cols_g)), na.rm = TRUE),
           CIMBI_ID    = as.factor(CIMBI_ID),
           scan        = as.factor(scan),
           Sex         = as.factor(Sex)) %>%
    select(CIMBI_ID, scan, Global_GCOR, psi_conc, FD, Age, Sex) %>%
    na.omit()

  global_summary <- tryCatch({
    gm <- suppressWarnings(mmrm(
      Global_GCOR ~ psi_conc + FD + Age + Sex + us(scan | CIMBI_ID),
      data = global_gcor
    ))
    if (gm$opt_details$convergence != 0) stop("us() did not converge")
    gm
  }, error = function(e) {
    suppressWarnings(mmrm(
      Global_GCOR ~ psi_conc + FD + Age + Sex + toeph(scan | CIMBI_ID),
      data = global_gcor
    ))
  })

  gcf <- summary(global_summary)$coefficients
  grow <- which(rownames(gcf) == "psi_conc")
  global_beta <- gcf[grow, "Estimate"]
  global_se   <- gcf[grow, "Std. Error"]
  global_t    <- gcf[grow, "t value"]
  global_p    <- gcf[grow, "Pr(>|t|)"]
  cat(sprintf("  Global GCOR ~ psi_conc: beta = %.4f, SE = %.4f, t = %.3f, p = %.4f\n",
              global_beta, global_se, global_t, global_p))

  # Compare to full-dataset betas
  full_betas <- regional_results_full %>%
    filter(Measure == "Regional_GCOR") %>%
    arrange(Region) %>%
    pull(Estimate)

  compare_df <- data.frame(
    Region       = seq_len(n_reg),
    beta_full    = full_betas,
    beta_mr001   = estimates_mr
  ) %>%
    left_join(network_mapping, by = "Region") %>%
    filter(!is.na(beta_full), !is.na(beta_mr001))

  r_cor <- cor(compare_df$beta_full, compare_df$beta_mr001)
  cat(sprintf("  Correlation between full and %s betas: r = %.3f\n", scanner_label, r_cor))

  # --- Plot 1: scatter of full vs mr001 betas, coloured by network ---
  p_scatter <- ggplot(compare_df, aes(x = beta_full, y = beta_mr001, color = Network)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "black", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("r = %.3f", r_cor),
             hjust = -0.2, vjust = 1.5, size = 5) +
    scale_color_manual(values = network_colors) +
    labs(x = "Beta — full cohort",
         y = sprintf("Beta — %s only", scanner_label),
         title = sprintf("Regional GCOR betas: full cohort vs %s", scanner_label)) +
    theme_bw() +
    theme(plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
          axis.title  = element_text(size = 14),
          axis.text   = element_text(size = 12),
          legend.text = element_text(size = 10))

  ggsave(file.path(fig_dir, sprintf("gcor_scanner_sensitivity_%s_scatter.png", scanner_label)),
         p_scatter, width = 7, height = 6, dpi = 300)

  # --- Plot 2: network-level mean betas, full vs mr001 ---
  net_summary <- compare_df %>%
    filter(Network %in% c("Default", "Cont", "DorsAttn", "SalVentAttn", "Vis", "SomMot")) %>%
    group_by(Network) %>%
    summarise(mean_full  = mean(beta_full),
              mean_mr001 = mean(beta_mr001),
              .groups = "drop") %>%
    pivot_longer(cols = c(mean_full, mean_mr001),
                 names_to = "cohort", values_to = "mean_beta") %>%
    mutate(
      cohort  = ifelse(cohort == "mean_full", "Full cohort", scanner_label),
      Network = factor(Network, levels = c("Default","Cont","DorsAttn",
                                           "SalVentAttn","Vis","SomMot"))
    )

  p_net <- ggplot(net_summary, aes(x = Network, y = mean_beta,
                                    fill = cohort, group = cohort)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("Full cohort" = "#4477AA", "mr001" = "#EE6677"),
                      name = NULL) +
    labs(x = NULL, y = "Mean beta (psi_conc)",
         title = sprintf("Network-level GCOR betas: full cohort vs %s", scanner_label)) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 35, hjust = 1, size = 12),
          axis.text.y  = element_text(size = 12),
          axis.title   = element_text(size = 14),
          plot.title   = element_text(size = 13, face = "bold", hjust = 0.5),
          legend.text  = element_text(size = 12))

  ggsave(file.path(fig_dir, sprintf("gcor_scanner_sensitivity_%s_networks.png", scanner_label)),
         p_net, width = 8, height = 5, dpi = 300)

  # --- Network-level beta summary table ---
  net_summary_full <- compare_df %>%
    group_by(Network) %>%
    summarise(
      n_regions      = n(),
      mean_beta_full = mean(beta_full,  na.rm = TRUE),
      mean_beta_mr   = mean(beta_mr001, na.rm = TRUE),
      sd_beta_full   = sd(beta_full,    na.rm = TRUE),
      sd_beta_mr     = sd(beta_mr001,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_beta_mr))

  cat("\n  --- Network-level mean betas (full vs scanner-subset) ---\n")
  print(as.data.frame(net_summary_full))

  # Save comparison table and network summary
  write_csv(compare_df,
            file.path(out_dir, sprintf("gcor_scanner_sensitivity_%s.csv", scanner_label)))
  write_csv(net_summary_full,
            file.path(out_dir, sprintf("gcor_scanner_sensitivity_%s_network_summary.csv", scanner_label)))

  cat(sprintf("\n  Plots and CSVs saved for %s sensitivity analysis.\n", scanner_label))
  return(list(
    compare_df        = compare_df,
    network_summary   = net_summary_full,
    r_full_vs_scanner = r_cor,
    global_beta       = global_beta,
    global_se         = global_se,
    global_t          = global_t,
    global_p          = global_p,
    plot_scatter      = p_scatter,
    plot_networks     = p_net
  ))
}

# -----------------------------------------------------------------------------
# TRANSMODAL vs UNIMODAL GCOR BETA PERMUTATION TEST
# -----------------------------------------------------------------------------

test_transmodal_vs_unimodal <- function(regional_results, measure = "Regional_GCOR", n_perm = 10000,
                                        direction = "greater") {
  # direction = "greater": test transmodal > unimodal (GCOR)
  # direction = "less":    test transmodal < unimodal (Directionality — more negative)
  cat(sprintf("\n=== TRANSMODAL vs UNIMODAL PERMUTATION TEST: %s ===\n", measure))

  transmodal_nets <- c("Default", "Cont", "DorsAttn", "SalVentAttn")
  unimodal_nets   <- c("Vis", "SomMot")
  # Limbic and SubC excluded

  measure_betas <- regional_results %>%
    filter(Measure == measure, !is.na(Estimate)) %>%
    filter(Network %in% c(transmodal_nets, unimodal_nets)) %>%
    mutate(group = ifelse(Network %in% transmodal_nets, "transmodal", "unimodal"))

  n_trans    <- sum(measure_betas$group == "transmodal")
  n_uni      <- sum(measure_betas$group == "unimodal")
  mean_trans <- mean(measure_betas$Estimate[measure_betas$group == "transmodal"])
  mean_uni   <- mean(measure_betas$Estimate[measure_betas$group == "unimodal"])
  obs_diff   <- mean_trans - mean_uni

  cat(sprintf("  Transmodal regions: %d  |  mean beta = %.6f\n", n_trans, mean_trans))
  cat(sprintf("  Unimodal regions:   %d  |  mean beta = %.6f\n", n_uni,   mean_uni))
  cat(sprintf("  Observed difference (transmodal - unimodal): %.6f\n", obs_diff))

  # Permutation test: shuffle group labels, recompute mean difference
  set.seed(42)
  all_betas  <- measure_betas$Estimate
  perm_diffs <- replicate(n_perm, {
    perm_labels <- sample(measure_betas$group)
    mean(all_betas[perm_labels == "transmodal"]) -
      mean(all_betas[perm_labels == "unimodal"])
  })

  # One-sided p-value
  if (direction == "greater") {
    p_value  <- mean(perm_diffs >= obs_diff)
    alt_label <- "transmodal > unimodal"
  } else {
    p_value  <- mean(perm_diffs <= obs_diff)
    alt_label <- "transmodal < unimodal"
  }

  cat(sprintf("  Permutation p-value (one-sided, %s): %.4f\n", alt_label, p_value))
  cat(sprintf("  (based on %d permutations)\n", n_perm))

  # Summary by network
  network_summary <- measure_betas %>%
    group_by(Network, group) %>%
    summarise(n = n(), mean_beta = mean(Estimate), sd_beta = sd(Estimate), .groups = "drop") %>%
    arrange(group, desc(mean_beta))
  cat("\nPer-network summary:\n")
  print(network_summary)

  # Plot: histogram of permutation distribution with observed value
  file_slug <- tolower(gsub("[^A-Za-z0-9]", "_", measure))
  p_hist <- ggplot(data.frame(diff = perm_diffs), aes(x = diff)) +
    geom_histogram(bins = 80, fill = "grey70", color = "white") +
    geom_vline(xintercept = obs_diff, color = "firebrick", linewidth = 1.2) +
    annotate("text", x = obs_diff, y = Inf,
             label = sprintf("Observed\np = %.4f", p_value),
             hjust = -0.1, vjust = 1.5, color = "firebrick", size = 5) +
    labs(x = sprintf("Mean(transmodal) \u2212 Mean(unimodal) beta  [%s]", alt_label),
         y = "Count",
         title = sprintf("Permutation Distribution: Transmodal vs Unimodal\n%s", measure)) +
    theme_bw() +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text  = element_text(size = 12))

  ggsave(file.path(fig_dir, sprintf("%s_transmodal_vs_unimodal_permtest.png", file_slug)),
         p_hist, width = 8, height = 5, dpi = 300)

  results <- data.frame(
    measure              = measure,
    mean_beta_transmodal = mean_trans,
    mean_beta_unimodal   = mean_uni,
    observed_difference  = obs_diff,
    n_transmodal         = n_trans,
    n_unimodal           = n_uni,
    n_permutations       = n_perm,
    p_value_one_sided    = p_value
  )
  write_csv(results, file.path(out_dir, sprintf("%s_transmodal_vs_unimodal_permtest.csv", file_slug)))

  return(list(results = results, network_summary = network_summary, plot = p_hist))
}

# -----------------------------------------------------------------------------
# 8b. NETWORK x NETWORK CONNECTIVITY HEATMAP
# -----------------------------------------------------------------------------

plot_network_connectivity_heatmap <- function(
  edge_file = file.path(project_dir, "outputs/2025-10-07/edge_GCOR/results/edge_results_mmrm.csv"),
  fig_dir   = NULL,
  out_dir   = NULL
) {
  cat("\n=== NETWORK x NETWORK CONNECTIVITY HEATMAP ===\n")

  if (!file.exists(edge_file)) stop("Edge results file not found: ", edge_file)

  edge_results <- read_csv(edge_file, show_col_types = FALSE)
  cat(sprintf("  Loaded %d edges\n", nrow(edge_results)))

  # Network display order: unimodal → transmodal, then SubC/Limbic
  net_order <- c("Default", "Cont", "SalVentAttn", "DorsAttn", "SomMot", "Vis", "Limbic", "SubC")
  net_order <- net_order[net_order %in% union(edge_results$from_network, edge_results$to_network)]

  # Mean beta per (from_network, to_network) pair — make symmetric by averaging both directions
  pair_summary <- edge_results %>%
    mutate(
      net_a = pmin(from_network, to_network),
      net_b = pmax(from_network, to_network),
      within = from_network == to_network
    ) %>%
    group_by(net_a, net_b, within) %>%
    summarise(mean_beta = mean(beta, na.rm = TRUE),
              n_edges   = n(),
              .groups   = "drop")

  # Expand to full symmetric matrix (both [A,B] and [B,A])
  mat_long <- bind_rows(
    pair_summary %>% rename(from_network = net_a, to_network = net_b),
    pair_summary %>% filter(net_a != net_b) %>% rename(from_network = net_b, to_network = net_a)
  ) %>%
    mutate(
      from_network = factor(from_network, levels = net_order),
      to_network   = factor(to_network,   levels = rev(net_order))
    )

  # Colour scale limits — symmetric around 0, clipped at 95th pct of abs values
  lim <- quantile(abs(mat_long$mean_beta), 0.95, na.rm = TRUE)

  # Label formatting
  mat_long <- mat_long %>%
    mutate(label = sprintf("%.1f", mean_beta))

  p_heat <- ggplot(mat_long, aes(x = from_network, y = to_network, fill = mean_beta)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 4.5, fontface = "bold", color = "black") +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#d6604d",
      midpoint = 0,
      limits   = c(-lim, lim),
      oob      = scales::squish,
      name     = "Mean \u03B2\n(PPL)"
    ) +
    scale_x_discrete(position = "bottom") +
    labs(
      x     = NULL,
      y     = NULL,
      title = "Network \u00D7 Network Connectivity: Mean \u03B2 (psilocin concentration)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y     = element_text(size = 11),
      legend.title    = element_text(size = 10),
      panel.grid      = element_blank()
    )

  # Save
  if (!is.null(fig_dir)) {
    heat_path <- file.path(fig_dir, "gcor_network_x_network_heatmap.png")
    ggsave(heat_path, p_heat, width = 7, height = 6, dpi = 300)
    cat(sprintf("  Heatmap saved to: %s\n", heat_path))
  }

  # Also save the summary table
  if (!is.null(out_dir)) {
    tbl_path <- file.path(out_dir, "gcor_network_x_network_mean_beta.csv")
    write_csv(pair_summary %>%
                mutate(connection_type = ifelse(within, "Within", "Between")) %>%
                select(net_a, net_b, connection_type, n_edges, mean_beta),
              tbl_path)
    cat(sprintf("  Table saved to: %s\n", tbl_path))
  }

  cat("\n  Network x Network mean beta:\n")
  print(pair_summary %>%
          mutate(connection_type = ifelse(within, "Within", "Between")) %>%
          select(net_a, net_b, connection_type, n_edges, mean_beta) %>%
          arrange(connection_type, desc(mean_beta)))

  invisible(list(plot = p_heat, data = pair_summary))
}

# -----------------------------------------------------------------------------
# 9. MAIN ANALYSIS PIPELINE
# -----------------------------------------------------------------------------

# Load data
cat("Loading and preparing data...\n")
combined_data <- load_and_prepare_data()


# Define variables - only absolute modularity
variables <- c(
  "insideout",
  "modularity_fc_abs",
  "modularity_forward_MI_abs",
  "global_GCOR"
)

var_names <- c(
  "insideout" = "INSIDEOUT",
  "modularity_fc_abs" = "Modularity (FC)",
  "modularity_forward_MI_abs" = "Modularity (Fwd MI)",
  "global_GCOR" = "Global GCOR"
)

# Check which variables exist
available_vars <- variables[variables %in% names(combined_data)]
cat(sprintf("\nFound %d/%d variables\n", length(available_vars), length(variables)))

if(length(available_vars) < length(variables)) {
  missing_vars <- setdiff(variables, available_vars)
  cat("Missing variables:\n")
  print(missing_vars)
}

# -----------------------------------------------------------------------------
# 10. TIMELINE ANALYSIS
# -----------------------------------------------------------------------------

cat("\nRunning timeline analysis...\n")
timeline_results <- list()

for(var in available_vars) {
  cat(sprintf("  Processing %s...\n", var))
  result <- analyze_timeline_variable(combined_data, var, var_names[[var]])
  
  ggsave(file.path(fig_dir, paste0("boxplot_", var, "_time.png")), 
         result$plot, width = 6, height = 6, dpi = 600)
  
  timeline_results[[var]] <- result$results
}

timeline_results_df <- do.call(rbind, timeline_results)
write_xlsx(timeline_results_df, file.path(fig_dir, "timeline_analysis_results.xlsx"))

# -----------------------------------------------------------------------------
# 11. REGRESSION ANALYSIS
# -----------------------------------------------------------------------------

cat("\nRunning regression analysis...\n")
predictors <- c("psi_conc")
regression_results <- list()
regression_models <- list()

for(var in available_vars) {
  for(pred in predictors) {
    cat(sprintf("  %s vs %s...\n", var, pred))
    
    tryCatch({
      result <- run_mixed_model(combined_data, var, pred, var_names)
      
      ggsave(file.path(fig_dir, paste0("regression_", var, "_vs_", pred, ".png")),
             result$plot, width = 4.5, height = 4, dpi = 300)
      
      regression_results[[paste(var, pred, sep = "_")]] <- result$results
      regression_models[[paste(var, pred, sep = "_")]] <- list(
        model = result$model,
        model_data = result$model_data
      )
    }, error = function(e) {
      warning(sprintf("Error in %s vs %s: %s", var, pred, e$message))
    })
  }
}

regression_results_df <- do.call(rbind, regression_results)
write_xlsx(regression_results_df, file.path(fig_dir, "regression_analysis_results.xlsx"))

# -----------------------------------------------------------------------------
# 12. NULL DISTRIBUTIONS VS GCOR
# -----------------------------------------------------------------------------

cat("\nPlotting null distributions vs GCOR...\n")
null_dist_plots <- plot_null_distributions_vs_gcor(combined_data, out_dir, fig_dir)

# -----------------------------------------------------------------------------
# 13. REGIONAL ANALYSIS
# -----------------------------------------------------------------------------

cat("\n=== RUNNING REGIONAL ANALYSIS ===\n")
regional_results <- run_regional_analysis()
create_regional_plots(regional_results, out_dir)
transmodal_test_gcor <- test_transmodal_vs_unimodal(regional_results, measure = "Regional_GCOR")
transmodal_test_dir  <- test_transmodal_vs_unimodal(regional_results, measure = "Directionality", direction = "less")
scanner_sensitivity  <- run_gcor_scanner_sensitivity(combined_data, regional_results, scanner_label = "N")
network_heatmap      <- plot_network_connectivity_heatmap(fig_dir = fig_dir, out_dir = out_dir)

# -----------------------------------------------------------------------------
# 14. HIERARCHY ANALYSIS
# -----------------------------------------------------------------------------

cat("\n=== RUNNING HIERARCHY ANALYSIS ===\n")
hierarchy_data <- calculate_hierarchy_positions()
hierarchy_results <- run_hierarchy_analysis(hierarchy_data)
hierarchy_plots <- create_hierarchy_plots(hierarchy_results, hierarchy_data$hierarchy_positions)
insideout_network_plot <- plot_regional_insideout_by_network()

# -----------------------------------------------------------------------------
# 15. CORRELATION HEATMAP
# -----------------------------------------------------------------------------

cat("\nCreating correlation heatmap...\n")
correlation_results <- create_correlation_heatmap(combined_data, fig_dir, out_dir)

# -----------------------------------------------------------------------------
# 16. CREATE SUPPLEMENTARY TABLES
# -----------------------------------------------------------------------------

# cat("\n=== CREATING SUPPLEMENTARY TABLES ===\n")
# (supplementary-table export not included in this repository) # source(file.path(project_dir, "code/create_supplementary_tables2.R"))

# create_supplementary_tables(
#   regression_results_df = regression_results_df,
#   regional_results = regional_results,
#   hierarchy_cor_summary = read_csv(file.path(out_dir, "hierarchy_correlation_summary.csv")),
#   global_cor_results = correlation_results$results,
#   out_dir = out_dir
# )

# -----------------------------------------------------------------------------
# 17. SUMMARY REPORTS
# -----------------------------------------------------------------------------

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Processed %d variables\n", length(available_vars)))
cat(sprintf("Created %d timeline plots\n", length(timeline_results)))
cat(sprintf("Created %d regression plots\n", length(regression_results)))
cat(sprintf("Created 2 hierarchy scatter plots\n"))
cat(sprintf("Created 1 hierarchy timeline plot\n"))
cat(sprintf("\nAll results saved to: %s\n", fig_dir))

# Final verification
cat("\n=== DATA INTEGRITY CHECK ===\n")
cat(sprintf("Total subjects: %d\n", length(unique(combined_data$CIMBI_ID))))
cat(sprintf("Total scans: %d\n", nrow(combined_data)))

# Export for MATLAB surface mapping
surface_export <- regional_results %>%
  filter(Measure %in% c("Directionality", "Regional_GCOR", "Forward_MI", "Reverse_MI")) %>%
  select(Region, Measure, Estimate, PValue) %>%
  arrange(Measure, Region)

surface_csv_path <- file.path(out_dir, 'R_regional_results_for_surfaces.csv')
write.csv(surface_export, surface_csv_path, row.names = FALSE)
cat(sprintf("\nSurface export saved to: %s\n", surface_csv_path))

cat("\n=== ALL ANALYSES COMPLETE ===\n")
