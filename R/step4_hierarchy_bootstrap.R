# ---- paths: set these environment variables (see README) ----
project_dir <- Sys.getenv("UPF_PROJECT_DIR")  # project root (data/, outputs/)
np2_dir     <- Sys.getenv("NP2_DATA_DIR")     # NeuroPharm2 derived data (restricted access)
roi_dir     <- Sys.getenv("ROI_DIR")          # Schaefer atlas files

# ============================================================================
# BOOTSTRAP HIERARCHY ANALYSIS
#
# For each of 100 bootstrap iterations (subjects sampled with replacement):
#   1. Estimate baseline hierarchy position per region from baseline scans only
#      (mean OUT - mean IN at psi_conc == 0 for bootstrapped subjects)
#   2. Estimate drug effect betas (psi_conc coefficient) per region using all
#      scans via MMRM: Score ~ psi_conc + FD + Age + Sex + scanner + us(scan | CIMBI_ID)
#   3. Fit lm(beta ~ hierarchy_position) across regions; save all LM parameters
#
# Outputs:
#   bootstrap_hierarchy_lm_params.xlsx  - LM parameters per bootstrap (sheets: IN, OUT, Combined)
#   bootstrap_hierarchy_regional_data.csv - all regional betas and hierarchy positions per bootstrap
# ============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(writexl)
require(mmrm)

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------

base_out_dir <- file.path(project_dir, "outputs")

date_folders <- list.dirs(base_out_dir, full.names = FALSE, recursive = FALSE)
date_folders <- date_folders[grepl("^\\d{4}-\\d{2}-\\d{2}$", date_folders)]
if(length(date_folders) == 0) stop("No date-specific output folders found in ", base_out_dir)

most_recent_date <- sort(date_folders, decreasing = TRUE)[2]  ## change this index if needed
out_dir <- file.path(base_out_dir, most_recent_date)

cat(sprintf("Using output folder: %s\n\n", out_dir))

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading forward and reverse MI scores...\n")
forward_scores <- read_csv(file.path(out_dir, "subject_forward_MI_scores.csv"), show_col_types = FALSE)
reverse_scores <- read_csv(file.path(out_dir, "subject_reverse_MI_scores.csv"), show_col_types = FALSE)

# Remove rows with missing psi_conc
forward_scores <- forward_scores %>% filter(!is.na(psi_conc) & !is.nan(psi_conc))
reverse_scores <- reverse_scores %>% filter(!is.na(psi_conc) & !is.nan(psi_conc))

# Load covariates
cat("Loading covariates...\n")
combined_data <- readxl::read_xlsx(file.path(out_dir, "combined_scanner_metrics.xlsx")) %>%
  mutate(CIMBI_ID = as.character(CIMBI_ID))

copbet_path  <- file.path(np2_dir, 'CopBET/')
copbet_files <- list.files(path = copbet_path, pattern = "^CopBETtbl.*\\.csv$", full.names = TRUE)
if(length(copbet_files) == 0) stop("FD data not found!")

fd_data <- read.csv(copbet_files[1]) %>%
  select(CIMBI, time_since_admin, FD) %>%
  rename(CIMBI_ID = CIMBI) %>%
  mutate(CIMBI_ID = as.character(CIMBI_ID))

combined_metadata <- combined_data %>%
  select(CIMBI_ID, time_since_admin, Age, Sex, scanner) %>%
  mutate(CIMBI_ID = as.character(CIMBI_ID))

covariates <- combined_metadata %>%
  left_join(fd_data, by = c("CIMBI_ID", "time_since_admin"))

# -----------------------------------------------------------------------------
# PREPROCESS SCORES (same time-window logic as main script)
# -----------------------------------------------------------------------------

preprocess_scores <- function(scores) {
  scores %>%
    mutate(time_window = case_when(
      time_since_admin < 3 ~ "predrug",
      time_since_admin >= 47  & time_since_admin < 95  ~ "upslope",
      time_since_admin >= 95  & time_since_admin < 140 ~ "peak",
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
      time_window == "predrug"   ~ 0,
      time_window == "upslope"   ~ 1,
      time_window == "peak"      ~ 2,
      time_window == "downslope" ~ 3,
      time_window == "return"    ~ 4
    ))
}

forward_proc <- preprocess_scores(forward_scores) %>%
  mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
  left_join(covariates, by = c("CIMBI_ID", "time_since_admin"))

reverse_proc <- preprocess_scores(reverse_scores) %>%
  mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
  left_join(covariates, by = c("CIMBI_ID", "time_since_admin"))

region_cols <- grep("^Region_", names(forward_proc), value = TRUE)
n_regions   <- 232
all_subjects <- unique(forward_proc$CIMBI_ID)
n_subjects   <- length(all_subjects)

cat(sprintf("Subjects: %d | Regions: %d\n\n", n_subjects, n_regions))

# -----------------------------------------------------------------------------
# HELPER: fit MMRM for all regions, return psi_conc betas
# -----------------------------------------------------------------------------

fit_all_regions <- function(boot_data) {
  boot_data <- boot_data %>%
    mutate(
      CIMBI_ID = as.factor(CIMBI_ID),
      scan     = as.factor(scan),
      Sex      = as.factor(Sex),
      scanner  = as.factor(scanner)
    )

  betas      <- rep(NA_real_, n_regions)
  partial_rs <- rep(NA_real_, n_regions)
  converged  <- rep(FALSE,    n_regions)
  cov_used   <- rep(NA_character_, n_regions)

  extract_psi_stats <- function(model, n_subj) {
    coef_table <- summary(model)$coefficients
    psi_row    <- which(rownames(coef_table) == "psi_conc")
    if(length(psi_row) == 0) return(list(beta = NA_real_, partial_r = NA_real_))
    beta      <- coef_table[psi_row, "Estimate"]
    t_val     <- coef_table[psi_row, "t value"]
    # Lipsitz et al. (2001) eq. 8: rho = sign(beta) * t / sqrt(N + t^2)
    partial_r <- sign(beta) * t_val / sqrt(n_subj + t_val^2)
    list(beta = beta, partial_r = partial_r)
  }

  for(reg in 1:n_regions) {
    col_name    <- sprintf("Region_%d", reg)
    region_data <- boot_data %>%
      mutate(Score = .data[[col_name]]) %>%
      select(CIMBI_ID, scan, Score, psi_conc, FD, Age, Sex, scanner) %>%
      na.omit()

    # 1st attempt: us() — same as original analysis
    tryCatch({
      model <- suppressWarnings(mmrm(
        Score ~ psi_conc + FD + Age + Sex + scanner + us(scan | CIMBI_ID),
        data = region_data
      ))
      if(model$opt_details$convergence == 0) {
        stats <- extract_psi_stats(model, length(unique(region_data$CIMBI_ID)))
        if(!is.na(stats$beta)) {
          betas[reg]      <- stats$beta
          partial_rs[reg] <- stats$partial_r
          converged[reg]  <- TRUE
          cov_used[reg]   <- "us"
        }
      }
    }, error = function(e) {})

    # 2nd attempt: toeph() fallback — bootstrapped duplicates can make us() ill-conditioned
    if(!converged[reg]) {
      tryCatch({
        model <- suppressWarnings(mmrm(
          Score ~ psi_conc + FD + Age + Sex + scanner + toeph(scan | CIMBI_ID),
          data = region_data
        ))
        if(model$opt_details$convergence == 0) {
          stats <- extract_psi_stats(model, length(unique(region_data$CIMBI_ID)))
          if(!is.na(stats$beta)) {
            betas[reg]      <- stats$beta
            partial_rs[reg] <- stats$partial_r
            converged[reg]  <- TRUE
            cov_used[reg]   <- "toeph"
          }
        }
      }, error = function(e) {})
    }

    # 3rd attempt: cs() fallback
    if(!converged[reg]) {
      tryCatch({
        model <- suppressWarnings(mmrm(
          Score ~ psi_conc + FD + Age + Sex + scanner + cs(scan | CIMBI_ID),
          data = region_data
        ))
        if(model$opt_details$convergence == 0) {
          stats <- extract_psi_stats(model, length(unique(region_data$CIMBI_ID)))
          if(!is.na(stats$beta)) {
            betas[reg]      <- stats$beta
            partial_rs[reg] <- stats$partial_r
            converged[reg]  <- TRUE
            cov_used[reg]   <- "cs"
          }
        }
      }, error = function(e) {})
    }
  }

  list(betas = betas, partial_rs = partial_rs, converged = converged, cov_used = cov_used)
}

# -----------------------------------------------------------------------------
# BOOTSTRAP
# -----------------------------------------------------------------------------

n_boot <- 3000
all_lm_params    <- vector("list", n_boot)
all_regional_data <- vector("list", n_boot)

set.seed(42)

for(b in 1:n_boot) {
  cat(sprintf("\n--- Bootstrap %d / %d ---\n", b, n_boot))

  # Sample subjects with replacement; give each copy a unique ID for MMRM
  boot_subjects <- sample(all_subjects, n_subjects, replace = TRUE)

  make_boot_data <- function(proc_data) {
    lapply(seq_along(boot_subjects), function(i) {
      proc_data %>%
        filter(CIMBI_ID == boot_subjects[i]) %>%
        mutate(CIMBI_ID = paste0(boot_subjects[i], "_", i))
    }) %>% bind_rows()
  }

  boot_forward <- make_boot_data(forward_proc)
  boot_reverse <- make_boot_data(reverse_proc)

  # --- Hierarchy position: baseline scans only ---
  baseline_fwd  <- boot_forward %>% filter(psi_conc == 0)
  baseline_rev  <- boot_reverse %>% filter(psi_conc == 0)
  boot_mean_OUT <- colMeans(baseline_fwd[, region_cols], na.rm = TRUE)
  boot_mean_IN  <- colMeans(baseline_rev[, region_cols], na.rm = TRUE)
  boot_hierarchy <- as.numeric(boot_mean_OUT - boot_mean_IN)

  # --- Betas: all scans via MMRM ---
  cat("  Fitting IN (Reverse MI) models...\n")
  res_in  <- fit_all_regions(boot_reverse)
  cat("  Fitting OUT (Forward MI) models...\n")
  res_out <- fit_all_regions(boot_forward)

  # --- Linear model: beta ~ hierarchy_position ---
  extract_lm_params <- function(betas, converged, measure) {
    valid <- converged & !is.na(betas) & !is.na(boot_hierarchy)
    fit   <- lm(betas[valid] ~ boot_hierarchy[valid])
    s     <- summary(fit)
    cf    <- s$coefficients
    slope <- cf[2, "Estimate"]
    rho   <- sign(slope) * sqrt(s$r.squared)
    data.frame(
      bootstrap     = b,
      measure       = measure,
      intercept     = cf[1, "Estimate"],
      slope         = slope,
      slope_SE      = cf[2, "Std. Error"],
      slope_t       = cf[2, "t value"],
      slope_p       = cf[2, "Pr(>|t|)"],
      r_squared     = s$r.squared,
      adj_r_squared = s$adj.r.squared,
      rho           = rho,
      F_statistic   = s$fstatistic[1],
      F_p_value     = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE),
      n_regions     = sum(valid)
    )
  }

  params_in  <- extract_lm_params(res_in$betas,  res_in$converged,  "IN")
  params_out <- extract_lm_params(res_out$betas, res_out$converged, "OUT")

  all_lm_params[[b]] <- bind_rows(params_in, params_out)

  all_regional_data[[b]] <- data.frame(
    bootstrap          = b,
    Region             = 1:n_regions,
    hierarchy_position = boot_hierarchy,
    beta_IN            = res_in$betas,
    partial_r_IN       = res_in$partial_rs,
    converged_IN       = res_in$converged,
    cov_used_IN        = res_in$cov_used,
    beta_OUT           = res_out$betas,
    partial_r_OUT      = res_out$partial_rs,
    converged_OUT      = res_out$converged,
    cov_used_OUT       = res_out$cov_used
  )

  n_toeph_in  <- sum(res_in$cov_used  == "toeph", na.rm = TRUE)
  n_toeph_out <- sum(res_out$cov_used == "toeph", na.rm = TRUE)
  n_cs_in     <- sum(res_in$cov_used  == "cs",    na.rm = TRUE)
  n_cs_out    <- sum(res_out$cov_used == "cs",    na.rm = TRUE)
  cat(sprintf("  IN  converged: %d (toeph: %d, cs: %d) | slope_IN  = %.4f\n",
              sum(res_in$converged),  n_toeph_in,  n_cs_in,  params_in$slope))
  cat(sprintf("  OUT converged: %d (toeph: %d, cs: %d) | slope_OUT = %.4f\n",
              sum(res_out$converged), n_toeph_out, n_cs_out, params_out$slope))
}

# -----------------------------------------------------------------------------
# SAVE OUTPUTS
# -----------------------------------------------------------------------------

cat("\n=== SAVING OUTPUTS ===\n")

lm_params_all <- bind_rows(all_lm_params)
regional_all  <- bind_rows(all_regional_data)

lm_params_IN  <- lm_params_all %>% filter(measure == "IN")
lm_params_OUT <- lm_params_all %>% filter(measure == "OUT")

# -----------------------------------------------------------------------------
# FISHER'S Z SUMMARY: atanh(rho) across bootstraps -> CI and p-value
# -----------------------------------------------------------------------------

fisher_summary <- function(rhos, measure_label) {
  z      <- atanh(rhos)
  mean_z <- mean(z, na.rm = TRUE)
  sd_z   <- sd(z,   na.rm = TRUE)
  data.frame(
    measure       = measure_label,
    n_bootstraps  = sum(!is.na(rhos)),
    mean_rho      = tanh(mean_z),
    CI_lower      = tanh(mean_z - 1.96 * sd_z),
    CI_upper      = tanh(mean_z + 1.96 * sd_z),
    p_value       = 2 * (1 - pnorm(abs(mean_z / sd_z))),
    mean_z        = mean_z,
    sd_z          = sd_z
  )
}

fisher_IN  <- fisher_summary(lm_params_IN$rho,  "IN")
fisher_OUT <- fisher_summary(lm_params_OUT$rho, "OUT")
fisher_all <- bind_rows(fisher_IN, fisher_OUT)

cat("\n=== FISHER'S Z SUMMARY ===\n")
print(fisher_all)

xlsx_path <- file.path(out_dir, "bootstrap_hierarchy_lm_params.xlsx")
write_xlsx(
  list(
    IN           = lm_params_IN,
    OUT          = lm_params_OUT,
    Combined     = lm_params_all,
    Fisher_summary = fisher_all
  ),
  xlsx_path
)

regional_path <- file.path(out_dir, "bootstrap_hierarchy_regional_data.csv")
write_csv(regional_all, regional_path)

cat(sprintf("LM parameters saved to:  %s\n", xlsx_path))
cat(sprintf("Regional data saved to:  %s\n", regional_path))
cat("\n=== BOOTSTRAP HIERARCHY ANALYSIS COMPLETE ===\n")
