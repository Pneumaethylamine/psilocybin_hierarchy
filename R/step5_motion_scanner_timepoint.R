# ---- paths: set these environment variables (see README) ----
project_dir <- Sys.getenv("UPF_PROJECT_DIR")  # project root (data/, outputs/)
np2_dir     <- Sys.getenv("NP2_DATA_DIR")     # NeuroPharm2 derived data (restricted access)
roi_dir     <- Sys.getenv("ROI_DIR")          # Schaefer atlas files

# ============================================================================
# PNAS revision analyses (responses to Advanced Science reviews)
#   B1  head motion (SPM realign FD, as used in the original models) by time
#       window, and its relation to psilocin / global metrics
#   B2  PPL x scanner interaction in the global MMRM models
#   B3  TOST equivalence test, predrug vs return
#   B4  number of pairs in predrug-vs-timepoint paired t-tests
# Data loading replicates step3_psilocin_mmrm_models.R
# (data folder outputs/2026-01-14, FD from first CopBET csv). Prints aggregates only.
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(readr); library(mmrm); library(writexl)
})

data_dir <- file.path(project_dir, "outputs/2026-01-14")
out_dir  <- file.path(file.path(project_dir, "outputs"),
                      format(Sys.Date()), "revision_analyses")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- load data
d <- read_xlsx(file.path(data_dir, "combined_scanner_metrics.xlsx")) %>%
  select(-any_of(c("modularity_fc", "modularity_forward_MI",
                   "modularity_fc_abs", "modularity_forward_MI_abs")))
mod_fc  <- read_csv(file.path(data_dir, "modularity_summary_fc.csv"),  show_col_types = FALSE) %>%
  mutate(CIMBI_ID = as.character(CIMBI_ID))
mod_fmi <- read_csv(file.path(data_dir, "modularity_summary_fmi.csv"), show_col_types = FALSE) %>%
  mutate(CIMBI_ID = as.character(CIMBI_ID))
d <- d %>% mutate(CIMBI_ID = as.character(CIMBI_ID)) %>%
  left_join(mod_fc %>% select(CIMBI_ID, time_since_admin, modularity_fc_abs = modularity_abs),
            by = c("CIMBI_ID", "time_since_admin")) %>%
  left_join(mod_fmi %>% select(CIMBI_ID, time_since_admin, modularity_forward_MI_abs = modularity_abs),
            by = c("CIMBI_ID", "time_since_admin"))

d <- d %>% mutate(time_window = case_when(
    time_since_admin < 3 ~ "predrug",
    time_since_admin >= 47 & time_since_admin < 95 ~ "upslope",
    time_since_admin >= 95 & time_since_admin < 140 ~ "peak",
    time_since_admin >= 140 & time_since_admin < 250 ~ "downslope",
    time_since_admin >= 282 ~ "return",
    TRUE ~ "other")) %>%
  filter(time_window != "other") %>%
  group_by(CIMBI_ID, time_window) %>% arrange(time_since_admin) %>% slice(1) %>% ungroup() %>%
  mutate(time_window = factor(time_window, levels = c("predrug","upslope","peak","downslope","return")),
         scan = factor(as.integer(time_window) - 1),
         Sex = factor(Sex), scanner = factor(scanner))

copbet <- sort(list.files(file.path(np2_dir, "CopBET/"), pattern = "^CopBETtbl.*\\.csv$",
                          full.names = TRUE))[1]
fd <- read.csv(copbet) %>% select(CIMBI, time_since_admin, FD) %>%
  rename(CIMBI_ID = CIMBI) %>% mutate(CIMBI_ID = as.character(CIMBI_ID))
d <- d %>% left_join(fd, by = c("CIMBI_ID", "time_since_admin")) %>%
  mutate(CIMBI_ID = factor(CIMBI_ID))

cat(sprintf("FD source: %s\nScans: %d, subjects: %d, scans with FD: %d\n",
            basename(copbet), nrow(d), n_distinct(d$CIMBI_ID), sum(!is.na(d$FD))))
cat("Scanner levels (scans):\n"); print(table(d$scanner))

vars <- c(insideout = "INSIDEOUT", global_GCOR = "Global GCOR",
          modularity_fc_abs = "Modularity (FC)", modularity_forward_MI_abs = "Modularity (Fwd MI)")

fit <- function(f, data) {
  m <- mmrm(as.formula(f), data = data)
  if (m$opt_details$convergence != 0) warning("non-convergence: ", f)
  m
}
coef_row <- function(m, term) {
  cf <- summary(m)$coefficients
  data.frame(term = term, estimate = cf[term, "Estimate"], se = cf[term, "Std. Error"],
             df = cf[term, "df"], t = cf[term, "t value"], p = cf[term, "Pr(>|t|)"])
}

# ------------------------------------------------------------------- B4 + B3
cat("\n================ B4: pairs per predrug-vs-timepoint t-test ================\n")
b4 <- list(); b3 <- list()
for (v in names(vars)) {
  pre <- d %>% filter(time_window == "predrug") %>% select(CIMBI_ID, pre = all_of(v))
  for (tw in c("upslope", "peak", "downslope", "return")) {
    td <- d %>% filter(time_window == tw) %>% select(CIMBI_ID, post = all_of(v)) %>%
      inner_join(pre, by = "CIMBI_ID") %>% filter(!is.na(pre), !is.na(post))
    tt <- t.test(td$post, td$pre, paired = TRUE)
    b4[[paste(v, tw)]] <- data.frame(metric = vars[[v]], comparison = paste("predrug vs", tw),
                                     n_pairs = nrow(td), t = unname(tt$statistic),
                                     df = unname(tt$parameter), p = tt$p.value)
  }
  # ---- B3: TOST predrug vs return (paired), two bound definitions
  ret <- d %>% filter(time_window == "return") %>% select(CIMBI_ID, post = all_of(v)) %>%
    inner_join(pre, by = "CIMBI_ID") %>% filter(!is.na(pre), !is.na(post))
  pk <- d %>% filter(time_window == "peak") %>% select(CIMBI_ID, post = all_of(v)) %>%
    inner_join(pre, by = "CIMBI_ID") %>% filter(!is.na(pre), !is.na(post))
  diff <- ret$post - ret$pre; n <- length(diff); se <- sd(diff) / sqrt(n)
  peak_effect <- abs(mean(pk$post - pk$pre))
  bounds <- c(dz_0.5 = 0.5 * sd(diff), half_peak_effect = 0.5 * peak_effect)
  for (b in names(bounds)) {
    delta <- bounds[[b]]
    p_lower <- pt((mean(diff) + delta) / se, df = n - 1, lower.tail = FALSE)
    p_upper <- pt((mean(diff) - delta) / se, df = n - 1, lower.tail = TRUE)
    ci90 <- mean(diff) + c(-1, 1) * qt(0.95, n - 1) * se
    b3[[paste(v, b)]] <- data.frame(metric = vars[[v]], bound_type = b, n_pairs = n,
      mean_diff = mean(diff), ci90_low = ci90[1], ci90_high = ci90[2],
      bound = delta, peak_vs_predrug_mean_diff = peak_effect, dz = mean(diff) / sd(diff),
      p_TOST = max(p_lower, p_upper), equivalent = max(p_lower, p_upper) < 0.05)
  }
}
b4 <- bind_rows(b4); rownames(b4) <- NULL
print(b4 %>% mutate(across(c(t, p), ~ signif(.x, 3))))
cat("\n================ B3: TOST predrug vs return (alpha = .05) ================\n")
b3 <- bind_rows(b3); rownames(b3) <- NULL
print(b3 %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ------------------------------------------------------------------------ B2
cat("\n================ B2: PPL x scanner interaction ================\n")
b2 <- list()
for (v in names(vars)) {
  md <- d %>% select(all_of(c(v, "psi_conc", "FD", "Age", "Sex", "scanner", "CIMBI_ID", "scan"))) %>% na.omit()
  m0 <- fit(paste(v, "~ psi_conc + FD + Age + Sex + scanner + us(scan | CIMBI_ID)"), md)
  m1 <- fit(paste(v, "~ psi_conc * scanner + FD + Age + Sex + us(scan | CIMBI_ID)"), md)
  cf <- summary(m1)$coefficients; int <- grep("^psi_conc:scanner", rownames(cf), value = TRUE)
  V <- vcov(m1); L <- setNames(rep(0, nrow(cf)), rownames(cf)); L[c("psi_conc", int)] <- 1
  slope2 <- sum(L * coef(m1)); se2 <- sqrt(as.numeric(t(L) %*% V %*% L))
  m2 <- fit(paste(v, "~ psi_conc * scanner + FD + Age + Sex + us(scan | CIMBI_ID)"),
            md %>% mutate(scanner = relevel(scanner, ref = "P")))
  p_other <- coef_row(m2, "psi_conc")$p
  b2[[v]] <- data.frame(metric = vars[[v]], n_scans = nrow(md),
    pooled_slope = coef(m0)[["psi_conc"]], pooled_p = coef_row(m0, "psi_conc")$p,
    slope_ref_scanner = cf["psi_conc", "Estimate"], p_ref = cf["psi_conc", "Pr(>|t|)"],
    slope_other_scanner = slope2, se_other = se2, p_other = p_other,
    interaction = cf[int, "Estimate"], interaction_se = cf[int, "Std. Error"],
    interaction_p = cf[int, "Pr(>|t|)"], interaction_term = int)
}
b2 <- bind_rows(b2); rownames(b2) <- NULL
cat(sprintf("Reference scanner = %s (N = MR2, n=13; P = MR1, n=15)\n", levels(d$scanner)[1]))
print(b2 %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ------------------------------------------------------------------------ B1
cat("\n================ B1: head motion (mean FD, SPM realign) ================\n")
fd_tw <- d %>% group_by(time_window) %>%
  summarise(n = sum(!is.na(FD)), mean_FD = mean(FD, na.rm = TRUE), sd_FD = sd(FD, na.rm = TRUE),
            median_FD = median(FD, na.rm = TRUE), max_FD = max(FD, na.rm = TRUE), .groups = "drop")
fd_tw_sc <- d %>% group_by(scanner, time_window) %>%
  summarise(n = sum(!is.na(FD)), mean_FD = mean(FD, na.rm = TRUE), sd_FD = sd(FD, na.rm = TRUE), .groups = "drop")
print(fd_tw %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(fd_tw_sc %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

md <- d %>% select(FD, psi_conc, Age, Sex, scanner, CIMBI_ID, scan) %>% na.omit()
m_fd <- fit("FD ~ psi_conc + Age + Sex + scanner + us(scan | CIMBI_ID)", md)
m_fd_tw <- fit("FD ~ scan + Age + Sex + scanner + us(scan | CIMBI_ID)", md)
fd_ppl <- coef_row(m_fd, "psi_conc")
cat("\nFD ~ PPL (MMRM, adjusted Age/Sex/scanner):\n"); print(signif(fd_ppl[-1], 3))
cat("\nFD ~ time window contrasts vs predrug:\n")
fd_tw_coef <- bind_rows(lapply(paste0("scan", 1:4), function(x) coef_row(m_fd_tw, x)))
print(fd_tw_coef %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

b1 <- list()
for (v in names(vars)) {
  md <- d %>% select(all_of(c(v, "psi_conc", "FD", "Age", "Sex", "scanner", "CIMBI_ID", "scan"))) %>% na.omit()
  m_full <- fit(paste(v, "~ psi_conc + FD + Age + Sex + scanner + us(scan | CIMBI_ID)"), md)
  m_nofd <- fit(paste(v, "~ psi_conc + Age + Sex + scanner + us(scan | CIMBI_ID)"), md)
  m_fdonly <- fit(paste(v, "~ FD + Age + Sex + scanner + us(scan | CIMBI_ID)"), md)
  b1[[v]] <- data.frame(metric = vars[[v]],
    PPL_slope_withFD = coef(m_full)[["psi_conc"]], PPL_p_withFD = coef_row(m_full, "psi_conc")$p,
    PPL_slope_noFD = coef(m_nofd)[["psi_conc"]], PPL_p_noFD = coef_row(m_nofd, "psi_conc")$p,
    FD_coef_adjPPL = coef(m_full)[["FD"]], FD_p_adjPPL = coef_row(m_full, "FD")$p,
    FD_coef_noPPL = coef(m_fdonly)[["FD"]], FD_p_noPPL = coef_row(m_fdonly, "FD")$p)
}
b1 <- bind_rows(b1); rownames(b1) <- NULL
cat("\nGlobal metrics: PPL effect with/without FD; FD effect with/without PPL:\n")
print(b1 %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

write_xlsx(list(B1_FD_by_window = fd_tw, B1_FD_by_scanner_window = fd_tw_sc,
                B1_FD_vs_PPL = fd_ppl, B1_FD_window_contrasts = fd_tw_coef,
                B1_metric_FD_sensitivity = b1, B2_PPL_x_scanner = b2,
                B3_TOST_predrug_return = b3, B4_ttest_pairs = b4),
           file.path(out_dir, "revision_analyses_results.xlsx"))
cat("\nSaved:", file.path(out_dir, "revision_analyses_results.xlsx"), "\n")
