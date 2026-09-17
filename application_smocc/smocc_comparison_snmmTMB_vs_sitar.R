# ------------------------------------------------------------------------------
# smocc_comparison_snmmTMB_vs_sitar.R
#
# Produces Table 1 and Figure 5 of Section 4
#
# Run after both fits:
#   Rscript smocc_model_snmmTMB.R
#   Rscript smocc_model_sitar.R
#   Rscript smocc_comparison_snmmTMB_vs_sitar.R
#
# Reads  results/smocc_snmmTMB_results.RDS
#        results/smocc_sitar_results.RDS
# Writes Figure5.pdf
#        results/smocc_comparison_summary.txt
# ------------------------------------------------------------------------------

library(TMB)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 0. Configuration
# ==============================================================================

# figure sizes
FIG_WIDTH_IN  <- 6.5
FIG_HEIGHT_IN <- 7.6
PT_AXIS_TITLE <- 11
PT_AXIS_TEXT  <- 10
PT_STRIP      <- 10
PT_LEGEND     <- 10
PT_TAG        <- 12

## Four subjects shown in the figure
main_subjects <- c(1, 2, 3, 4)

col_simultaneous <- "grey65"
fill_simultaneous <- "grey75"
alpha_simultaneous <- 0.22
lty_tmb <- "solid"
lty_sit <- "22"

## Method colours
col_tmb <- "#0072B2"
col_sit <- "#D55E00"
col_values <- c(snmmTMB = col_tmb, sitar = col_sit)
lw_main <- 0.8

fig_theme <- theme_bw(base_size = PT_AXIS_TEXT) +
  theme(
    axis.title      = element_text(size = PT_AXIS_TITLE),
    axis.text       = element_text(size = PT_AXIS_TEXT, colour = "black"),
    strip.text      = element_text(size = PT_STRIP, colour = "black"),
    strip.background = element_rect(fill = "grey92", colour = "grey40"),
    legend.text     = element_text(size = PT_LEGEND),
    legend.title    = element_blank(),
    legend.key.width = unit(1.6, "lines"),
    panel.grid.minor = element_blank(),
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.tag      = element_text(size = PT_TAG, face = "bold")
  )

# ==============================================================================
# 1. Load
# ==============================================================================

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

src_dir  <- "../src"
dll_name <- "snmmTMB_likelihood_smocc_application"

f_tmb <- file.path(results_dir, "smocc_snmmTMB_results.RDS")
f_sit <- file.path(results_dir, "smocc_sitar_results.RDS")

if (!file.exists(f_tmb)) {
  stop("missing ", f_tmb,
       " - run smocc_model_snmmTMB.R first")
}
if (!file.exists(f_sit)) stop("missing ", f_sit, " - run smocc_model_sitar.R first")

tmb_all <- readRDS(f_tmb)
sit <- readRDS(f_sit)

tmb <- tmb_all$snmmTMB

stopifnot(identical(tmb_all$subject_ids, sit$subject_ids))
if (!isTRUE(all.equal(tmb_all$x_grid, sit$x_grid))) {
  stop("the two fits used different reporting grids; set `ngrid` equal in both scripts")
}

subject_ids <- tmb_all$subject_ids

method_levels <- c("snmmTMB", "sitar")
lty_values <- c(snmmTMB = lty_tmb, sitar = lty_sit)
stopifnot(all(main_subjects %in% tmb_all$subject_ids))

# ==============================================================================
# 2. Curve data
# ==============================================================================

pop <- bind_rows(
  tmb$curves$population %>% mutate(method = "snmmTMB"),
  sit$curves$population %>% mutate(method = "sitar")
) %>%
  mutate(method = factor(method, levels = method_levels))

subj <- bind_rows(
  tmb$curves$subjects %>% mutate(method = "snmmTMB"),
  sit$curves$subjects %>% mutate(method = "sitar")
) %>%
  mutate(method = factor(method, levels = method_levels))

# ==============================================================================
# 3. Figure 5, panel (a)
# ==============================================================================

p_pop <- ggplot(pop, aes(x = x)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, group = method, fill = method),
              alpha = alpha_simultaneous, show.legend = FALSE) +
  geom_line(aes(y = lwrS, linetype = method, colour = method),
            linewidth = 0.35, show.legend = FALSE) +
  geom_line(aes(y = uprS, linetype = method, colour = method),
            linewidth = 0.35, show.legend = FALSE) +
  geom_line(aes(y = y, linetype = method, colour = method), linewidth = lw_main) +
  scale_linetype_manual(values = lty_values, name = NULL) +
  scale_colour_manual(values = col_values, name = NULL) +
  scale_fill_manual(values = col_values, guide = "none") +
  labs(x = "Age (weeks)", y = "Height (cm)") +
  fig_theme +
  theme(legend.position = "top")

# ==============================================================================
# 4. Figure 5, panel (b)
# ==============================================================================

pop_tmb <- tmb$curves$population %>% select(x, y_tmb = y, lwrS_tmb = lwrS, uprS_tmb = uprS)
pop_sit <- sit$curves$population %>% select(x, y_sit = y)

diff_df <- inner_join(pop_tmb, pop_sit, by = "x") %>%
  mutate(diff = y_tmb - y_sit,
         halfwidth = (uprS_tmb - lwrS_tmb) / 2)
stopifnot(nrow(diff_df) == nrow(pop_tmb))

max_abs_diff <- max(abs(diff_df$diff))
mean_halfwidth <- mean(diff_df$halfwidth)
frac_within_band <- mean(abs(diff_df$diff) <= diff_df$halfwidth)

p_diff <- ggplot(diff_df, aes(x = x)) +
  geom_ribbon(aes(ymin = -halfwidth, ymax = halfwidth),
              fill = fill_simultaneous, alpha = alpha_simultaneous) +
  geom_line(aes(y = -halfwidth), colour = col_simultaneous, linewidth = 0.35) +
  geom_line(aes(y = halfwidth), colour = col_simultaneous, linewidth = 0.35) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.3) +
  geom_line(aes(y = diff), colour = "black", linewidth = 0.8) +
  labs(x = "Age (weeks)", y = "Difference (cm)") +
  fig_theme

# ==============================================================================
# 5. Figure 5, panel (c)
# ==============================================================================

raw <- tmb_all$bootstrap_inputs$smocc_200 %>%
  mutate(subject = as.integer(as.character(id))) %>%
  filter(subject %in% main_subjects) %>%
  select(subject, age, hgt)

subj_fig <- subj %>%
  filter(subject %in% main_subjects) %>%
  mutate(subject = factor(subject, levels = sort(unique(subject))))

raw <- raw %>% mutate(subject = factor(subject, levels = levels(subj_fig$subject)))

subj_labeller <- as_labeller(function(x) paste("Subject", x))

p_subj <- ggplot(subj_fig, aes(x = x)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, group = method, fill = method),
              alpha = alpha_simultaneous, show.legend = FALSE) +
  geom_line(aes(y = lwrS, linetype = method, colour = method),
            linewidth = 0.3, show.legend = FALSE) +
  geom_line(aes(y = uprS, linetype = method, colour = method),
            linewidth = 0.3, show.legend = FALSE) +
  geom_line(aes(y = y, linetype = method, colour = method), linewidth = lw_main) +
  geom_point(data = raw, aes(x = age, y = hgt), shape = 21, colour = "black",
             fill = "white", size = 2.2, stroke = 0.5) +
  facet_wrap(~ subject, ncol = 2, scales = "fixed", labeller = subj_labeller) +
  scale_linetype_manual(values = lty_values, name = NULL) +
  scale_colour_manual(values = col_values, name = NULL) +
  scale_fill_manual(values = col_values, guide = "none") +
  labs(x = "Age (weeks)", y = "Height (cm)") +
  fig_theme

# ==============================================================================
# 6. Figure 5
# ==============================================================================

figure5 <- (p_pop / p_diff / p_subj) +
  plot_layout(heights = c(1.0, 0.6, 1.5), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top")

ggsave("Figure5.pdf", figure5, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN)

cat(sprintf("max |snmmTMB - sitar| = %.4f cm, mean snmmTMB simultaneous half-width = %.4f cm, within-band fraction = %.3f
",
            max_abs_diff, mean_halfwidth, frac_within_band))

# ==============================================================================
# 7. Effective degrees of freedom of the snmmTMB spline
# ==============================================================================

bi <- tmb_all$bootstrap_inputs
dll_path <- file.path(src_dir, dll_name)
if (!file.exists(dynlib(dll_path))) TMB::compile(paste0(dll_path, ".cpp"))
dyn.load(dynlib(dll_path))

obj_edf <- MakeADFun(data = bi$Data, parameters = bi$parList,
                     random = c("b_intercept", "b_shift", "vpos"),
                     DLL = dll_name, silent = TRUE)
obj_edf$par <- bi$opt_par
invisible(obj_edf$fn(obj_edf$par))

lambda_hat <- exp(as.numeric(bi$opt_par["log_lambda"]))
random_names <- names(obj_edf$env$last.par.best)[obj_edf$env$random]
vpos_pos <- which(random_names == "vpos")
stopifnot(length(vpos_pos) == length(bi$dpos))

A_diag <- diag(solve(as.matrix(obj_edf$env$spHess(random = TRUE))))
edf_tmb <- length(vpos_pos) - sum(A_diag[vpos_pos] * lambda_hat * bi$dpos) + ncol(bi$U0)

message(sprintf(
  "snmmTMB: lambda = %.4f, EDF = %.3f (of %d basis functions: %d penalized + %d null-space)",
  lambda_hat, edf_tmb, bi$K, length(bi$dpos), ncol(bi$U0)))
message(sprintf("sitar: df = %d, chosen by BIC over %s, then held fixed",
                sit$best_df, paste(range(sit$bic_table$spline_df), collapse = "-")))

# ==============================================================================
# 8. Parameter comparison
# ==============================================================================


tmb_term <- function(tab, term) {
  r <- tab[tab$term == term, ]
  if (nrow(r) != 1) return(c(NA_real_, NA_real_, NA_real_))
  c(r$estimate_natural, r$ci_lower_natural, r$ci_upper_natural)
}


sit_term <- function(tab, term) {
  if (is.null(tab)) return(c(NA_real_, NA_real_, NA_real_))
  r <- tab[tab$term == term, ]
  if (nrow(r) != 1) return(c(NA_real_, NA_real_, NA_real_))
  c(r$estimate, r$ci_lower, r$ci_upper)
}

sitar_lookup <- function(name, source) {
  switch(source,
         fixed = sit_term(sit$fixed_ci, name),
         re    = sit_term(sit$re_ci, name),
         sigma = sit_term(sit$sigma_ci, name),
         none  = c(NA_real_, NA_real_, NA_real_))
}


param_map <- tribble(
  ~quantity,                  ~tmb_name,             ~sitar_name,  ~sitar_source,
  "Intercept b0",             "beta_intercept",      "a",          "fixed",
  "Sex (male) intercept b1",  "beta_intercept_sex",  "a.sex",      "fixed",
  "Sex (male) scale b2",      "beta_amplitude_sex",  NA,           "none",
  "GA, time shift b3",        "beta_shift_ga",       "b.ga",       "fixed",
  "sd(b1)",                   "log_sd_b_intercept",  "sd(a)",      "re",
  "sd(b2)",                   "log_sd_b_shift",      "sd(b)",      "re",
  "rho",                      "transf_rho",          "cor(a,b)",   "re",
  "sigma",                    "log_sigma",           "sigma",      "sigma"
)

param_table <- bind_rows(lapply(seq_len(nrow(param_map)), function(i) {
  m <- param_map[i, ]
  a <- tmb_term(tmb$fixed, m$tmb_name)
  b <- sitar_lookup(m$sitar_name, m$sitar_source)
  tibble(quantity = m$quantity,
         tmb_est = a[1], tmb_lo = a[2], tmb_hi = a[3],
         sit_est = b[1], sit_lo = b[2], sit_hi = b[3])
}))


param_table <- bind_rows(
  param_table,
  tibble(quantity = "Spline df",
         tmb_est = edf_tmb, tmb_lo = NA_real_, tmb_hi = NA_real_,
         sit_est = as.numeric(sit$best_df), sit_lo = NA_real_, sit_hi = NA_real_)
)

# ==============================================================================
# 9. Runtime
# ==============================================================================

sit_setup <- if (is.null(sit$time_df_search)) NA_real_ else sit$time_df_search

runtime <- tibble(
  method = c("snmmTMB", "sitar"),
  fit_seconds = c(tmb$time_fit, sit$time_fit),
  setup_seconds = c(tmb$time_starting_values, sit_setup),
  inference_seconds = c(0, sit$time_bootstrap),
  total_seconds = c(
    tmb$time_fit + tmb$time_starting_values,
    sit$time_fit + (if (is.na(sit_setup)) 0 else sit_setup) + sit$time_bootstrap
  )
)

# ==============================================================================
# 10. Summary
# ==============================================================================

#' estimate (lower, upper), or the estimate alone when there is no interval
cell <- function(est, lo, hi, digits = 2) {
  if (is.na(est)) return("--")
  if (is.na(lo) || is.na(hi)) return(formatC(est, format = "f", digits = digits))
  sprintf("%s (%s, %s)",
          formatC(est, format = "f", digits = digits),
          formatC(lo,  format = "f", digits = digits),
          formatC(hi,  format = "f", digits = digits))
}

tmb_cells <- mapply(cell, param_table$tmb_est, param_table$tmb_lo, param_table$tmb_hi)
sit_cells <- mapply(cell, param_table$sit_est, param_table$sit_lo, param_table$sit_hi)

con <- file(file.path(results_dir, "smocc_comparison_summary.txt"), open = "wt")
w <- function(...) writeLines(sprintf(...), con)

wq <- max(nchar(param_table$quantity), nchar("Quantity"))
w1 <- max(nchar(tmb_cells), nchar("snmmTMB"))
w2 <- max(nchar(sit_cells), nchar("sitar"))
fmt <- paste0("%-", wq, "s  %-", w1, "s  %-", w2, "s")

w("SMOCC application: snmmTMB vs sitar")
w("%d observations on %d subjects", tmb_all$config$nobs, tmb_all$config$nGroup)
w("")
w(fmt, "Quantity", "snmmTMB", "sitar")
w(fmt, strrep("-", wq), strrep("-", w1), strrep("-", w2))
for (i in seq_len(nrow(param_table))) {
  w(fmt, param_table$quantity[i], tmb_cells[i], sit_cells[i])
}
w("")
w("Runtime (seconds)")
w("%-9s %10s %10s %10s %10s", "", "fit", "setup", "inference", "total")
for (i in seq_len(nrow(runtime))) {
  w("%-9s %10.2f %10.2f %10.2f %10.2f",
    runtime$method[i], runtime$fit_seconds[i], runtime$setup_seconds[i],
    runtime$inference_seconds[i], runtime$total_seconds[i])
}
close(con)

writeLines(readLines(file.path(results_dir, "smocc_comparison_summary.txt")))

cat("\nWrote:\n",
    "  Figure5.pdf\n",
    "  ", file.path(results_dir, "smocc_comparison_summary.txt"), "\n", sep = "")
