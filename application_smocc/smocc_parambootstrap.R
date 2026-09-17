# ------------------------------------------------------------------------------
# smocc_parambootstrap.R
#
# Parametric bootstrap for the SMOCC application (manuscript Section 4.1), produces Figure 6
#
# `run_bootstrap <- FALSE` skips the computation and redraws the figure from the saved results.
#
# Run after the fit:
#   Rscript smocc_model_snmmTMB.R
#   Rscript smocc_parambootstrap.R
#
# Reads  results/smocc_snmmTMB_results.RDS
# Writes results/smocc_bootstrap_results.RDS
#        Figure6.pdf
# ------------------------------------------------------------------------------

# ==============================================================================
# 0. Configuration
# ==============================================================================

n_boot <- 500                 # bootstrap replicates

run_bootstrap <- FALSE         # FALSE = only redraw the figure from the saved RDS

seed_base <- 20260807

## figure sizes
FIG_WIDTH_IN  <- 6.5
FIG_HEIGHT_IN <- 4

PT_AXIS_TITLE <- 11
PT_AXIS_TEXT  <- 10
PT_LABEL      <- 10
PT_LEGEND     <- 10


label_map <- c(
  beta_intercept      = "beta[0]",
  beta_intercept_sex  = "beta[1]",
  beta_amplitude_sex  = "beta[2]",
  beta_shift_ga       = "beta[3]",
  log_sd_b_intercept  = "log~sigma[b1]",
  log_sd_b_shift      = "log~sigma[b2]",
  transf_rho          = "tilde(rho)",
  log_sigma           = "log~sigma"
)


category_map <- c(
  beta_intercept      = "Regression coefficient",
  beta_intercept_sex  = "Regression coefficient",
  beta_amplitude_sex  = "Regression coefficient",
  beta_shift_ga       = "Regression coefficient",
  log_sd_b_intercept  = "Variance parameter",
  log_sd_b_shift      = "Variance parameter",
  transf_rho          = "Variance parameter",
  log_sigma           = "Variance parameter"
)
shape_values <- c(`Regression coefficient` = 21, `Variance parameter` = 24)

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

results_file <- file.path(results_dir, "smocc_bootstrap_results.RDS")
inputs_file <- file.path(results_dir, "smocc_snmmTMB_results.RDS")
src_dir  <- "../src"
dll_name <- "snmmTMB_likelihood_smocc_application"

# ==============================================================================
# 1. Packages and inputs
# ==============================================================================

library(TMB)
library(splines)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(ggrepel)

if (!file.exists(inputs_file)) {
  stop("missing ", inputs_file,
       " - run smocc_model_snmmTMB.R first")
}

fit_all <- readRDS(inputs_file)
bi <- fit_all$bootstrap_inputs
if (is.null(bi)) {
  stop(inputs_file, " has no `bootstrap_inputs`; re-run the fitting script")
}

dll_path <- file.path(src_dir, dll_name)
if (!file.exists(dynlib(dll_path))) TMB::compile(paste0(dll_path, ".cpp"))
dyn.load(dynlib(dll_path))


#' Make TMB parameter names unique, by appending an index to the elements of
#' vector parameters
#' @param nm character vector of parameter names, as in names(obj$par)
#' @return character vector of the same length, with unique entries
make_unique_par_names <- function(nm) {
  idx <- ave(seq_along(nm), nm, FUN = seq_along)
  ifelse(nm %in% nm[duplicated(nm)], paste0(nm, "_", idx), nm)
}

theta_hat <- bi$opt_par
param_names <- make_unique_par_names(names(theta_hat))
names(theta_hat) <- param_names
se_hat <- sqrt(diag(bi$cov_fixed))
names(se_hat) <- param_names
n_par <- length(param_names)

smocc_200 <- bi$smocc_200
subjects_df <- bi$subjects_df
group_idx <- as.integer(smocc_200$id)
n_obs <- nrow(smocc_200)
nGroup <- bi$nGroup

## per-observation subject covariates
sex_obs <- subjects_df$sex[group_idx]
ga_obs <- subjects_df$ga[group_idx]

beta_intercept     <- as.numeric(theta_hat["beta_intercept"])
beta_intercept_sex <- as.numeric(theta_hat["beta_intercept_sex"])
beta_amplitude_sex <- as.numeric(theta_hat["beta_amplitude_sex"])
beta_shift_ga      <- as.numeric(theta_hat["beta_shift_ga"])
sd_b_intercept     <- exp(as.numeric(theta_hat["log_sd_b_intercept"]))
sd_b_shift         <- exp(as.numeric(theta_hat["log_sd_b_shift"]))
rho_hat            <- tanh(as.numeric(theta_hat["transf_rho"]))
sigma_resid        <- exp(as.numeric(theta_hat["log_sigma"]))

knot_left <- bi$knot_vec[1]
knot_right <- bi$knot_vec[length(bi$knot_vec)]

sqrt_1m_rho2 <- sqrt(1 - rho_hat^2)

#' Simulate one bootstrap response vector from the fitted model
#'
#' @return list(y, b_intercept, b_shift)
simulate_bootstrap_y <- function() {

  z1 <- rnorm(nGroup)
  z2 <- rnorm(nGroup)
  b_intercept <- sd_b_intercept * z1
  b_shift <- sd_b_shift * (rho_hat * z1 + sqrt_1m_rho2 * z2)

  u <- smocc_200$age + ga_obs * beta_shift_ga + b_shift[group_idx]
  v <- (u - bi$a_hat) / bi$s_hat
  v <- pmin(pmax(v, knot_left), knot_right - 1e-8 * (knot_right - knot_left))

  B <- splineDesign(bi$knot_vec, v, ord = bi$degree + 1,
                    outer.ok = TRUE)[, 1:bi$Kint, drop = FALSE]
  h <- as.numeric(sweep(B, 2, bi$Bmean, FUN = "-") %*% bi$c_hat)

  mu <- beta_intercept + beta_intercept_sex * sex_obs + b_intercept[group_idx] +
    exp(beta_amplitude_sex * sex_obs) * h

  list(y = mu + rnorm(n_obs, 0, sigma_resid),
       b_intercept = b_intercept, b_shift = b_shift)
}

# ==============================================================================
# 3. Bootstrap loop
# ==============================================================================

Data_boot <- bi$Data
Data_boot$spline_ci <- as.integer(0)
Data_boot$subj_flag <- as.integer(rep(0, nGroup))
Params_start <- bi$parList
randoms <- c("b_intercept", "b_shift", "vpos")

if (run_bootstrap) {

  if (file.exists(results_file)) {
    store <- readRDS(results_file)
    message("Resuming: ", sum(store$done), " / ", n_boot, " replicates already on disk.")
  } else {
    store <- list(
      estimates = matrix(NA_real_, n_boot, n_par, dimnames = list(NULL, param_names)),
      ses       = matrix(NA_real_, n_boot, n_par, dimnames = list(NULL, param_names)),
      convergence = rep(NA_integer_, n_boot),
      pdHess = rep(NA, n_boot),
      time = rep(NA_real_, n_boot),
      done = rep(FALSE, n_boot),
      theta_hat = theta_hat,
      se_hat = se_hat,
      rho_true = rho_hat,
      config = list(n_boot = n_boot, seed_base = seed_base)
    )
  }

  if (nrow(store$estimates) < n_boot) {
    pad <- n_boot - nrow(store$estimates)
    store$estimates <- rbind(store$estimates, matrix(NA_real_, pad, n_par, dimnames = list(NULL, param_names)))
    store$ses <- rbind(store$ses, matrix(NA_real_, pad, n_par, dimnames = list(NULL, param_names)))
    store$convergence <- c(store$convergence, rep(NA_integer_, pad))
    store$pdHess <- c(store$pdHess, rep(NA, pad))
    store$time <- c(store$time, rep(NA_real_, pad))
    store$done <- c(store$done, rep(FALSE, pad))
  }

  for (b in seq_len(n_boot)) {

    if (isTRUE(store$done[b])) next

    ## one seed per replicate
    set.seed(seed_base + b)

    cat(sprintf("[%s] bootstrap replicate %d / %d\n",
                format(Sys.time(), "%H:%M:%S"), b, n_boot))

    sim <- simulate_bootstrap_y()
    Data_b <- Data_boot
    Data_b$y <- as.numeric(sim$y)

    t0 <- Sys.time()

    obj <- tryCatch(
      MakeADFun(data = Data_b, parameters = Params_start, random = randoms,
                DLL = dll_name, silent = TRUE),
      error = function(e) { message("  MakeADFun failed: ", e$message); NULL }
    )
    if (is.null(obj)) { store$done[b] <- TRUE; saveRDS(store, results_file); next }

    opt <- tryCatch(
      nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e4, iter.max = 1e4)),
      error = function(e) { message("  nlminb failed: ", e$message); NULL }
    )
    if (is.null(opt)) { store$done[b] <- TRUE; saveRDS(store, results_file); next }

    obj$par <- opt$par
    
    rep_b <- tryCatch(sdreport(obj, skip.delta.method = TRUE),
                      error = function(e) { message("  sdreport failed: ", e$message); NULL })
    if (is.null(rep_b)) { store$done[b] <- TRUE; saveRDS(store, results_file); next }

    store$estimates[b, ] <- as.numeric(opt$par)
    store$ses[b, ] <- as.numeric(sqrt(diag(rep_b$cov.fixed)))
    store$convergence[b] <- opt$convergence
    store$pdHess[b] <- isTRUE(rep_b$pdHess)
    store$time[b] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    store$done[b] <- TRUE

    saveRDS(store, results_file)
  }
}

# ==============================================================================
# 4. Usable replicates
# ==============================================================================

if (!file.exists(results_file)) {
  stop("missing ", results_file,
       " - set run_bootstrap <- TRUE to produce it")
}
store <- readRDS(results_file)

ok <- store$done &
  !is.na(store$convergence) & store$convergence == 0 &
  !is.na(store$pdHess) & store$pdHess &
  rowSums(!is.finite(store$estimates)) == 0 &
  rowSums(!is.finite(store$ses)) == 0

est <- store$estimates[ok, , drop = FALSE]
ses <- store$ses[ok, , drop = FALSE]
B_ok <- nrow(est)

message(sprintf("Usable replicates: %d / %d (%.1f%%). Median fit time %.1f s.",
                B_ok, sum(store$done), 100 * B_ok / max(1, sum(store$done)),
                median(store$time, na.rm = TRUE)))
if (B_ok < 30) warning("fewer than 30 usable replicates; the diagnostics are noisy")


# ==============================================================================
# 5. Figure 6
# ==============================================================================

plot_params <- param_names[1:8]
stopifnot(identical(plot_params, names(label_map)))

plot_df <- tibble(
  parameter = plot_params,
  theta_hat = as.numeric(store$theta_hat[plot_params]),
  bootstrap_mean = colMeans(est[, plot_params, drop = FALSE], na.rm = TRUE),
  bootstrap_sd = apply(est[, plot_params, drop = FALSE], 2, sd, na.rm = TRUE),
  mean_asymptotic_se = colMeans(ses[, plot_params, drop = FALSE], na.rm = TRUE)
) %>%
  mutate(
    label = label_map[parameter],
    category = factor(category_map[parameter],
                      levels = c("Regression coefficient",
                                 "Variance parameter"))
  )

plot_df <- plot_df %>%
  mutate(
    theta_hat = ifelse(parameter == "transf_rho", -theta_hat, theta_hat),
    bootstrap_mean = ifelse(parameter == "transf_rho", -bootstrap_mean, bootstrap_mean)
  )

print(as.data.frame(plot_df), digits = 4)

# ==============================================================================
# 6. Shared layers
# ==============================================================================

point_layer <- geom_point(aes(shape = category), colour = "black", fill = "grey30",
                          size = 3.2, stroke = 0.7, alpha = 0.9)
abline_layer <- geom_abline(intercept = 0, slope = 1, colour = "grey55",
                            linetype = "dashed", linewidth = 0.7)

fig_theme <- theme_bw(base_size = PT_AXIS_TEXT) +
  theme(
    axis.title   = element_text(size = PT_AXIS_TITLE),
    axis.text    = element_text(size = PT_AXIS_TEXT, colour = "black"),
    legend.text  = element_text(size = PT_LEGEND),
    legend.title = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title    = element_blank(),
    plot.subtitle = element_blank()
  )

## labels are nudged off the diagonal
NUDGE_DELTA <- 0.35
nudge_side <- rep(c(1, -1), length.out = nrow(plot_df))

# ==============================================================================
# 7. Panel (a)
# ==============================================================================

pseudo_log <- scales::pseudo_log_trans(sigma = 0.1, base = exp(1))
pseudo_breaks <- c(-1, 0, 1, 5, 20, 70)

label_layer_a <- geom_text_repel(
  aes(label = label), parse = TRUE, size = PT_LABEL / .pt,
  nudge_x = -nudge_side * NUDGE_DELTA, nudge_y = nudge_side * NUDGE_DELTA,
  min.segment.length = 0, seed = 1, max.overlaps = Inf,
  box.padding = 0.2, point.padding = 0.2, force = 0.5, force_pull = 2
)

p_est <- ggplot(plot_df, aes(theta_hat, bootstrap_mean)) +
  abline_layer + point_layer + label_layer_a +
  scale_x_continuous(trans = pseudo_log, breaks = pseudo_breaks,
                     expand = expansion(mult = 0.14)) +
  scale_y_continuous(trans = pseudo_log, breaks = pseudo_breaks,
                     expand = expansion(mult = 0.14)) +
  scale_shape_manual(values = shape_values, drop = FALSE) +
  labs(x = "Original model estimate", y = "Bootstrap mean") +
  fig_theme

# ==============================================================================
# 8. Panel (b)
# ==============================================================================

label_layer_b <- geom_text_repel(
  aes(label = label), parse = TRUE, size = PT_LABEL / .pt,
  nudge_x = -nudge_side * 0.15, nudge_y = nudge_side * 0.15,
  min.segment.length = 0, seed = 1, max.overlaps = Inf,
  box.padding = 0.2, point.padding = 0.2, force = 0.5, force_pull = 2
)

p_se <- ggplot(plot_df, aes(mean_asymptotic_se, bootstrap_sd)) +
  abline_layer + point_layer + label_layer_b +
  scale_x_log10(expand = expansion(mult = 0.16)) +
  scale_y_log10(expand = expansion(mult = 0.16)) +
  scale_shape_manual(values = shape_values, drop = FALSE) +
  labs(x = "Average asymptotic SE", y = "Bootstrap SE") +
  fig_theme

# ==============================================================================
# 9. Assemble
# ==============================================================================

figure6 <- (p_est + p_se) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        plot.tag = element_text(size = 12, face = "bold"))

ggsave("Figure6.pdf", figure6, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN)

cat("\nWrote Figure6.pdf\n")
