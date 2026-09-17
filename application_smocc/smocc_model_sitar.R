# ------------------------------------------------------------------------------
# smocc_model_sitar.R
#
# SMOCC height application of Section 4 fitted with sitar (Cole et al. 2010)
#
# sitar fits  mu_ij = a_i + f(age_ij - b_i),  with a_i a subject size shift and
# b_i a subject tempo shift, and f a natural cubic spline whose degrees of freedom are fixed in advance. 
# This script chooses df by BIC over a grid, then refits at the chosen value and
# bootstraps the fitted curve, since sitar supplies no band of its own.
#
#
#   Rscript smocc_model_sitar.R
#
# Reads  data/smocc_200.csv
# Writes results/smocc_sitar_results.RDS
# ------------------------------------------------------------------------------

# ==============================================================================
# 0. Configuration
# ==============================================================================

## Subjects whose individual curves are reported
subject_ids <- c(1, 2, 3, 4)

## Candidate degrees of freedom for the mean curve
spline_df_grid <- 3:20

ngrid <- 30            # grid points for reported curves
n_boot <- 200          # parametric-bootstrap replicates
boot_seed <- 0

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

results_file <- file.path(results_dir, "smocc_sitar_results.RDS")

# ==============================================================================
# 1. Packages and data
# ==============================================================================

library(sitar)
library(nlme)
library(MASS)      # mvrnorm()
library(dplyr)
library(tibble)


smocc_200 <- read.csv("data/smocc_200.csv", stringsAsFactors = FALSE)

smocc_200 <- smocc_200[rowSums(is.na(smocc_200)) == 0, ]
smocc_200$sex <- as.integer(as.factor(smocc_200$sex)) - 1
smocc_200$id <- as.factor(smocc_200$id)
levels(smocc_200$id) <- as.character(seq_len(nlevels(smocc_200$id)))
smocc_200$bw <- smocc_200$bw / 1000
smocc_200$ga <- -(smocc_200$ga - 40) # the sign is flipped in this because of the way sitar works, to match the sign of the coefficients and correlations
smocc_200$age <- smocc_200$age * 52

smocc_200 <- smocc_200[smocc_200$age < 130, ]
smocc_200 <- smocc_200[smocc_200$hgt > 40, ]

smocc_200 <- smocc_200 %>%
  mutate(id = as.integer(factor(id, levels = unique(id))))
smocc_200$id <- as.factor(smocc_200$id)

nGroup <- nlevels(smocc_200$id)
x_grid <- seq(min(smocc_200$age), max(smocc_200$age), length.out = ngrid)

# ==============================================================================
# 2. Choose the mean curve's df by BIC
# ==============================================================================

sitar_bic_table <- tibble(spline_df = spline_df_grid, BIC = NA_real_)

t_df0 <- Sys.time()
for (i in seq_along(spline_df_grid)) {
  spline_df <- spline_df_grid[i]
  fit_i <- tryCatch(
    sitar(x = age, y = hgt, id = id, data = smocc_200, df = spline_df,
          a.formula = ~ 1 + sex, b.formula = ~ -1 + ga, c.formula = ~ 1,
          random = "a+b"),
    error = function(e) NULL
  )
  if (!is.null(fit_i)) sitar_bic_table$BIC[i] <- BIC(fit_i)
}
time_df_search <- as.numeric(difftime(Sys.time(), t_df0, units = "secs"))

print(sitar_bic_table)
if (all(is.na(sitar_bic_table$BIC))) stop("every sitar df candidate failed to fit")
sitar_best_df <- sitar_bic_table$spline_df[which.min(sitar_bic_table$BIC)]
message(sprintf("Selected sitar df = %d by BIC (search took %.1f s)",
                sitar_best_df, time_df_search))

# ==============================================================================
# 3. Final fit
# ==============================================================================

spline_df <- sitar_best_df

t0 <- Sys.time()
sitar_fit <- sitar(x = age, y = hgt, id = id, data = smocc_200, df = spline_df,
                   a.formula = ~ 1 + sex, b.formula = ~ -1 + ga, c.formula = ~ 1,
                   random = "a+b")
time_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(sprintf("sitar final fit: %.2f s", time_fit))

print(summary(sitar_fit))
print(VarCorr(sitar_fit))

# ==============================================================================
# 4. Parameter estimates with confidence intervals
# ==============================================================================

sitar_tTable <- as.data.frame(summary(sitar_fit)$tTable)
sitar_tTable <- rownames_to_column(sitar_tTable, "term")
print(sitar_tTable)

sitar_intervals <- tryCatch(intervals(sitar_fit), error = function(e) {
  message("intervals() failed: ", e$message)
  NULL
})

if (!is.null(sitar_intervals)) {
  sitar_fixed_ci <- as.data.frame(sitar_intervals$fixed) %>%
    rownames_to_column("term") %>%
    rename(ci_lower = lower, estimate = `est.`, ci_upper = upper) %>%
    as_tibble()
  sitar_re_ci <- as.data.frame(sitar_intervals$reStruct$id) %>%
    rownames_to_column("term") %>%
    rename(ci_lower = lower, estimate = `est.`, ci_upper = upper) %>%
    as_tibble()
  sitar_sigma_ci <- tibble(
    term = "sigma",
    ci_lower = as.numeric(sitar_intervals$sigma[["lower"]]),
    estimate = as.numeric(sitar_intervals$sigma[["est."]]),
    ci_upper = as.numeric(sitar_intervals$sigma[["upper"]])
  )
} else {
  z975 <- qnorm(0.975)
  sitar_fixed_ci <- sitar_tTable %>%
    transmute(term,
              estimate = Value,
              ci_lower = Value - z975 * Std.Error,
              ci_upper = Value + z975 * Std.Error) %>%
    as_tibble()
  sitar_re_ci <- NULL
  sitar_sigma_ci <- NULL
}

print(sitar_fixed_ci)
print(sitar_re_ci)

## Random-effect covariance
vc <- VarCorr(sitar_fit)
sd_a <- as.numeric(vc["a", "StdDev"])
sd_b <- as.numeric(vc["b", "StdDev"])
corr_ab <- suppressWarnings(as.numeric(vc["b", "Corr"]))
if (is.na(corr_ab)) corr_ab <- 0
D_hat <- matrix(c(sd_a^2, corr_ab * sd_a * sd_b,
                  corr_ab * sd_a * sd_b, sd_b^2), 2, 2)
sigma_hat <- sitar_fit$sigma

message(sprintf("sd(a) = %.3f, sd(b) = %.3f, cor(a,b) = %.3f, sigma = %.3f",
                sd_a, sd_b, corr_ab, sigma_hat))

## Predicted random effects
sitar_ranef <- as.data.frame(ranef(sitar_fit)) %>%
  rownames_to_column("subject") %>%
  as_tibble()

# ==============================================================================
# 5. Fitted curves (point estimates)
#
# Population curve at sex = 0, ga = 0
# Subject curves with that subject's own sex, ga and ranef
# ==============================================================================

grid_pop <- data.frame(age = x_grid, sex = 0, ga = 0, id = smocc_200$id[1])

subject_covariates <- smocc_200 %>%
  mutate(id_int = as.integer(id)) %>%
  distinct(id_int, sex, ga) %>%
  arrange(id_int) %>%
  filter(id_int %in% subject_ids)

grid_subj_list <- lapply(seq_len(nrow(subject_covariates)), function(k) {
  s <- subject_covariates[k, ]
  data.frame(age = x_grid, sex = s$sex, ga = s$ga,
             id = factor(as.character(s$id_int), levels = levels(smocc_200$id)))
})
names(grid_subj_list) <- as.character(subject_covariates$id_int)

sitar_mu_pop <- as.numeric(suppressMessages(
  predict(sitar_fit, newdata = grid_pop, level = 0)))
sitar_mu_subj <- lapply(grid_subj_list, function(g) {
  as.numeric(suppressMessages(predict(sitar_fit, newdata = g, level = 1)))
})

# ==============================================================================
# 6. Parametric bootstrap for curve standard errors and bands
#
# sitar returns no confidence band for the fitted curve, so one is simulated:
# draw new random effects and residuals from the fitted model, refit, and use
# the spread of the refitted curves as the sampling distribution.
# ==============================================================================

set.seed(boot_seed)

err_pop <- matrix(NA_real_, nrow = n_boot, ncol = ngrid)
err_subj <- array(NA_real_, dim = c(n_boot, length(grid_subj_list), ngrid))
boot_fixed <- matrix(NA_real_, nrow = n_boot, ncol = length(fixef(sitar_fit)),
                     dimnames = list(NULL, names(fixef(sitar_fit))))

t_boot0 <- Sys.time()

for (b in seq_len(n_boot)) {

  # new subject random effects from the fitted covariance
  re_star <- mvrnorm(nGroup, mu = c(0, 0), Sigma = D_hat)
  colnames(re_star) <- c("a", "b")

  # truth for this replicate
  fit_star <- sitar_fit
  fit_star$coefficients$random$id <- re_star

  # subject-level predictions at the observed design points
  mu_star <- tryCatch(as.numeric(suppressMessages(predict(fit_star, level = 1))),
                      error = function(e) NULL)
  if (is.null(mu_star)) next

  # and on the reporting grid
  true_subj_b <- tryCatch(
    lapply(grid_subj_list, function(g) {
      as.numeric(suppressMessages(predict(fit_star, newdata = g, level = 1)))
    }),
    error = function(e) NULL
  )
  if (is.null(true_subj_b)) next

  # simulated response
  data_star <- smocc_200
  data_star$hgt <- mu_star + rnorm(nrow(smocc_200), 0, sigma_hat)

  # refit
  refit <- tryCatch(
    sitar(x = age, y = hgt, id = id, data = data_star, df = spline_df,
          a.formula = ~ 1 + sex, b.formula = ~ -1 + ga, c.formula = ~ 1,
          random = "a+b"),
    error = function(e) NULL
  )
  if (is.null(refit)) next


  est_pop_b <- tryCatch(
    as.numeric(suppressMessages(predict(refit, newdata = grid_pop, level = 0))),
    error = function(e) NULL)
  if (is.null(est_pop_b)) next
  err_pop[b, ] <- est_pop_b - sitar_mu_pop

  for (k in seq_along(grid_subj_list)) {
    est_k <- tryCatch(
      as.numeric(suppressMessages(predict(refit, newdata = grid_subj_list[[k]], level = 1))),
      error = function(e) NULL)
    if (!is.null(est_k)) err_subj[b, k, ] <- est_k - true_subj_b[[k]]
  }

  boot_fixed[b, ] <- fixef(refit)[colnames(boot_fixed)]

  if (b %% 25 == 0) message(sprintf("  bootstrap %d / %d", b, n_boot))
}

time_bootstrap <- as.numeric(difftime(Sys.time(), t_boot0, units = "secs"))
n_boot_ok <- sum(rowSums(is.na(err_pop)) == 0)
message(sprintf("Bootstrap: %d / %d replicates usable, %.1f s",
                n_boot_ok, n_boot, time_bootstrap))

#' Pointwise and simultaneous 95% bands from a matrix of bootstrap errors.
#' @param mu_hat point estimate on the grid
#' @param err    n_boot x ngrid matrix
#' @return tibble with x, y, se and the two sets of bounds
boot_bands <- function(mu_hat, err) {
  err <- err[rowSums(is.na(err)) == 0, , drop = FALSE]   # drop failed replicates
  se_point <- apply(err, 2, sd)                          # pointwise se
  # simultaneous critical value
  crit <- as.numeric(quantile(apply(abs(sweep(err, 2, se_point, FUN = "/")), 1, max),
                              0.95, na.rm = TRUE))
  tibble(
    x = x_grid, y = mu_hat, se = se_point,
    lwrP = mu_hat - 1.96 * se_point, uprP = mu_hat + 1.96 * se_point,
    lwrS = mu_hat - crit * se_point, uprS = mu_hat + crit * se_point
  )
}

sitar_curve_pop <- boot_bands(sitar_mu_pop, err_pop)

sitar_curve_subj <- bind_rows(lapply(seq_along(grid_subj_list), function(k) {
  boot_bands(sitar_mu_subj[[k]], err_subj[, k, , drop = TRUE]) %>%
    mutate(subject = as.integer(names(grid_subj_list)[k]), .before = 1)
}))

# ==============================================================================
# 7. Save
# ==============================================================================

saveRDS(
  list(
    label = "sitar",
    bic_table = sitar_bic_table,
    best_df = sitar_best_df,
    time_df_search = time_df_search,
    time_fit = time_fit,
    time_bootstrap = time_bootstrap,
    n_boot = n_boot,
    n_boot_ok = n_boot_ok,
    tTable = sitar_tTable,
    fixed_ci = sitar_fixed_ci,
    re_ci = sitar_re_ci,
    sigma_ci = sitar_sigma_ci,
    ranef = sitar_ranef,
    D_hat = D_hat,
    corr_ab = corr_ab,
    sigma_hat = sigma_hat,
    boot_fixed = boot_fixed,
    x_grid = x_grid,
    subject_ids = subject_ids,
    curves = list(population = sitar_curve_pop, subjects = sitar_curve_subj)
  ),
  file = results_file
)
message("Wrote ", results_file)
