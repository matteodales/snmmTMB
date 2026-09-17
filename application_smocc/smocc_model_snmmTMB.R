# ------------------------------------------------------------------------------
# smocc_model_snmmTMB.R
#
# SMOCC height application of Section 4 fitted with snmmTMB:
#
#   hgt_ij = beta_intercept + beta_intercept_sex * sex_i + b_intercept_i
#            + exp(beta_amplitude_sex * sex_i)
#              * f(age_ij + ga_i * beta_shift_ga + b_shift_i) + eps_ij
#
#   ( b_intercept_i , b_shift_i ) ~ N2( 0, Sigma ),
#   Sigma = [ sd1^2              rho*sd1*sd2 ]
#           [ rho*sd1*sd2        sd2^2       ],   rho = tanh(transf_rho)
#
#
#   Rscript smocc_model_snmmTMB.R
#
# Reads  data/smocc_200.csv
# Writes results/smocc_snmmTMB_results.RDS
# ------------------------------------------------------------------------------

# ==============================================================================
# 0. Configuration
# ==============================================================================

## Subjects whose individual curves are reported (same set as the manuscript).
subject_ids <- c(1, 2, 3, 4)

## Fixed-grid centering
center_grid_range <- c(0, 1)
center_grid_size <- 1000

ngrid <- 30           # grid points for reported curves
Nsims_band <- 10000   # draws for the simultaneous-band critical value
set.seed(2026)

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

results_file <- file.path(results_dir, "smocc_snmmTMB_results.RDS")
src_dir  <- "../src"
dll_name <- "snmmTMB_likelihood_smocc_application"

# ==============================================================================
# 1. Packages and data
# ==============================================================================

library(TMB)
library(MASS)      # mvrnorm()
library(dplyr)
library(tibble)
library(nlme)
library(splines)   # splineDesign() for the fixed-grid centering constants

smocc_200 <- read.csv("data/smocc_200.csv", stringsAsFactors = FALSE)

smocc_200 <- smocc_200[rowSums(is.na(smocc_200)) == 0, ]
smocc_200$sex <- as.integer(as.factor(smocc_200$sex)) - 1
smocc_200$id <- as.factor(smocc_200$id)
levels(smocc_200$id) <- as.character(seq_len(nlevels(smocc_200$id)))
smocc_200$bw <- smocc_200$bw / 1000
smocc_200$ga <- (smocc_200$ga - 40)
smocc_200$age <- smocc_200$age * 52

smocc_200 <- smocc_200[smocc_200$age < 130, ]
smocc_200 <- smocc_200[smocc_200$hgt > 40, ]

smocc_200 <- smocc_200 %>%
  mutate(id = as.integer(factor(id, levels = unique(id))))
smocc_200$id <- as.factor(smocc_200$id)

subjects_df <- smocc_200 %>%
  mutate(id_int = as.integer(id)) %>%
  distinct(id_int, sex, ga) %>%
  arrange(id_int)

nGroup <- max(as.integer(smocc_200$id))
stopifnot(nrow(subjects_df) == nGroup)

# ==============================================================================
# 2. Knot vector and penalty matrix
# ==============================================================================

degree <- 3
K <- 15
Kint <- K - 1
n_internal <- K - (degree + 1)

knot_lowlim <- 0
knot_uplim <- 1
internal_knots <- seq(knot_lowlim, knot_uplim, length.out = n_internal + 2)

dx_left  <- diff(internal_knots)[1]
dx_right <- diff(internal_knots)[length(internal_knots) - 1]

knot_vec <- c(
  internal_knots[1] - dx_left * (degree:1),
  internal_knots,
  internal_knots[length(internal_knots)] + dx_right * (1:degree)
)

D_pen <- diff(diag(K), differences = 2)[, 1:(K - 1)]
P <- t(D_pen) %*% D_pen
P <- P * qr(P)$rank / sum(diag(P))

ev <- eigen(P, symmetric = TRUE)
pos_idx <- which(ev$values > 1e-12)
zero_idx <- which(ev$values <= 1e-12)

Upos <- if (length(pos_idx) > 0) ev$vectors[, pos_idx, drop = FALSE] else matrix(0, nrow = Kint, ncol = 0)
dpos <- if (length(pos_idx) > 0) ev$values[pos_idx] else numeric(0)
U0   <- if (length(zero_idx) > 0) ev$vectors[, zero_idx, drop = FALSE] else matrix(0, nrow = Kint, ncol = 0)

x_grid <- seq(min(smocc_200$age), max(smocc_200$age), length.out = ngrid)

# ==============================================================================
# 3. Fixed-grid centering constants
# ==============================================================================

edges <- seq(center_grid_range[1], center_grid_range[2],
             length.out = center_grid_size + 1)
v_grid <- (head(edges, -1) + tail(edges, -1)) / 2
Bmean <- colMeans(splineDesign(knot_vec, v_grid, ord = degree + 1,
                               outer.ok = TRUE)[, 1:Kint, drop = FALSE])

# ==============================================================================
# 4. Compile
# ==============================================================================

for (f in c("starting_points", dll_name)) {
  dll <- file.path(src_dir, f)
  if (!file.exists(dynlib(dll))) TMB::compile(paste0(dll, ".cpp"))
  dyn.load(dynlib(dll))
}

# ==============================================================================
# 5. Starting values
# ==============================================================================

t_start0 <- Sys.time()

Data0 <- list(
  y = as.numeric(smocc_200$hgt),
  x = (as.numeric(smocc_200$age) - min(smocc_200$age)) / (max(smocc_200$age) - min(smocc_200$age)),
  knots = as.numeric(knot_vec),
  degree = as.integer(degree),
  K = as.integer(K),
  Upos = as.matrix(Upos),
  U0 = as.matrix(U0),
  dpos = as.numeric(dpos),
  spline_ci = as.integer(1),
  x_grid = as.numeric(x_grid)
)

Params0 <- list(
  beta1 = as.numeric(0),
  log_sigma = log(1),
  vpos = rep(0, length(dpos)),
  gamma0 = rep(0, ncol(U0)),
  log_lambda = log(1)
)

obj0 <- MakeADFun(data = Data0, parameters = Params0,
                  random = if (length(Params0$vpos) > 0) "vpos" else NULL,
                  DLL = "starting_points", silent = TRUE)
opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr, control = list(eval.max = 1e4, iter.max = 1e4))
obj0$par <- opt0$par
rep0 <- sdreport(obj0, getJointPrecision = TRUE)
parList0 <- obj0$env$parList()

c_start <- rep0$value[grep("c", names(rep0$value))]
m_start <- as.numeric(obj0$report()$m)

spline_start <- function(t) {
  sweep(
    as.matrix(splineDesign(
      knot_vec,
      (t - min(smocc_200$age)) / (max(smocc_200$age) - min(smocc_200$age)),
      4, outer.ok = TRUE
    ))[, 1:Kint, drop = FALSE],
    2, m_start, FUN = "-"
  ) %*% c_start
}

nlmeobj <- nlme(
  hgt ~ b1 + exp(b2) * spline_start(age + b3),
  fixed = list(b1 ~ 1 + sex, b2 ~ -1 + sex, b3 ~ -1 + ga),
  random = pdLogChol(b1 + b3 ~ 1),
  data = smocc_200,
  groups = ~id,
  verbose = FALSE,
  start = c(mean(smocc_200$hgt), rep(0, 3)),
  control = list(maxIter = 25)
)

time_starting_values <- as.numeric(difftime(Sys.time(), t_start0, units = "secs"))
message(sprintf("Starting values took %.1f s", time_starting_values))

vc_start <- VarCorr(nlmeobj)
sd_b1_start <- as.numeric(vc_start[1, "StdDev"])
sd_b3_start <- as.numeric(vc_start[2, "StdDev"])
sigma_start <- as.numeric(vc_start[3, "StdDev"])

corr_start <- if ("Corr" %in% colnames(vc_start)) {
  suppressWarnings(as.numeric(vc_start[2, "Corr"]))
} else NA_real_
if (is.na(corr_start)) corr_start <- 0

corr_start <- max(min(corr_start, 0.95), -0.95)

starting_values <- list(
  beta_intercept = as.numeric(nlmeobj$coefficients$fixed["b1.(Intercept)"]),
  beta_intercept_sex = as.numeric(nlmeobj$coefficients$fixed["b1.sex"]),
  beta_amplitude_sex = as.numeric(nlmeobj$coefficients$fixed["b2.sex"]),
  beta_shift_ga = as.numeric(nlmeobj$coefficients$fixed["b3.ga"]),
  b_intercept = as.numeric(nlmeobj$coefficients$random$id[, "b1.(Intercept)"]),
  b_shift = as.numeric(nlmeobj$coefficients$random$id[, "b3.(Intercept)"]),
  log_sd_b_intercept = log(sd_b1_start),
  log_sd_b_shift = log(sd_b3_start),
  transf_rho = atanh(corr_start),
  log_sigma = log(sigma_start),
  vpos = as.numeric(parList0$vpos),
  gamma0 = as.numeric(parList0$gamma0),
  log_lambda = as.numeric(parList0$log_lambda)
)

# ==============================================================================
# 6. Fit
# ==============================================================================

subj_flag_vec <- rep(0L, nGroup)
subj_flag_vec[subject_ids] <- 1L

Data <- list(
  y = as.numeric(smocc_200$hgt),
  age = as.numeric(smocc_200$age),
  age_max = max(smocc_200$age),
  age_min = min(smocc_200$age),
  ga_max = max(smocc_200$ga),
  ga_min = min(smocc_200$ga),
  sex = as.numeric(smocc_200$sex),
  ga = as.numeric(smocc_200$ga),
  sex_subj = as.numeric(subjects_df$sex),
  ga_subj = as.numeric(subjects_df$ga),
  group = as.integer(smocc_200$id),
  nGroup = as.integer(nGroup),
  subj_flag = as.integer(subj_flag_vec),
  knots = as.numeric(knot_vec),
  degree = as.integer(degree),
  K = as.integer(K),
  Bmean = as.numeric(Bmean),
  Upos = as.matrix(Upos),
  U0 = as.matrix(U0),
  dpos = as.numeric(dpos),
  spline_ci = as.integer(1),
  age_grid = as.numeric(x_grid)
)

randoms <- c("b_intercept", "b_shift", "vpos")

message("Fitting ...")

t_fit0 <- Sys.time()
obj_un <- MakeADFun(data = Data, parameters = starting_values, random = randoms,
                    DLL = dll_name, silent = TRUE)
opt_un <- nlminb(obj_un$par, obj_un$fn, obj_un$gr,
                 control = list(eval.max = 1e4, iter.max = 1e4))
obj_un$par <- opt_un$par
sdr_un <- sdreport(obj_un, getJointPrecision = TRUE)
time_fit <- as.numeric(difftime(Sys.time(), t_fit0, units = "secs"))

message(sprintf("  convergence = %d, nll = %.4f, %.1f s",
                opt_un$convergence, opt_un$objective, time_fit))

rho_hat <- as.numeric(obj_un$report(obj_un$env$last.par.best)$rho)

# ==============================================================================
# 7. Parameter estimates with Wald confidence intervals
# ==============================================================================

z975 <- qnorm(0.975)
s <- summary(sdr_un)

keep <- c("beta_intercept", "beta_intercept_sex", "beta_amplitude_sex",
          "beta_shift_ga", "log_sd_b_intercept", "log_sd_b_shift",
          "transf_rho", "log_sigma", "log_lambda", "rho")
keep <- keep[keep %in% rownames(s)]
est <- s[keep, "Estimate"]
se  <- s[keep, "Std. Error"]

fixed_un <- tibble(
  term = keep,
  estimate = as.numeric(est),
  se = as.numeric(se),
  ci_lower = as.numeric(est) - z975 * as.numeric(se),
  ci_upper = as.numeric(est) + z975 * as.numeric(se)
)

is_log  <- startsWith(fixed_un$term, "log_")
is_tanh <- startsWith(fixed_un$term, "transf_")

fixed_un$estimate_natural <- ifelse(is_log, exp(fixed_un$estimate),
                             ifelse(is_tanh, tanh(fixed_un$estimate), fixed_un$estimate))
fixed_un$se_natural <- ifelse(is_log, fixed_un$se * exp(fixed_un$estimate),
                       ifelse(is_tanh, fixed_un$se * (1 - tanh(fixed_un$estimate)^2), fixed_un$se))
fixed_un$ci_lower_natural <- ifelse(is_log, exp(fixed_un$ci_lower),
                             ifelse(is_tanh, tanh(fixed_un$ci_lower), fixed_un$ci_lower))
fixed_un$ci_upper_natural <- ifelse(is_log, exp(fixed_un$ci_upper),
                             ifelse(is_tanh, tanh(fixed_un$ci_upper), fixed_un$ci_upper))

fixed_un$z <- fixed_un$estimate / fixed_un$se
fixed_un$p <- 2 * (1 - pnorm(abs(fixed_un$z)))


random_un <- bind_rows(lapply(c("b_intercept", "b_shift"), function(nm) {
  rows <- rownames(s) == nm
  tibble(
    subject = seq_len(sum(rows)),
    term = nm,
    estimate = as.numeric(s[rows, "Estimate"]),
    se = as.numeric(s[rows, "Std. Error"])
  )
})) %>%
  mutate(ci_lower = estimate - z975 * se,
         ci_upper = estimate + z975 * se)

print(as.data.frame(fixed_un), digits = 4)

## The estimated covariance matrix on the natural scale
sd1_hat <- exp(as.numeric(opt_un$par["log_sd_b_intercept"]))
sd2_hat <- exp(as.numeric(opt_un$par["log_sd_b_shift"]))
Sigma_hat <- matrix(c(sd1_hat^2, rho_hat * sd1_hat * sd2_hat,
                      rho_hat * sd1_hat * sd2_hat, sd2_hat^2), 2, 2,
                    dimnames = list(c("b_intercept", "b_shift"),
                                    c("b_intercept", "b_shift")))
print(Sigma_hat)

# ==============================================================================
# 8. Fitted curves with pointwise and simultaneous 95% bands
# ==============================================================================

curve_bands <- function(mu_hat, Cov, Nsims = Nsims_band) {
  Cov <- (Cov + t(Cov)) / 2
  eig <- eigen(Cov, symmetric = TRUE, only.values = TRUE)
  if (!all(eig$values > 0)) {
    e <- eigen(Cov, symmetric = TRUE)
    e$values[e$values < 1e-10] <- 1e-10
    Cov <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  }
  se_point <- sqrt(diag(Cov))
  sims <- mvrnorm(Nsims, mu = rep(0, ncol(Cov)), Sigma = Cov)
  crit <- as.numeric(quantile(apply(abs(sweep(sims, 2, se_point, FUN = "/")), 1, max), 0.95))

  tibble(
    x = x_grid, y = mu_hat, se = se_point,
    lwrP = mu_hat - z975 * se_point, uprP = mu_hat + z975 * se_point,
    lwrS = mu_hat - crit * se_point, uprS = mu_hat + crit * se_point
  )
}


idx_pop <- grep("^h_grid$", names(sdr_un$value))
population <- curve_bands(sdr_un$value[idx_pop], sdr_un$cov[idx_pop, idx_pop])

idx_subj <- grep("^mu_sel$", names(sdr_un$value))
subjects <- bind_rows(lapply(seq_along(subject_ids), function(k) {
  take <- idx_subj[(ngrid * (k - 1) + 1):(ngrid * k)]
  curve_bands(sdr_un$value[take], sdr_un$cov[take, take]) %>%
    mutate(subject = subject_ids[k], .before = 1)
}))

curves_un <- list(population = population, subjects = subjects)

# ==============================================================================
# 9. Save
# ==============================================================================

rep_un <- obj_un$report(obj_un$env$last.par.best)

bootstrap_inputs <- list(
  opt_par = opt_un$par,
  cov_fixed = sdr_un$cov.fixed,
  parList = obj_un$env$parList(),
  starting_values = starting_values,
  c_hat = as.numeric(rep_un$c),
  a_hat = as.numeric(rep_un$a),
  s_hat = as.numeric(rep_un$s),
  rho_hat = rho_hat,
  Data = Data,
  smocc_200 = smocc_200,
  subjects_df = subjects_df,
  knot_vec = knot_vec, degree = degree, K = K, Kint = Kint,
  Upos = Upos, U0 = U0, dpos = dpos, Bmean = Bmean,
  x_grid = x_grid, nGroup = nGroup
)

out <- list(
  snmmTMB = list(
    label = "snmmTMB (unstructured)",
    centering = "fixed-grid",
    covariance = "unstructured",
    convergence = opt_un$convergence,
    pdHess = isTRUE(sdr_un$pdHess),
    nll = opt_un$objective,
    time_fit = time_fit,
    time_starting_values = time_starting_values,
    time_total = time_fit + time_starting_values,
    rho = rho_hat,
    Sigma = Sigma_hat,
    fixed = fixed_un,
    random = random_un,
    curves = curves_un
  ),
  bootstrap_inputs = bootstrap_inputs,
  x_grid = x_grid,
  subject_ids = subject_ids,
  subjects_df = subjects_df,
  config = list(K = K, degree = degree, knot_vec = knot_vec, Bmean = Bmean,
                center_grid_range = center_grid_range,
                center_grid_size = center_grid_size,
                nGroup = nGroup, nobs = nrow(smocc_200))
)

saveRDS(out, file = results_file)
message("Wrote ", results_file)
