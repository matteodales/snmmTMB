# ------------------------------------------------------------------------------
# sinecurve_server_simulation_t_residuals_snmmTMB.R
#
# snmmTMB script for the sine-curve simulation with t-distributed residuals of
# Online Appendix B (Supplementary Figure 1)
#
# Identical to ../../../simulation_sinecurve/sinecurve_server_simulation_snmmTMB.R
# except for the residual draw in simulate_dataset(): the errors come from a
# scaled t distribution with err_df = 3 degrees of freedom instead of a normal,
# rescaled so that Var(eps_ij) = sigma^2 as in the Gaussian run. The fitted model
# is unchanged and still assumes Gaussian residuals.
#
# seed_base is the same as in the Gaussian run: simulate_dataset() draws the
# random effects before the residuals, so the true b_i of every replication are
# identical to the Gaussian run and only the errors differ.
#
# This script runs the full study and is meant to run on a computing cluster
# One iteration of the Gaussian model can be fitted and plotted by running
# ../../../simulation_sinecurve/sinecurve_example.R
#
# Usage, one process per job:
#
#   Rscript sinecurve_server_compile.R                       # once, first
#   Rscript sinecurve_server_simulation_t_residuals_snmmTMB.R <job> <njobs>
#
# with <job> in 1..<njobs>; sinecurve_server_run_t_residuals_snmmTMB.sh launches all of the jobs at once.
# Each job writes only its own results/<output_prefix>_job<j>_of<njobs>.RDS
# ------------------------------------------------------------------------------

# ==============================================================================
# 0. Command-line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript sinecurve_server_simulation_t_residuals_snmmTMB.R <job> [<njobs>]")
}

itnum <- as.integer(args[1])
njobs <- if (length(args) >= 2) as.integer(args[2]) else 20L

if (is.na(itnum) || is.na(njobs) || njobs < 1L || itnum < 1L || itnum > njobs) {
  stop("Bad arguments: need 1 <= job <= njobs, got job = ", args[1],
       ", njobs = ", if (length(args) >= 2) args[2] else njobs)
}

# ==============================================================================
# 1. Run configuration
# ==============================================================================

iterations <- 500            # replications per setting
settings_to_run <- 1:16
seed_base <- 20260804          # same as the Gaussian run, so the random effects are paired

err_df <- 3                    # degrees of freedom of the residual t; must be > 2
                               # for the variance to exist

Nsims <- 10000                 # draws for the simultaneous-band critical value

src_dir <- "../../../src"
results_dir <- "results"
output_prefix <- "sinecurve_server_t_residuals_snmmTMB"
save_every <- 1L

# ==============================================================================
# 2. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm() for multivariate-normal random-effect draws
  library(dplyr)    # bind_rows()
  library(tibble)   # tibble() in place of data.frame()
  library(TMB)      # automatic differentiation / Laplace approximation
  library(nlme)     # nlme() for starting values
  library(splines)  # splineDesign() for the starting-value GAM and Bmean
})

# ==============================================================================
# 3. Knot vector and penalty matrix
# ==============================================================================

degree <- 3
K <- 15                 # number of basis functions
Kint <- K - 1
n_internal <- K - (degree + 1)

knot_lowlim <- -1
knot_uplim <- 1
internal_knots <- seq(knot_lowlim, knot_uplim, length.out = n_internal + 2)

if (length(internal_knots) > 1) {
  dx_left  <- diff(internal_knots)[1]
  dx_right <- diff(internal_knots)[length(internal_knots) - 1]
} else {
  dx_left <- dx_right <- 1
}

knot_vec <- c(
  internal_knots[1] - dx_left * (degree:1),
  internal_knots,
  internal_knots[length(internal_knots)] + dx_right * (1:degree)
)

D <- diff(diag(K), differences = 2)[, 1:(K - 1)] # drop one column to enforce the sum-to-zero constraint on the smooth
P <- t(D) %*% D
P <- P * qr(P)$rank / sum(diag(P))

# decompose into null and active space of the penalty
ev <- eigen(P, symmetric = TRUE)
pos_idx <- which(ev$values > 1e-12)
zero_idx <- which(ev$values <= 1e-12)

Upos <- if (length(pos_idx) > 0) ev$vectors[, pos_idx, drop = FALSE] else matrix(0, nrow = K, ncol = 0)
dpos <- if (length(pos_idx) > 0) ev$values[pos_idx] else numeric(0)
U0   <- if (length(zero_idx) > 0) ev$vectors[, zero_idx, drop = FALSE] else matrix(0, nrow = K, ncol = 0)

t_grid <- seq(0, 1, length.out = 50)  # grid the curves are reported on
ngrid <- length(t_grid)

# ==============================================================================
# 4. Fixed-grid centering constants
# ==============================================================================

center_grid_range <- c(-1, 1)
center_grid_size <- 1000

# Column means of the B-spline basis over a fixed grid

fixed_grid_center <- function(rng, ngrid = 1000) {
  edges <- seq(rng[1], rng[2], length.out = ngrid + 1)
  u_grid <- (head(edges, -1) + tail(edges, -1)) / 2
  colMeans(splineDesign(knot_vec, u_grid, ord = degree + 1,
                        outer.ok = TRUE)[, 1:Kint, drop = FALSE])
}

Bmean <- fixed_grid_center(center_grid_range, center_grid_size)

# ==============================================================================
# 5. Coverage helper functions
# ==============================================================================

#' Simultaneous coverage indicator
#'
#' @param x   vector of true values
#' @param upr,lwr  vectors of upper/lower confidence bounds, same length as x
#' @return TRUE if every element of x falls within its band (simultaneous
#'   coverage for one replication), FALSE otherwise.
inCI <- function(x, upr, lwr) {
  all(x >= lwr & x <= upr)
}

#' Pointwise coverage proportion
#'
#' @param x   vector of true values
#' @param upr,lwr  vectors of upper/lower confidence bounds, same length as x
#' @return proportion of elements of x falling within their band (pointwise
#'   coverage for one replication).
meanCI <- function(x, upr, lwr) {
  mean(x >= lwr & x <= upr)
}

# ==============================================================================
# 6. Simulation grid
#
# 16 settings = 2 (n) x 2 (m) x 2 (sigma) x 2 (D, low/high variance)
# ==============================================================================

n_vals <- c(10, 50)
m_vals <- c(10, 20)
sigma_vals <- c(0.4, 1)
D_diag_vals <- list(
  c(0.25, 0.16, 0.04),
  c(1, 0.25, 0.16)
)

param_grid <- expand.grid(
  n = n_vals,
  m = m_vals,
  sigma = sigma_vals,
  D_idx = seq_along(D_diag_vals)
)

alpha_true <- 5

# ==============================================================================
# 7. Load the compiled likelihoods
# ==============================================================================

cpp_start <- file.path(src_dir, "starting_points.cpp")
cpp_model <- file.path(src_dir, "snmmTMB_likelihood_sinecurve_simulation.cpp")

#' dyn.load() a TMB model, refusing to compile it here
#'
#' @param cpp_file path to the .cpp source
#' @return the absolute path to the loaded shared library, invisibly
load_model <- function(cpp_file) {
  dll <- dynlib(sub("\\.cpp$", "", cpp_file))
  if (!file.exists(dll)) {
    stop("Shared library not found: ", dll)
  }
  dyn.load(dll)
  invisible(normalizePath(dll))
}

load_model(cpp_start)
load_model(cpp_model)

model_name <- sub("\\.cpp$", "", basename(cpp_model))
start_name <- sub("\\.cpp$", "", basename(cpp_start))

# ==============================================================================
# 8. One replication
# ==============================================================================
#'
#' @param setting_idx row of param_grid
#' @param iter replication number within the setting
#' @return integer seed
seed_for <- function(setting_idx, iter) {
  seed_base + 977L * as.integer(setting_idx) + as.integer(iter)
}

#' Simulate one sine-curve dataset with t-distributed residuals
#'
#' Data-generating model (Online Appendix B, first variant):
#'
#'   y_ij = alpha + b1_i + 2*exp(b2_i) * sin(2*pi*(t_ij - expit(b3_i))) + eps_ij
#'   b_i = (b1_i, b2_i, b3_i) ~ N(0, sigma^2 * D),  D = diag(D_diag)
#'   eps_ij = sigma * T_ij / sqrt(nu / (nu - 2)),  T_ij ~ t_nu,  nu = err_df
#'
#' so that Var(eps_ij) = sigma^2 as in the Gaussian run. The random effects are
#' drawn first, exactly as in the Gaussian run, so they are identical to it for
#' the same seed.
#'
#' @param n,m subjects and observations per subject
#' @param sigma residual sd; D_diag length-3 random-effect variance multipliers
#' @param alpha population intercept
#' @return list(dat = tibble, b_mat = n x 3 matrix of true random effects,
#'   Sigma_b = the random-effect covariance used)

simulate_dataset <- function(n, m, sigma, D_diag, alpha) {
  Sigma_b <- sigma^2 * diag(D_diag)
  b_mat <- mvrnorm(n = n, mu = rep(0, 3), Sigma = Sigma_b)
  colnames(b_mat) <- c("b1", "b2", "b3")

  dat_list <- vector("list", n)
  for (i in 1:n) {
    b1 <- b_mat[i, "b1"]
    b2 <- b_mat[i, "b2"]
    b3 <- b_mat[i, "b3"]

    t <- seq(0, 1, length.out = m)

    amp <- 2 * exp(b2)
    phase <- exp(b3) / (1 + exp(b3))
    mu_ij <- alpha + b1 + amp * sin(2 * pi * (t - phase))
    y_ij  <- mu_ij + sigma * rt(m, df = err_df) / sqrt(err_df / (err_df - 2))

    dat_list[[i]] <- tibble(
      subject = rep(i, m),
      j = 1:m,
      t = t,
      y = y_ij,
      mu = mu_ij,
      b1 = rep(b1, m),
      b2 = rep(b2, m),
      b3 = rep(b3, m)
    )
  }

  list(dat = bind_rows(dat_list), b_mat = b_mat, Sigma_b = Sigma_b)
}

#' Run one replication
#'
#' @param setting_idx row of param_grid
#' @param iter replication number within the setting
#' @return named list of metrics, or NULL on failure
run_replication <- function(setting_idx, iter) {

  on.exit(if (exists("spline_start", envir = .GlobalEnv, inherits = FALSE)) {
    rm("spline_start", envir = .GlobalEnv)
  }, add = TRUE)

  set.seed(seed_for(setting_idx, iter),
           kind = "Mersenne-Twister", normal.kind = "Inversion")

  n <- param_grid$n[setting_idx]
  m <- param_grid$m[setting_idx]
  sigma <- param_grid$sigma[setting_idx]
  D_diag <- D_diag_vals[[param_grid$D_idx[setting_idx]]]
  subj_flag_vec <- c(1, 1, rep(0, n - 2))

  sim <- simulate_dataset(n, m, sigma, D_diag, alpha_true)
  dat <- sim$dat
  b_mat <- sim$b_mat
  Sigma_b <- sim$Sigma_b

  out <- tryCatch({

    # ------------------------------------------------------------------------
    # Starting values (manuscript Section 2.3.2): fit a simple penalized-spline
    # GAM (no random effects, starting_points.cpp) with the same knot vector
    # and penalty matrix as the final model, then use its fitted curve inside
    # an nlme() fit to get starting values for the fixed/random effects.
    # ------------------------------------------------------------------------
    Params0 <- list(
      beta1 = as.numeric(0), # intercept
      log_sigma = log(1),
      vpos = rep(0, length(dpos)),
      gamma0 = rep(0, ncol(U0)),
      log_lambda = log(1)
    )

    Data0 <- list(
      y = dat$y,
      x = dat$t - 0.5,
      knots = as.numeric(knot_vec),
      degree = as.integer(degree),
      K = as.integer(K),
      Upos = as.matrix(Upos),
      U0 = as.matrix(U0),
      dpos = as.numeric(dpos),
      spline_ci = as.integer(1),
      x_grid = t_grid - 0.5
    )

    obj0 <- MakeADFun(data = Data0, parameters = Params0, random = c("vpos"),
                      DLL = start_name, silent = TRUE)

    opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
    obj0$par <- opt0$par
    rep0 <- sdreport(obj0)
    parList0 <- obj0$env$parList()

    c_hat <- rep0$value[grep("c", names(rep0$value))]
    m_var <- as.numeric(obj0$report()$m)

    #' Fitted starting-value curve, centered on its sample mean.
    #'
    #' Must be visible from .GlobalEnv: otherwise it fails whenever the
    #' nlme() call sits inside a function, as it does here.

    spline_start <- function(t) {
      sweep(as.matrix(splineDesign(knot_vec, t, 4, outer.ok = TRUE))[, 1:Kint, drop = FALSE],
            2, m_var, FUN = "-") %*% c_hat
    }

    assign("spline_start", spline_start, envir = .GlobalEnv)  # removed by the on.exit() above

    ## run nlme model with spline_start
    nlmeobj <- nlme(y ~ b1 + exp(b2) * spline_start(t - exp(b3) / (1 + exp(b3))),
                    fixed = list(b1 ~ 1),
                    random = pdDiag(b1 + b2 + b3 ~ 1),
                    data = dat[, 1:4],
                    groups = ~subject,
                    verbose = FALSE,
                    start = mean(dat$y),
                    control = list(returnObject = TRUE, tolerance = .01))

    # ------------------------------------------------------------------------
    # Fit the full SNMM via TMB
    # ------------------------------------------------------------------------
    Data <- list(
      y = dat$y,
      t = dat$t,
      group = as.integer(dat$subject),
      nGroup = as.integer(n),
      subj_flag = as.integer(subj_flag_vec),
      knots = as.numeric(knot_vec),
      degree = as.integer(degree),
      K = as.integer(K),
      Upos = as.matrix(Upos),
      U0 = as.matrix(U0),
      dpos = as.numeric(dpos),
      spline_ci = as.integer(1),
      t_grid = as.numeric(t_grid)
    )

    # the likelihood takes the centering constants as data
    Data$Bmean <- as.numeric(Bmean)

    Params <- list(
      alpha = nlmeobj$coefficients$fixed["b1"],
      b1 = as.numeric(nlmeobj$coefficients$random$subject[, "b1"]),
      b2 = as.numeric(nlmeobj$coefficients$random$subject[, "b2"]),
      b3 = as.numeric(nlmeobj$coefficients$random$subject[, "b3"]),
      log_sd_b1 = log(as.numeric(VarCorr(nlmeobj)["b1", "StdDev"])),
      log_sd_b2 = log(as.numeric(VarCorr(nlmeobj)["b2", "StdDev"])),
      log_sd_b3 = log(as.numeric(VarCorr(nlmeobj)["b3", "StdDev"])),
      log_sigma = log(nlmeobj$sigma),
      vpos = as.numeric(parList0$vpos),
      gamma0 = as.numeric(parList0$gamma0),
      log_lambda = as.numeric(obj0$par["log_lambda"])
    )

    randoms <- c("b1", "b2", "b3", "vpos")

    start_time <- Sys.time()

    obj <- MakeADFun(data = Data, parameters = Params, random = randoms,
                     DLL = model_name, silent = TRUE)

    # Optimization
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 1e4, iter.max = 1e4))

    obj$par <- opt$par
    rep <- sdreport(obj)

    end_time  <- Sys.time()
    time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # fixed parameters
    true_par <- c(alpha_true, log(sqrt(diag(Sigma_b))), log(sigma))
    bias_vec <- true_par - opt$par[1:5]

    parList <- obj$env$parList()
    # b1_hat, b2_hat, b3_hat
    b1_hat <- if ("b1" %in% names(parList)) parList$b1 else rep(0, n)
    b2_hat <- if ("b2" %in% names(parList)) parList$b2 else rep(0, n)
    b3_hat <- if ("b3" %in% names(parList)) parList$b3 else rep(0, n)

    se <- summary(rep)
    se_b1 <- se[rownames(se) == "b1", "Std. Error"]
    se_b2 <- se[rownames(se) == "b2", "Std. Error"]
    se_b3 <- se[rownames(se) == "b3", "Std. Error"]

    b1_mse <- sqrt(sum((b_mat[, "b1"] - b1_hat)**2) / length(b1_hat))
    b2_mse <- sqrt(sum((b_mat[, "b2"] - b2_hat)**2) / length(b2_hat))
    b3_mse <- sqrt(sum((b_mat[, "b3"] - b3_hat)**2) / length(b3_hat))

    b1_bias <- sum((b_mat[, "b1"] - b1_hat)) / length(b1_hat)
    b2_bias <- sum((b_mat[, "b2"] - b2_hat)) / length(b2_hat)
    b3_bias <- sum((b_mat[, "b3"] - b3_hat)) / length(b3_hat)

    cover_b1 <- mean(b_mat[, "b1"] >= (b1_hat - 1.96 * se_b1) &
                       b_mat[, "b1"] <= (b1_hat + 1.96 * se_b1))
    cover_b2 <- mean(b_mat[, "b2"] >= (b2_hat - 1.96 * se_b2) &
                       b_mat[, "b2"] <= (b2_hat + 1.96 * se_b2))
    cover_b3 <- mean(b_mat[, "b3"] >= (b3_hat - 1.96 * se_b3) &
                       b_mat[, "b3"] <= (b3_hat + 1.96 * se_b3))

    bias_vec <- c(bias_vec, b1_bias, b2_bias, b3_bias)

    # ------------------------------------------------------------------------
    # Population curve coverage
    # ------------------------------------------------------------------------
    mu_true <- alpha_true + 2 * sin(2 * pi * (t_grid - 0.5)) # true population curve

    idx <- grep("h_grid", names(rep$value))
    mu_hat <- rep$value[idx]
    Cov_mu <- rep$cov[idx, idx]

    # symmetric
    Cov_mu <- (Cov_mu + t(Cov_mu)) / 2

    # simulation for simultaneous coverage critical value
    mu_sims <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu)), Sigma = Cov_mu)

    se_point <- sqrt(diag(Cov_mu)) # pointwise SE
    absDev <- abs(sweep(mu_sims, 1, se_point, FUN = "/"))

    # take maximum absolute deviation over the grid for each simulation
    max_dev <- apply(absDev, 1, max)
    crit <- quantile(max_dev, 0.95) # 95% simultaneous critical value

    # save
    pointwise_coverage <- meanCI(mu_true, upr = mu_hat + 1.96 * se_point, lwr = mu_hat - 1.96 * se_point)
    simultaneous_coverage <- inCI(mu_true, upr = mu_hat + crit * se_point, lwr = mu_hat - crit * se_point)
    pointwise_CI_length <- mean(2 * 1.96 * se_point)
    simultaneous_CI_length <- mean(2 * crit * se_point)

    # ------------------------------------------------------------------------
    # Subject curve coverage (for the two subjects flagged in subj_flag_vec)
    # ------------------------------------------------------------------------
    sim_cov_subj <- c(0, 0)
    point_cov_subj <- c(0, 0)

    CI_length_point_subj <- c(0, 0)
    CI_length_sim_subj <- c(0, 0)

    for (subject_index in 1:2) {

      # true individual curve
      b1 <- b_mat[subject_index, "b1"]
      b2 <- b_mat[subject_index, "b2"]
      b3 <- b_mat[subject_index, "b3"]

      amp <- 2 * exp(b2)
      phase <- exp(b3) / (1 + exp(b3)) # logistic(b3)

      mu_true_subj <- alpha_true + b1 + amp * sin(2 * pi * (t_grid - phase))

      # obtain estimated curve and CI
      idx <- grep("mu_sel", names(rep$value))[(ngrid * (subject_index - 1) + 1):(ngrid * subject_index)]
      mu_hat_subj <- rep$value[idx]
      Cov_mu_subj <- rep$cov[idx, idx]

      Cov_mu_subj <- (Cov_mu_subj + t(Cov_mu_subj)) / 2

      mu_sims_subj <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu_subj)), Sigma = Cov_mu_subj)

      se_point_subj <- sqrt(diag(Cov_mu_subj)) # pointwise SE
      absDev_subj <- abs(sweep(mu_sims_subj, 1, se_point_subj, FUN = "/"))

      # Take maximum absolute deviation over the grid for each simulation
      max_dev_subj <- apply(absDev_subj, 1, max)
      crit_subj <- quantile(max_dev_subj, 0.95) # 95% simultaneous critical value

      sim_cov_subj[subject_index] <- inCI(mu_true_subj, upr = mu_hat_subj + as.numeric(crit_subj) * se_point_subj, lwr = mu_hat_subj - as.numeric(crit_subj) * se_point_subj)
      point_cov_subj[subject_index] <- meanCI(mu_true_subj, upr = mu_hat_subj + 1.96 * se_point_subj, lwr = mu_hat_subj - 1.96 * se_point_subj)
      CI_length_sim_subj[subject_index] <- mean(2 * as.numeric(crit_subj) * se_point_subj)
      CI_length_point_subj[subject_index] <- mean(2 * 1.96 * se_point_subj)
    }

    subject_pointwise_coverage <- mean(point_cov_subj)
    subject_simultaneous_coverage <- mean(sim_cov_subj)
    subject_pointwise_CI_length <- mean(CI_length_point_subj)
    subject_simultaneous_CI_length <- mean(CI_length_sim_subj)

    list(
      pointwise_coverage = pointwise_coverage,
      simultaneous_coverage = simultaneous_coverage,
      subject_pointwise_coverage = subject_pointwise_coverage,
      subject_simultaneous_coverage = subject_simultaneous_coverage,
      pointwise_CI_length = pointwise_CI_length,
      simultaneous_CI_length = simultaneous_CI_length,
      subject_pointwise_CI_length = subject_pointwise_CI_length,
      subject_simultaneous_CI_length = subject_simultaneous_CI_length,
      b1_mse = b1_mse,
      b2_mse = b2_mse,
      b3_mse = b3_mse,
      b1_coverage = cover_b1,
      b2_coverage = cover_b2,
      b3_coverage = cover_b3,
      bias_vec = bias_vec,
      b1 = b1_hat,
      b2 = b2_hat,
      b3 = b3_hat,
      time = time_diff
    )

  }, error = function(e) {
    message("Setting ", setting_idx, " iteration ", iter,
            " skipped due to error: ", conditionMessage(e))
    NULL
  })

  out
}

# ==============================================================================
# 9. Task list for this job
# ==============================================================================

tasks <- expand.grid(iter = seq_len(iterations), setting_idx = settings_to_run)
tasks <- tasks[, c("setting_idx", "iter")]
tasks$job <- ((seq_len(nrow(tasks)) - 1L) %% njobs) + 1L

my_tasks <- tasks[tasks$job == itnum, c("setting_idx", "iter")]
rownames(my_tasks) <- NULL
ntask <- nrow(my_tasks)

if (ntask == 0L) stop("Job ", itnum, " of ", njobs, " has no work to do.")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
job_file <- file.path(results_dir,
                      sprintf("%s_job%02d_of%02d.RDS", output_prefix, itnum, njobs))

# ==============================================================================
# 10. Resume
# ==============================================================================

records <- vector("list", ntask)

if (file.exists(job_file)) {
  prev <- readRDS(job_file)
  same_run <- identical(prev$method, "snmmTMB") &&
    identical(as.integer(prev$njobs), njobs) &&
    identical(as.integer(prev$job), itnum) &&
    identical(as.integer(prev$iterations), as.integer(iterations)) &&
    identical(prev$seed_base, seed_base) &&
    identical(as.numeric(prev$err_df), as.numeric(err_df)) &&
    identical(prev$tasks, my_tasks)
  if (!same_run) {
    stop("Existing ", job_file, " was written with a different configuration ",
         "(njobs / iterations / seed_base / err_df / settings). Move it aside or fix the configuration.")
  }
  records <- prev$records
  message("Resuming ", job_file, ": ",
          sum(!vapply(records, is.null, logical(1))), "/", ntask, " already done.")
}

#' Write this job's file to a temporary file
save_job <- function() {
  tmp <- paste0(job_file, ".tmp")
  saveRDS(list(
    method = "snmmTMB",
    job = itnum,
    njobs = njobs,
    iterations = iterations,
    settings_to_run = settings_to_run,
    seed_base = seed_base,
    err_df = err_df,
    param_grid = param_grid,
    tasks = my_tasks,
    records = records,
    finished_at = Sys.time(),
    sessionInfo = utils::sessionInfo()
  ), tmp)
  invisible(file.rename(tmp, job_file))
}

# ==============================================================================
# 11. Run
# ==============================================================================

message("snmmTMB job ", itnum, "/", njobs, ": ", ntask, " replications ",
        "(", iterations, " iterations x ", length(settings_to_run), " settings, ",
        "seed_base = ", seed_base, ", err_df = ", err_df, ")")

job_start <- Sys.time()

for (k in seq_len(ntask)) {

  if (!is.null(records[[k]])) next

  s <- my_tasks$setting_idx[k]
  it <- my_tasks$iter[k]

  res <- run_replication(s, it)

  records[[k]] <- list(setting_idx = s, iter = it, ok = !is.null(res), metrics = res)

  message(sprintf("[%s] job %d/%d  task %d/%d  setting %2d  iter %3d  %s  (%.1f min elapsed)",
                  format(Sys.time(), "%H:%M:%S"), itnum, njobs, k, ntask, s, it,
                  if (is.null(res)) "FAILED" else sprintf("ok, %.1fs", res$time),
                  as.numeric(difftime(Sys.time(), job_start, units = "mins"))))

  if (k %% save_every == 0L || k == ntask) save_job()
}

save_job()

n_ok <- sum(vapply(records, function(r) isTRUE(r$ok), logical(1)))
message(sprintf("DONE snmmTMB job %d/%d: %d/%d replications succeeded in %.1f min -> %s",
                itnum, njobs, n_ok, ntask,
                as.numeric(difftime(Sys.time(), job_start, units = "mins")), job_file))
