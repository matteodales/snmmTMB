# ------------------------------------------------------------------------------
# bellcurve_server_simulation_assist.R
#
# assist script for the bell-curve simulation of manuscript Section 3 (Figure 4)
#
# This script runs the full study and is meant to run on a computing cluster
# One iteration of the same model can be fitted and plotted by running bellcurve_example.R
#
# Usage, one process per job:
#
#   Rscript bellcurve_server_simulation_assist.R <job> <njobs>
#
# with <job> in 1..<njobs>; bellcurve_server_run_assist.sh launches all of the jobs at once.
# Each job writes only its own results/<output_prefix>_job<j>_of<njobs>.RDS
# ------------------------------------------------------------------------------

# ==============================================================================
# 0. Command-line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript bellcurve_server_simulation_assist.R <job> [<njobs>]")
}

itnum <- as.integer(args[1])
njobs <- if (length(args) >= 2) as.integer(args[2]) else 20L

if (is.na(itnum) || is.na(njobs) || njobs < 1L || itnum < 1L || itnum > njobs) {
  stop("Bad arguments: need 1 <= job <= njobs, got job = ", args[1],
       ", njobs = ", if (length(args) >= 2) args[2] else njobs)
}

# ==============================================================================
# 1. Run configuration
# ===============================================================================

iterations <- 500            # replications per setting
settings_to_run <- 1:8
seed_base <- 20260812         
Nsims <- 10000                 # draws for the simultaneous-band critical value

results_dir <- "results"
output_prefix <- "bellcurve_server_assist"
save_every <- 1L               # save the job file every N replications

# ==============================================================================
# 2. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm() for multivariate-normal random-effect draws
  library(dplyr)    # bind_rows()
  library(tibble)   # tibble() in place of data.frame()
  library(nlme)     # VarCorr(), pdSymm(), pdMatrix() on the fit's nlmeObj
  library(assist)   # snm(), cubic(), intervals()
})

# ==============================================================================
# 3. Coverage helpers
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

#' Reconstructed simultaneous band for an snm() population curve
#'
# We obtain the critical value by simulation as the other two methods do.
#
# All required quantities are taken from snm_fit$forCI object, which is used for assist's own pstd
#'
#' @param fit an snm() fit
#' @param u_new evaluation grid, on f's own scale (t_grid / 10 + 1/2 here)
#' @param nsim draws for the critical value
#' @param level coverage level
#' @param eig_tol keep kernel eigen-directions above eig_tol * the largest one
#' @return list(f_rec, se_rec, crit, edf, edf_target, nbasis, lam)
assist_recon_band <- function(fit, u_new, nsim = Nsims, level = 0.95,
                              eig_tol = 1e-12) {

  fc <- fit$forCI
  ro <- fc$rkpk.obj

  u_obs <- fc$data[[1]]$u        # u_i = t_scaled_ij - b2_hat_i, as assist used them
  d1    <- fc$delta1[, 1]        # d mu / d f: all ones for this model

  Vinv <- crossprod(fc$weight)
  lam  <- 10^ro$nlaht
  sig2 <- ro$varht

  K_obs <- fc$q[[1]] / tcrossprod(d1)   # undo the (unit) scaling
  K_new <- cubic(u_new, u_obs)

  eg   <- eigen((K_obs + t(K_obs)) / 2, symmetric = TRUE)
  keep <- eg$values > eig_tol * eg$values[1]
  Uk   <- eg$vectors[, keep, drop = FALSE]
  Lk   <- eg$values[keep]
  nk   <- length(Lk)
  if (nk < 2) stop("kernel eigenbasis collapsed to ", nk, " directions")

  Xt     <- cbind(1, d1 * cbind(u_obs))       # model space is ~u - 1: one column
  Xt_new <- cbind(1, cbind(u_new))
  Xc     <- d1 * (Uk %*% diag(sqrt(Lk), nrow = nk))
  Xc_new <- K_new %*% (Uk %*% diag(1 / sqrt(Lk), nrow = nk))

  X     <- cbind(Xt, Xc)
  X_new <- cbind(Xt_new, Xc_new)
  M     <- ncol(Xt)
  P     <- diag(c(rep(0, M), rep(1, nk)))     # zero on the model space

  XtVX     <- t(X) %*% Vinv %*% X
  Ainv     <- solve(XtVX + lam * P)
  beta_hat <- as.numeric(Ainv %*% t(X) %*% Vinv %*% fc$y)
  Cov_beta <- sig2 * Ainv

  Cov_curve <- X_new %*% Cov_beta %*% t(X_new)
  Cov_curve <- (Cov_curve + t(Cov_curve)) / 2
  se_rec    <- sqrt(pmax(diag(Cov_curve), 0))

  sims <- mvrnorm(nsim, mu = rep(0, length(se_rec)), Sigma = Cov_curve)
  r    <- apply(abs(sweep(sims, 2, se_rec, FUN = "/")), 1, max)

  list(f_rec      = as.numeric(X_new %*% beta_hat),
       se_rec     = se_rec,
       crit       = as.numeric(quantile(r, level)),
       edf        = sum(diag(Ainv %*% XtVX)),
       edf_target = ro$df,
       nbasis     = nk,
       lam        = lam)
}

# ==============================================================================
# 4. Simulation grid
#
# 8 settings = 2 (n) x 2 (m) x 2 (sigma), with D = [2 1; 1 2], i.e. a
# correlation rho = 0.5 between the two random effects
# ==============================================================================

n_vals <- c(10, 50)
m_vals <- c(10, 20)
sigma_vals <- c(0.2, 0.4)
D_diag <- c(2, 2)              # diagonal of D
rho_vals <- c(0.5)             # correlation between b1 and b2

param_grid <- expand.grid(
  n = n_vals,
  m = m_vals,
  sigma = sigma_vals,
  rho = rho_vals
)

t_grid <- seq(-2, 2, length.out = 50)  # grid the curves are reported on
alpha_true <- 1

# ==============================================================================
# 5. One replication
# ==============================================================================

#' @param setting_idx row of param_grid
#' @param iter replication number within the setting
#' @return integer seed
seed_for <- function(setting_idx, iter) {
  seed_base + 977L * as.integer(setting_idx) + as.integer(iter)
}

#' Simulate one bell-curve dataset
#' @param n,m subjects and observations per subject
#' @param sigma residual sd; D_diag length-2 diagonal of D; rho correlation
#' @param alpha population intercept
#' @return list(dat = tibble, b_mat = n x 2 matrix of true random effects,
#'   Sigma_b = the random-effect covariance used)
simulate_dataset <- function(n, m, sigma, D_diag, rho, alpha) {
  cov_b12 <- rho * sqrt(D_diag[1] * D_diag[2])
  Sigma_b <- sigma^2 * matrix(c(D_diag[1], cov_b12,
                                cov_b12,   D_diag[2]), nrow = 2, byrow = TRUE)
  b_mat <- mvrnorm(n = n, mu = rep(0, 2), Sigma = Sigma_b)
  colnames(b_mat) <- c("b1", "b2")

  dat_list <- vector("list", n)
  for (i in 1:n) {
    b1 <- b_mat[i, "b1"]
    b2 <- b_mat[i, "b2"]

    t <- seq(-2, 2, length.out = m)

    mu_ij <- alpha + b1 + exp(-0.5 * (t - b2)^2)
    y_ij  <- mu_ij + rnorm(m, mean = 0, sd = sigma)

    dat_list[[i]] <- tibble(
      subject = rep(i, m),
      j = 1:m,
      t = t,
      y = y_ij,
      mu = mu_ij,
      b1 = rep(b1, m),
      b2 = rep(b2, m)
    )
  }

  list(dat = bind_rows(dat_list), b_mat = b_mat, Sigma_b = Sigma_b)
}

#' Run one replication end to end
#' @param setting_idx row of param_grid
#' @param iter replication number within the setting
#' @return named list of metrics, or NULL on failure
run_replication <- function(setting_idx, iter) {


  on.exit(if (exists("start_val", envir = .GlobalEnv, inherits = FALSE)) {
    rm("start_val", envir = .GlobalEnv)
  }, add = TRUE)

  # snm() assigns its spline closure to `f`
  on.exit({
    if (exists("f", envir = .GlobalEnv, inherits = FALSE) &&
        is.function(get("f", envir = .GlobalEnv))) {
      rm("f", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed_for(setting_idx, iter),
           kind = "Mersenne-Twister", normal.kind = "Inversion")

  n <- param_grid$n[setting_idx]
  m <- param_grid$m[setting_idx]
  sigma <- param_grid$sigma[setting_idx]
  rho <- param_grid$rho[setting_idx]

  sim <- simulate_dataset(n, m, sigma, D_diag, rho, alpha_true)
  dat <- sim$dat
  b_mat <- sim$b_mat
  Sigma_b <- sim$Sigma_b


  mu_true <- alpha_true + exp(-0.5 * t_grid^2)

  # drop b1/b2/mu/j and keep only what snm() needs; cubic() expects its
  # argument roughly in [0,1], so t in [-2,2] is mapped to t / 10 + 1/2
  dat_fit <- as.data.frame(dat[, c("subject", "t", "y")])
  dat_fit$subject <- as.factor(dat_fit$subject)
  dat_fit$t_scaled <- dat_fit$t / 10 + 1 / 2

  start_val <- mean(dat_fit$y)
  assign("start_val", start_val, envir = .GlobalEnv)  # removed by the on.exit() above

  out <- tryCatch({

    start_time <- Sys.time()

    # general (correlated) random-effect covariance, initialised at Sigma_b
    snm_fit <- snm(
      y ~ b1 + f(t_scaled - b2),
      func   = f(u) ~ list(~u - 1, cubic(u)),
      fixed  = list(b1 ~ 1),
      random = pdSymm(b1 + b2 ~ 1, value = Sigma_b),
      data   = dat_fit,
      groups = ~subject,
      start  = start_val,
      verbose = FALSE
    )

    end_time  <- Sys.time()
    time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))

    true_par <- c(alpha_true, log(sqrt(diag(Sigma_b))), log(sigma))
    est_par <- c(snm_fit$coefficients$fixed,
                 log(as.numeric(as.matrix(VarCorr(snm_fit$nlmeObj))[, "StdDev"])))
    bias_vec <- true_par - est_par

    # b1_hat, b2_hat
    b1_hat <- snm_fit$coefficients$random$subject[, "b1"]
    b2_hat <- snm_fit$coefficients$random$subject[, "b2"]

    b1_mse <- sqrt(sum((b_mat[, "b1"] - b1_hat)**2) / length(b1_hat))
    b2_mse <- sqrt(sum((b_mat[, "b2"] - b2_hat)**2) / length(b2_hat))

    b1_bias <- sum((b_mat[, "b1"] - b1_hat)) / length(b1_hat)
    b2_bias <- sum((b_mat[, "b2"] - b2_hat)) / length(b2_hat)

    bias_vec <- c(bias_vec, b1_bias, b2_bias)

    # estimated correlation, read off the pdSymm structure (no standard error)
    rho_hat <- tryCatch({
      vc <- nlme::pdMatrix(snm_fit$nlmeObj$modelStruct$reStruct[[1]]) *
            snm_fit$nlmeObj$sigma^2
      r <- vc[1, 2] / sqrt(vc[1, 1] * vc[2, 2])
      if (is.finite(r)) r else NA_real_
    }, error = function(e) {
      tryCatch(as.numeric(as.matrix(VarCorr(snm_fit$nlmeObj))[2, "Corr"]),
               error = function(e2) NA_real_)
    })
    rho_bias <- rho - rho_hat

    # Population curve coverage
    bci <- intervals(snm_fit, newdata = data.frame(u = t_grid / 10 + 1 / 2))
    pred_vec <- as.numeric(as.matrix(bci$fit)[, 1]) + snm_fit$coefficients$fixed
    pstd_vec <- as.numeric(as.matrix(bci$pstd))

    band <- assist_recon_band(snm_fit, u_new = t_grid / 10 + 1 / 2,
                              nsim = Nsims, level = 0.95)
    crit <- band$crit
    se_sim <- band$se_rec

    ci_df <- tibble(
      t = t_grid,
      mu_est = pred_vec,
      lower = pred_vec - qnorm(0.975) * pstd_vec,
      upper = pred_vec + qnorm(0.975) * pstd_vec,
      lower_sim = pred_vec - crit * se_sim,
      upper_sim = pred_vec + crit * se_sim
    )

    list(
      pointwise_coverage = meanCI(mu_true, upr = ci_df$upper, lwr = ci_df$lower),
      simultaneous_coverage = inCI(mu_true, upr = ci_df$upper_sim,
                                   lwr = ci_df$lower_sim),
      pointwise_CI_length = mean(ci_df$upper - ci_df$lower),
      simultaneous_CI_length = mean(ci_df$upper_sim - ci_df$lower_sim),
      simultaneous_crit = crit,
      b1_mse = b1_mse,
      b2_mse = b2_mse,
      rho_true = rho,
      rho_hat = rho_hat,
      rho_bias = rho_bias,
      bias_vec = bias_vec,
      b1 = b1_hat,
      b2 = b2_hat,
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
# 6. Task list for this job
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
# 7. Resume
# ==============================================================================

records <- vector("list", ntask)

if (file.exists(job_file)) {
  prev <- readRDS(job_file)
  same_run <- identical(prev$method, "assist") &&
    identical(as.integer(prev$njobs), njobs) &&
    identical(as.integer(prev$job), itnum) &&
    identical(as.integer(prev$iterations), as.integer(iterations)) &&
    identical(prev$seed_base, seed_base) &&
    identical(prev$tasks, my_tasks)
  if (!same_run) {
    stop("Existing ", job_file, " was written with a different configuration ",
         "(njobs / iterations / seed_base / settings). Move it aside or fix the configuration.")
  }
  records <- prev$records
  message("Resuming ", job_file, ": ",
          sum(!vapply(records, is.null, logical(1))), "/", ntask, " already done.")
}

#' Write this job's file to a temporary file
save_job <- function() {
  tmp <- paste0(job_file, ".tmp")
  saveRDS(list(
    method = "assist",
    job = itnum,
    njobs = njobs,
    iterations = iterations,
    settings_to_run = settings_to_run,
    seed_base = seed_base,
    param_grid = param_grid,
    tasks = my_tasks,
    records = records,
    finished_at = Sys.time(),
    sessionInfo = utils::sessionInfo()
  ), tmp)
  invisible(file.rename(tmp, job_file))
}

# ==============================================================================
# 8. Run
# ==============================================================================

message("assist job ", itnum, "/", njobs, ": ", ntask, " replications ",
        "(", iterations, " iterations x ", length(settings_to_run), " settings, ",
        "seed_base = ", seed_base, ")")

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
message(sprintf("DONE assist job %d/%d: %d/%d replications succeeded in %.1f min -> %s",
                itnum, njobs, n_ok, ntask,
                as.numeric(difftime(Sys.time(), job_start, units = "mins")), job_file))
