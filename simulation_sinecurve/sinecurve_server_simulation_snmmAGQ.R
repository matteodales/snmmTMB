# ------------------------------------------------------------------------------
# sinecurve_server_simulation_snmmAGQ.R
#
# snmmAGQ script for the sine-curve simulation of manuscript Section 3 (Figures 1 and 2)
#
# This script runs the full study and is meant to run on a computing cluster
# One iteration of the same model can be fitted and plotted by running sinecurve_example.R
#
# Usage, one process per job:
#
#   Rscript sinecurve_server_simulation_snmmAGQ.R <job> <njobs>
#
# with <job> in 1..<njobs>; sinecurve_server_run_snmmAGQ.sh launches all of the jobs at once.
# Each job writes only its own results/<output_prefix>_job<j>_of<njobs>.RDS.
#
# The estimator is sourced from ../src/snmmAGQ_functions_sinecurve.R. The number
# of interior knots is selected by AIC over nk in 0..6 for every dataset
# (fit_select_nk()).
# ------------------------------------------------------------------------------

# ==============================================================================
# 0. Command-line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript sinecurve_server_simulation_snmmAGQ.R <job> [<njobs>]")
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
seed_base <- 20260804
Nsims <- 10000                 # draws for the simultaneous-band critical value

src_dir <- "../src"
results_dir <- "results"
output_prefix <- "sinecurve_server_snmmAGQ"
save_every <- 1L               # flush the job file every N replications

# ==============================================================================
# 2. Packages and estimator
# ==============================================================================

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm() for multivariate-normal random-effect draws
  library(dplyr)    # bind_rows()
  library(tibble)   # tibble() in place of data.frame()
})

# fit_select_nk(), neglogLik(), subj_curve_fun(), make_pd(); loads splines,
# nlme, statmod and numDeriv
source(file.path(src_dir, "snmmAGQ_functions_sinecurve.R"))

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

# ==============================================================================
# 4. Simulation grid
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

t_grid <- seq(0, 1, length.out = 50)  # grid the curves are reported on
alpha_true <- 5

# ==============================================================================
# 5. One replication
# ==============================================================================

#' @param setting_idx row of param_grid
#' @param iter replication number within the setting
#' @return integer seed
seed_for <- function(setting_idx, iter) {
  seed_base + 977L * as.integer(setting_idx) + as.integer(iter)
}

#' Simulate one sine-curve dataset
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
    y_ij  <- mu_ij + rnorm(m, mean = 0, sd = sigma)

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

#' Run one replication end to end
#' @param setting_idx row of param_grid
#' @param iter replication number within the setting
#' @return named list of metrics, or NULL on failure
run_replication <- function(setting_idx, iter) {

  set.seed(seed_for(setting_idx, iter),
           kind = "Mersenne-Twister", normal.kind = "Inversion")

  n <- param_grid$n[setting_idx]
  m <- param_grid$m[setting_idx]
  sigma <- param_grid$sigma[setting_idx]
  D_diag <- D_diag_vals[[param_grid$D_idx[setting_idx]]]

  sim <- simulate_dataset(n, m, sigma, D_diag, alpha_true)
  dat <- sim$dat
  b_mat <- sim$b_mat
  Sigma_b <- sim$Sigma_b

  # per-subject responses and observation times, as bspline.snm() takes them
  subj_idx <- split(seq_len(nrow(dat)), dat$subject)
  y_list <- lapply(subj_idx, function(idx) dat$y[idx])
  t_list <- lapply(subj_idx, function(idx) dat$t[idx])

  out <- tryCatch({

    start_time <- Sys.time()

    # Fit at every nk in nk_grid
    fit <- fit_select_nk(dat, y_list, t_list)

    a <- fit$a; alpha <- fit$alpha; knotseq <- fit$knotseq; col_means <- fit$col_means
    theta <- fit$theta; modes <- fit$modes; hessians <- fit$hessians; cache <- fit$cache

    # ------------------------------------------------------------------------
    # Standard errors (Equation 2.16 of Elmi et al. 2011): G is the numerical
    # Hessian of the marginal log-likelihood in the fixed effects (alpha, a)
    # ------------------------------------------------------------------------
    G <- numDeriv::hessian(function(par){
      alpha_par <- par[1]; a_par <- par[-1]
      neglogLik(cache, a_par, alpha_par, theta)
    }, c(alpha, a))
    Ginv <- solve(G)

    end_time  <- Sys.time()
    time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # variance parameters
    true_par <- c(sqrt(diag(Sigma_b)), sigma)
    bias_vec <- true_par - theta

    # b1_hat, b2_hat, b3_hat: the posterior modes
    b1_hat <- modes[,1]; b2_hat <- modes[,2]; b3_hat <- modes[,3]

    b1_mse <- sqrt(sum((b_mat[,'b1']-b1_hat)^2)/n)
    b2_mse <- sqrt(sum((b_mat[,'b2']-b2_hat)^2)/n)
    b3_mse <- sqrt(sum((b_mat[,'b3']-b3_hat)^2)/n)

    b1_bias <- mean(b_mat[,'b1']-b1_hat)
    b2_bias <- mean(b_mat[,'b2']-b2_hat)
    b3_bias <- mean(b_mat[,'b3']-b3_hat)
    alpha_bias <- alpha_true - alpha
    bias_vec <- c(bias_vec, b1_bias, b2_bias, b3_bias, alpha_bias)

    se_b_all <- t(sapply(hessians, function(h) sqrt(diag(solve(h)))))
    se_b1 <- se_b_all[,1]; se_b2 <- se_b_all[,2]; se_b3 <- se_b_all[,3]

    cover_b1 <- mean(b_mat[,"b1"] >= (b1_hat-1.96*se_b1) & b_mat[,"b1"] <= (b1_hat+1.96*se_b1))
    cover_b2 <- mean(b_mat[,"b2"] >= (b2_hat-1.96*se_b2) & b_mat[,"b2"] <= (b2_hat+1.96*se_b2))
    cover_b3 <- mean(b_mat[,"b3"] >= (b3_hat-1.96*se_b3) & b_mat[,"b3"] <= (b3_hat+1.96*se_b3))

    # ------------------------------------------------------------------------
    # Population curve coverage
    # ------------------------------------------------------------------------
    mu_true <- alpha_true + 2 * sin(2 * pi * (t_grid - 0.5)) # true population curve

    X_pop <- cbind(1, sweep(splines::splineDesign(knotseq, t_grid - 0.5, 4, outer.ok = TRUE), 2, col_means, FUN = "-"))
    mu_hat <- as.numeric(X_pop %*% c(alpha, a))

    # the population curve does not depend on b_i, so its covariance is X Ginv X'
    Cov_mu <- make_pd(X_pop %*% Ginv %*% t(X_pop))
    se_point <- sqrt(diag(Cov_mu)) # pointwise SE

    # simulation for simultaneous coverage critical value
    mu_sims <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu)), Sigma = Cov_mu)
    absDev <- abs(sweep(mu_sims, 2, se_point, FUN = "/"))
    crit <- quantile(apply(absDev, 1, max), 0.95) # 95% simultaneous critical value

    pointwise_coverage <- meanCI(mu_true, upr = mu_hat + 1.96*se_point, lwr = mu_hat - 1.96*se_point)
    simultaneous_coverage <- inCI(mu_true, upr = mu_hat + crit*se_point, lwr = mu_hat - crit*se_point)
    pointwise_CI_length <- mean(2*1.96*se_point)
    simultaneous_CI_length <- mean(2*crit*se_point)

    # ------------------------------------------------------------------------
    # Subject curve coverage (subjects 1 and 2)
    # ------------------------------------------------------------------------
    sim_cov_subj <- c(0,0); point_cov_subj <- c(0,0)
    CI_length_point_subj <- c(0,0); CI_length_sim_subj <- c(0,0)

    for(subject_index in 1:2){

      # true individual curve
      b1s <- b_mat[subject_index,"b1"]; b2s <- b_mat[subject_index,"b2"]; b3s <- b_mat[subject_index,"b3"]
      amp_s <- 2 * exp(b2s); phase_s <- exp(b3s) / (1 + exp(b3s))
      mu_true_subj <- alpha_true + b1s + amp_s * sin(2 * pi * (t_grid - phase_s))

      # estimated curve and CI
      b_hat_i <- modes[subject_index,]
      X_i <- cbind(1, exp(b_hat_i[2]) * sweep(splines::splineDesign(knotseq, t_grid - expit(b_hat_i[3]), 4, outer.ok=TRUE), 2, col_means, FUN = "-"))
      Z_i <- numDeriv::jacobian(function(b) subj_curve_fun(b, alpha, a, knotseq, col_means), b_hat_i)

      Var_mu_i <- make_pd(X_i %*% Ginv %*% t(X_i) + Z_i %*% solve(hessians[[subject_index]]) %*% t(Z_i))
      se_i <- sqrt(diag(Var_mu_i))

      mu_sims_i <- mvrnorm(Nsims, mu = rep(0, ncol(Var_mu_i)), Sigma = Var_mu_i)
      absDev_i <- abs(sweep(mu_sims_i, 2, se_i, FUN = "/"))
      crit_i <- quantile(apply(absDev_i, 1, max), 0.95)

      mu_hat_subj <- subj_curve_fun(b_hat_i, alpha, a, knotseq, col_means)

      sim_cov_subj[subject_index] <- inCI(mu_true_subj, upr = mu_hat_subj + as.numeric(crit_i)*se_i, lwr = mu_hat_subj - as.numeric(crit_i)*se_i)
      point_cov_subj[subject_index] <- meanCI(mu_true_subj, upr = mu_hat_subj + 1.96*se_i, lwr = mu_hat_subj - 1.96*se_i)
      CI_length_sim_subj[subject_index] <- mean(2*as.numeric(crit_i)*se_i)
      CI_length_point_subj[subject_index] <- mean(2*1.96*se_i)
    }

    list(
      pointwise_coverage = pointwise_coverage,
      simultaneous_coverage = simultaneous_coverage,
      subject_pointwise_coverage = mean(point_cov_subj),
      subject_simultaneous_coverage = mean(sim_cov_subj),
      pointwise_CI_length = pointwise_CI_length,
      simultaneous_CI_length = simultaneous_CI_length,
      subject_pointwise_CI_length = mean(CI_length_point_subj),
      subject_simultaneous_CI_length = mean(CI_length_sim_subj),
      b1_mse = b1_mse, b2_mse = b2_mse, b3_mse = b3_mse,
      b1_coverage = cover_b1, b2_coverage = cover_b2, b3_coverage = cover_b3,
      bias_vec = bias_vec,
      b1 = b1_hat, b2 = b2_hat, b3 = b3_hat,
      nk_selected = fit$nk_selected,   # AIC-selected number of interior knots
      nk_aic = fit$nk_aic,             # AIC at every nk in nk_grid (NA if that fit failed)
      nk_ok = fit$nk_ok,               # number of candidates that converged
      time = time_diff,
      boundary_any = any(fit$boundary) # a variance component collapsed to zero (degenerate fit)
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
  same_run <- identical(prev$method, "snmmAGQ") &&
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
    method = "snmmAGQ",
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

message("snmmAGQ job ", itnum, "/", njobs, ": ", ntask, " replications ",
        "(", iterations, " iterations x ", length(settings_to_run), " settings, ",
        "seed_base = ", seed_base, ")")

job_start <- Sys.time()

for (k in seq_len(ntask)) {

  if (!is.null(records[[k]])) next

  s <- my_tasks$setting_idx[k]
  it <- my_tasks$iter[k]

  # Crash guard, gets overwritten if the run finishes
  records[[k]] <- list(setting_idx = s, iter = it, ok = FALSE, metrics = NULL,
                       crashed = TRUE)
  save_job()

  res <- run_replication(s, it)

  records[[k]] <- list(setting_idx = s, iter = it, ok = !is.null(res), metrics = res)

  message(sprintf("[%s] job %d/%d  task %d/%d  setting %2d  iter %3d  %s  (%.1f min elapsed)",
                  format(Sys.time(), "%H:%M:%S"), itnum, njobs, k, ntask, s, it,
                  if (is.null(res)) "FAILED" else
                    sprintf("ok, %.1fs%s", res$time,
                            if (isTRUE(res$boundary_any)) " [BOUNDARY]" else ""),
                  as.numeric(difftime(Sys.time(), job_start, units = "mins"))))

  if (k %% save_every == 0L || k == ntask) save_job()
}

save_job()

n_ok <- sum(vapply(records, function(r) isTRUE(r$ok), logical(1)))
n_bnd <- sum(vapply(records, function(r) isTRUE(r$metrics$boundary_any), logical(1)))
message(sprintf("DONE snmmAGQ job %d/%d: %d/%d replications succeeded (%d degenerate) in %.1f min -> %s",
                itnum, njobs, n_ok, ntask, n_bnd,
                as.numeric(difftime(Sys.time(), job_start, units = "mins")), job_file))
