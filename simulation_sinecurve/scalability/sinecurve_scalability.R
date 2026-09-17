# ------------------------------------------------------------------------------
# sinecurve_scalability.R
#
# Scalability study for the sine-curve simulation of manuscript Section 3 (Figure 3)
#
# Times snmmTMB, assist and snmmAGQ on the same simulated datasets, varying one quantity at a time:
#
#   subjects       n in 10, 20, 40, 60, 80, 100
#   observations   m in 5, 10, 20, 40, 80
#   basis          K in 8, 10, 15, 20, 30, 45     (snmmTMB and snmmAGQ only)
#
#
# This script runs the full study and is meant to run on a computing cluster
#
# Usage, one process per job:
#
#   Rscript sinecurve_server_compile.R                  # once, first, from ../
#   Rscript sinecurve_scalability.R <job> <njobs>
#
# with <job> in 1..<njobs>; sinecurve_server_run_scalability.sh launches all of the jobs at once.
# Each job writes only its own results/<output_prefix>_job<j>_of<njobs>.RDS
#
# To measure timing equally, so every job is launched on one core
# with a single-threaded BLAS, and deals tasks so that machine load is the 
# same on average for every cell and method.
# ------------------------------------------------------------------------------

# ==============================================================================
# 0. Command-line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript sinecurve_scalability.R <job> [<njobs>]")
}

itnum <- as.integer(args[1])
njobs <- if (length(args) >= 2) as.integer(args[2]) else 50L

if (is.na(itnum) || is.na(njobs) || njobs < 1L || itnum < 1L || itnum > njobs) {
  stop("Bad arguments: need 1 <= job <= njobs, got job = ", args[1],
       ", njobs = ", if (length(args) >= 2) args[2] else njobs)
}

# ==============================================================================
# 1. Run configuration
# ==============================================================================

iterations <- 50               # replications per cell
seed_base <- 20260825

src_dir <- "../../src"
results_dir <- "results"
output_prefix <- "sinecurve_scalability"
save_every <- 1L

# data-generating scenario
sigma_true <- 1
D_diag <- c(1, 0.25, 0.16)
alpha_true <- 5

# three sweeps
n_grid_subjects <- c(10, 20, 40, 60, 80, 100)
m_grid_observations <- c(5, 10, 20, 40, 80)
K_grid_basis <- c(8, 10, 15, 20, 30, 45)

n_base <- 20L
m_base <- 10L
K_base <- 15L                  # snmmTMB basis functions, as in the main simulation

# ==============================================================================
# 2. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm() for multivariate-normal random-effect draws
  library(dplyr)    # bind_rows()
  library(tibble)   # tibble() in place of data.frame()
  library(TMB)      # automatic differentiation / Laplace approximation
  library(nlme)     # nlme() for the snmmTMB starting values
  library(splines)  # splineDesign() for the starting-value GAM and Bmean
  library(assist)   # snm(), lspline(), alogit(), intervals()
})

# fit_select_nk(), neglogLik() and the knot grid nk_grid
source(file.path(src_dir, "snmmAGQ_functions_sinecurve.R"))

# single-threaded BLAS
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

# ==============================================================================
# 3. Knot vector, penalty matrix and centering constants
# ==============================================================================

degree <- 3

knot_lowlim <- -1
knot_uplim <- 1

t_grid <- seq(0, 1, length.out = 50)  # grid the population curve is reported on

# fixed-grid centering constants
center_grid_range <- c(-1, 1)
center_grid_size <- 1000

#' B-spline basis, penalty eigendecomposition and centering constants
#'
#' @param K number of basis functions
#' @return list(K, Kint = K - 1, knot_vec, Upos, dpos, U0, Bmean)
make_basis <- function(K) {

  K <- as.integer(K)
  Kint <- K - 1L
  n_internal <- K - (degree + 1L)

  internal_knots <- seq(knot_lowlim, knot_uplim, length.out = n_internal + 2)
  dx_left  <- diff(internal_knots)[1]
  dx_right <- diff(internal_knots)[length(internal_knots) - 1]

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

  Upos <- if (length(pos_idx) > 0) ev$vectors[, pos_idx, drop = FALSE] else matrix(0, nrow = Kint, ncol = 0)
  dpos <- if (length(pos_idx) > 0) ev$values[pos_idx] else numeric(0)
  U0   <- if (length(zero_idx) > 0) ev$vectors[, zero_idx, drop = FALSE] else matrix(0, nrow = Kint, ncol = 0)

  # column means of the B-spline basis over the fixed grid
  edges <- seq(center_grid_range[1], center_grid_range[2], length.out = center_grid_size + 1)
  u_grid <- (head(edges, -1) + tail(edges, -1)) / 2
  Bmean <- colMeans(splineDesign(knot_vec, u_grid, ord = degree + 1,
                                 outer.ok = TRUE)[, 1:Kint, drop = FALSE])

  list(K = K, Kint = Kint, knot_vec = knot_vec,
       Upos = Upos, dpos = dpos, U0 = U0, Bmean = Bmean)
}

# one basis per distinct K
K_all <- sort(unique(c(K_base, K_grid_basis)))
basis_by_K <- lapply(K_all, make_basis)
names(basis_by_K) <- as.character(K_all)

# ==============================================================================
# 4. Simulation grid
# ==============================================================================

cell_grid <- bind_rows(
  tibble(sweep = "subjects", n = as.integer(n_grid_subjects), m = m_base, K = K_base),
  tibble(sweep = "observations", n = n_base, m = as.integer(m_grid_observations), K = K_base),
  tibble(sweep = "basis", n = n_base, m = m_base, K = as.integer(K_grid_basis))
) %>%
  mutate(cell = row_number(), .before = 1) %>%
  # snmmAGQ interior knots: NA = selected by AIC over nk_grid, as in the main simulation
  mutate(nk = if_else(sweep == "basis", K - 4L, NA_integer_))

# ==============================================================================
# 5. Load the compiled likelihoods
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
# 6. One replication
# ==============================================================================
#'
#' @param n,m subjects and observations per subject
#' @param iter replication number within the cell
#' @return integer seed; depends on (n, m, iter)
seed_for <- function(n, m, iter) {
  seed_base + 977L * as.integer(n) + 7919L * as.integer(m) + as.integer(iter)
}

#' Simulate one sine-curve dataset
#'
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

#' Time snmmTMB on one dataset
#'
#' @param dat simulated dataset
#' @param n number of subjects
#' @param basis one element of basis_by_K
#' @return list(ok, time_total); time_total is NA on failure
fit_snmmTMB <- function(dat, n, basis) {

  on.exit(if (exists("spline_start", envir = .GlobalEnv, inherits = FALSE)) {
    rm("spline_start", envir = .GlobalEnv)
  }, add = TRUE)

  out <- tryCatch({

    start_time <- Sys.time()

    Params0 <- list(
      beta1 = as.numeric(0), # intercept
      log_sigma = log(1),
      vpos = rep(0, length(basis$dpos)),
      gamma0 = rep(0, ncol(basis$U0)),
      log_lambda = log(1)
    )

    Data0 <- list(
      y = dat$y,
      x = dat$t - 0.5,
      knots = as.numeric(basis$knot_vec),
      degree = as.integer(degree),
      K = as.integer(basis$K),
      Upos = as.matrix(basis$Upos),
      U0 = as.matrix(basis$U0),
      dpos = as.numeric(basis$dpos),
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
      sweep(as.matrix(splineDesign(basis$knot_vec, t, 4, outer.ok = TRUE))[, 1:basis$Kint, drop = FALSE],
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
      subj_flag = as.integer(c(1, 1, rep(0, n - 2))),
      knots = as.numeric(basis$knot_vec),
      degree = as.integer(degree),
      K = as.integer(basis$K),
      Upos = as.matrix(basis$Upos),
      U0 = as.matrix(basis$U0),
      dpos = as.numeric(basis$dpos),
      spline_ci = as.integer(1),
      t_grid = as.numeric(t_grid)
    )

    # the likelihood takes the centering constants as data
    Data$Bmean <- as.numeric(basis$Bmean)

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

    obj <- MakeADFun(data = Data, parameters = Params, random = randoms,
                     DLL = model_name, silent = TRUE)

    # Optimization
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 1e4, iter.max = 1e4))

    obj$par <- opt$par
    rep <- sdreport(obj)

    end_time <- Sys.time()

    list(ok = TRUE, time_total = as.numeric(difftime(end_time, start_time, units = "secs")))

  }, error = function(e) {
    message("    snmmTMB failed: ", conditionMessage(e))
    list(ok = FALSE, time_total = NA_real_)
  })

  out
}

#' Time assist on one dataset
#' 
#' @param dat simulated dataset
#' @return list(ok, time_total); time_total is NA on failure
fit_assist <- function(dat) {

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


  dat_fit <- as.data.frame(dat[, c("subject", "t", "y")])
  dat_fit$subject <- as.factor(dat_fit$subject)

  assign("start_val", mean(dat_fit$y), envir = .GlobalEnv)  # removed by the on.exit() above

  out <- tryCatch({

    start_time <- Sys.time()

    snm_fit <- snm(
      y ~ b1 + exp(b2) * f(t - alogit(b3)),
      func   = f(u) ~ list(~sin(2*pi*u)+cos(2*pi*u)-1, lspline(u, type = "sine0")),
      fixed  = list(b1 ~ 1),
      random = pdDiag(b1 + b2 + b3 ~ 1),
      data   = dat_fit,
      groups = ~subject,
      start  = start_val,
      verbose = FALSE
    )

    # Population curve standard errors
    bci <- intervals(snm_fit, newdata = data.frame(u = t_grid - 0.5))

    end_time <- Sys.time()

    list(ok = TRUE, time_total = as.numeric(difftime(end_time, start_time, units = "secs")))

  }, error = function(e) {
    message("    assist failed: ", conditionMessage(e))
    list(ok = FALSE, time_total = NA_real_)
  })

  out
}

#' Time snmmAGQ on one dataset
#'
#' @param dat simulated dataset
#' @param nk number of interior knots (K = nk + 4 basis functions), or NA to
#'   select it by AIC over nk_grid
#' @return list(ok, time_total); time_total is NA on failure
fit_snmmAGQ <- function(dat, nk) {


  subj_idx <- split(seq_len(nrow(dat)), dat$subject)
  y_list <- lapply(subj_idx, function(idx) dat$y[idx])
  t_list <- lapply(subj_idx, function(idx) dat$t[idx])

  nk_values <- if (is.na(nk)) nk_grid else as.integer(nk)

  out <- tryCatch({

    start_time <- Sys.time()

    # Fit at every nk in nk_values
    fit <- fit_select_nk(dat, y_list, t_list, nk_values = nk_values)

    a <- fit$a; alpha <- fit$alpha; theta <- fit$theta; cache <- fit$cache

    G <- numDeriv::hessian(function(par){
      alpha_par <- par[1]; a_par <- par[-1]
      neglogLik(cache, a_par, alpha_par, theta)
    }, c(alpha, a))
    try(solve(G), silent = TRUE)

    end_time <- Sys.time()

    list(ok = TRUE, time_total = as.numeric(difftime(end_time, start_time, units = "secs")))

  }, error = function(e) {
    message("    snmmAGQ failed: ", conditionMessage(e))
    list(ok = FALSE, time_total = NA_real_)
  })

  out
}

#' Run one replication: simulate one dataset and time every method on it
#'
#' @param cell_idx row of cell_grid
#' @param iter replication number within the cell
#' @return tibble with one row per method timed (ok = FALSE and time_total = NA
#'   where the method failed)
run_replication <- function(cell_idx, iter) {

  cr <- cell_grid[cell_idx, ]

  set.seed(seed_for(cr$n, cr$m, iter),
           kind = "Mersenne-Twister", normal.kind = "Inversion")

  sim <- simulate_dataset(cr$n, cr$m, sigma_true, D_diag, alpha_true)
  dat <- sim$dat

  # assist has no K to vary
  todo <- c("snmmTMB", "assist", "snmmAGQ")
  if (cr$sweep == "basis") todo <- setdiff(todo, "assist")

  rows <- lapply(todo, function(mth) {

    rec <- switch(
      mth,
      snmmTMB = fit_snmmTMB(dat, cr$n, basis_by_K[[as.character(cr$K)]]),
      assist  = fit_assist(dat),
      snmmAGQ = fit_snmmAGQ(dat, cr$nk),
      stop("unknown method: ", mth)
    )

    tibble(
      sweep = cr$sweep,
      cell = cr$cell,
      n = cr$n,
      m = cr$m,
      K = cr$K,
      iter = iter,
      method = mth,
      ok = rec$ok,
      time_total = rec$time_total
    )
  })

  bind_rows(rows)
}

# ==============================================================================
# 7. Task list for this job
# ==============================================================================

tasks <- expand.grid(iter = seq_len(iterations), cell = cell_grid$cell)
tasks <- tasks[, c("cell", "iter")]
tasks$job <- ((seq_len(nrow(tasks)) - 1L) %% njobs) + 1L

my_tasks <- tasks[tasks$job == itnum, c("cell", "iter")]

# Shuffle the order in which this job runs its tasks, so that the jobs do not
# all work on the same cell at the same time
set.seed(seed_base + 104729L * itnum,
         kind = "Mersenne-Twister", normal.kind = "Inversion")
my_tasks <- my_tasks[sample.int(nrow(my_tasks)), , drop = FALSE]

rownames(my_tasks) <- NULL
ntask <- nrow(my_tasks)

if (ntask == 0L) stop("Job ", itnum, " of ", njobs, " has no work to do.")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
job_file <- file.path(results_dir,
                      sprintf("%s_job%02d_of%02d.RDS", output_prefix, itnum, njobs))

# ==============================================================================
# 8. Resume
# ==============================================================================

config <- list(
  iterations = iterations,
  seed_base = seed_base,
  sigma_true = sigma_true,
  D_diag = D_diag,
  alpha_true = alpha_true,
  n_grid_subjects = n_grid_subjects,
  m_grid_observations = m_grid_observations,
  K_grid_basis = K_grid_basis,
  n_base = n_base,
  m_base = m_base,
  K_base = K_base
)

records <- vector("list", ntask)

if (file.exists(job_file)) {
  prev <- readRDS(job_file)
  same_run <- identical(as.integer(prev$njobs), njobs) &&
    identical(as.integer(prev$job), itnum) &&
    identical(prev$config, config) &&
    identical(prev$tasks, my_tasks)
  if (!same_run) {
    stop("Existing ", job_file, " was written with a different configuration ",
         "(njobs / iterations / seed_base / grid). Move it aside or fix the configuration.")
  }
  records <- prev$records
  message("Resuming ", job_file, ": ",
          sum(!vapply(records, is.null, logical(1))), "/", ntask, " already done.")
}

#' Write this job's file to a temporary file
save_job <- function() {
  tmp <- paste0(job_file, ".tmp")
  saveRDS(list(
    job = itnum,
    njobs = njobs,
    config = config,
    cell_grid = cell_grid,
    tasks = my_tasks,
    records = records,
    finished_at = Sys.time(),
    sessionInfo = utils::sessionInfo()
  ), tmp)
  invisible(file.rename(tmp, job_file))
}

# ==============================================================================
# 9. Run
# ==============================================================================

message("scalability job ", itnum, "/", njobs, ": ", ntask, " replications ",
        "(", iterations, " iterations x ", nrow(cell_grid), " cells, ",
        "seed_base = ", seed_base, ")")

job_start <- Sys.time()

for (k in seq_len(ntask)) {

  if (!is.null(records[[k]])) next

  ci <- my_tasks$cell[k]
  it <- my_tasks$iter[k]
  cr <- cell_grid[ci, ]

  res <- run_replication(ci, it)

  records[[k]] <- res

  message(sprintf("[%s] job %d/%d  task %d/%d  %-12s cell %2d  n=%3d m=%3d K=%2d  iter %2d  %d/%d ok, %.1fs  (%.1f min elapsed)",
                  format(Sys.time(), "%H:%M:%S"), itnum, njobs, k, ntask,
                  cr$sweep, cr$cell, cr$n, cr$m, cr$K, it,
                  sum(res$ok), nrow(res), sum(res$time_total, na.rm = TRUE),
                  as.numeric(difftime(Sys.time(), job_start, units = "mins"))))

  if (k %% save_every == 0L || k == ntask) save_job()
}

save_job()

n_ok <- sum(vapply(records, function(r) sum(r$ok), integer(1)))
n_all <- sum(vapply(records, nrow, integer(1)))
message(sprintf("DONE scalability job %d/%d: %d/%d fits succeeded in %.1f min -> %s",
                itnum, njobs, n_ok, n_all,
                as.numeric(difftime(Sys.time(), job_start, units = "mins")), job_file))
