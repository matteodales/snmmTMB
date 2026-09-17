# ------------------------------------------------------------------------------
# bellcurve_edf_diagnostic.R
#
# This script runs the effective degrees of freedom (EDF) comparison of Online Appendix C
#
# Run from this directory:
#
#   Rscript bellcurve_edf_diagnostic.R
#
# Writes results/bellcurve_edf_diagnostic_results.RDS   per-replication EDFs
#        Supplementary_Figure3.pdf   
#
# Set run_simulation <- FALSE to skip the fits and redraw the figure from the
# saved RDS.
#
# RUNTIME: about 25 minutes on a laptop for the 50 replications
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm()
  library(dplyr)
  library(tidyr)    # pivot_longer()
  library(tibble)
  library(ggplot2)
  library(TMB)      # snmmTMB
  library(nlme)     # starting values
  library(splines)  # splineDesign()
  library(assist)   # snm()
})

# startingvalues(), bspline.snm(), fit_select_nk(), nk_grid, knots.equispaced()
source("../../src/snmmAGQ_functions_bellcurve.R")

run_simulation <- FALSE        # FALSE: read results_file and only draw the figure
results_file <- "results/bellcurve_edf_diagnostic_results.RDS"


# ==============================================================================
# 1. Settings
#
# Data-generating model (manuscript Section 3):
#
#   y_ij = alpha + b1_i + exp(-0.5 * (t_j - b2_i)^2) + eps_ij
#   t_j  = m equally spaced points in [-2,2]
#   b_i  = (b1_i, b2_i) ~ N(0, sigma^2 * D),  D = [2 1; 1 2]
#   eps_ij ~ N(0, sigma^2)
#
# The values below are setting 5 of the manuscript grid
# ==============================================================================

n      <- 10                  # subjects
m      <- 10                  # observations per subject
sigma  <- 0.4                 # residual sd
D_diag <- c(2, 2)             # diagonal of D
rho    <- 0.5                 # correlation between b1 and b2
alpha  <- 1                   # population intercept

K      <- 15                  # snmmTMB basis functions
degree <- 3                   # cubic B-splines
                              # snmmAGQ: interior knots selected by AIC over
                              # nk_grid (0:6, set in the sourced functions file)

iterations <- 50              # replications

t_grid <- seq(-2, 2, length.out = 50)   # grid the curves are reported on
ngrid  <- length(t_grid)

cov_b12 <- rho * sqrt(D_diag[1] * D_diag[2])
Sigma_b <- sigma^2 * matrix(c(D_diag[1], cov_b12,
                              cov_b12,   D_diag[2]), nrow = 2, byrow = TRUE)


# ==============================================================================
# 2. Basis, penalty and centering constants (snmmTMB)
# ==============================================================================

Kint <- K - 1
n_internal <- K - (degree + 1)

internal_knots <- seq(0, 1, length.out = n_internal + 2)
dx_left  <- diff(internal_knots)[1]
dx_right <- diff(internal_knots)[length(internal_knots) - 1]
knot_vec <- c(internal_knots[1] - dx_left * (degree:1),
              internal_knots,
              internal_knots[length(internal_knots)] + dx_right * (1:degree))

D <- diff(diag(K), differences = 2)[, 1:Kint]
P <- t(D) %*% D
P <- P * qr(P)$rank / sum(diag(P))

ev <- eigen(P, symmetric = TRUE)
pos_idx  <- which(ev$values > 1e-12)
zero_idx <- which(ev$values <= 1e-12)
Upos <- ev$vectors[, pos_idx, drop = FALSE]
dpos <- ev$values[pos_idx]
U0   <- ev$vectors[, zero_idx, drop = FALSE]

edges  <- seq(0, 1, length.out = 1001)
u_grid <- (head(edges, -1) + tail(edges, -1)) / 2
Bmean  <- colMeans(splineDesign(knot_vec, u_grid, ord = degree + 1,
                                outer.ok = TRUE)[, 1:Kint, drop = FALSE])


edf_ceiling <- length(dpos) + ncol(U0)


# ==============================================================================
# 3. EDF helpers
# ==============================================================================

#' EDF of snmmTMB's shape spline from the joint Hessian
#'
#' `obj` must be evaluated at its optimum (obj$par <- opt$par; obj$fn(opt$par))
#' so that obj$env$last.par.best is the joint mode. Errors if the optimizer did
#' not reach a mode, the Hessian is not positive definite, or the EDF falls
#' outside its admissible range
#'
#' @param obj fitted TMB ADFun object with random = c("b1", "b2", "vpos")
#' @param lambda estimated smoothing parameter, exp(log_lambda)
#' @return the EDF

compute_edf_tmb <- function(obj, lambda) {
  random_names <- names(obj$env$last.par.best)[obj$env$random]
  vpos_pos <- which(random_names == "vpos")
  p <- length(vpos_pos)

  H <- as.matrix(obj$env$spHess(obj$env$last.par.best, random = TRUE))
  min_eig <- min(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(min_eig) || min_eig <= 0) {
    stop(sprintf("joint Hessian not positive definite (min eigenvalue %.3g)", min_eig))
  }

  A_diag <- diag(chol2inv(chol(H)))
  edf <- p - sum(A_diag[vpos_pos] * lambda * dpos) + ncol(U0)

  if (edf < ncol(U0) || edf > p + ncol(U0)) {
    stop(sprintf("EDF %.3g outside [%d, %d]", edf, ncol(U0), p + ncol(U0)))
  }
  edf
}

#' EDF of snmmAGQ's shape spline: the rank of its centered B-spline basis
#'
#' @param fit a fit_select_nk() fit
#' @param dat the dataset it was fitted to
#' @return the EDF
compute_edf_agq <- function(fit, dat) {
  u_all <- dat$t - fit$modes[, 2][dat$subject]
  B <- splineDesign(fit$knotseq, u_all, 4, outer.ok = TRUE)
  qr(sweep(B, 2, colMeans(B), FUN = "-"))$rank
}


# ==============================================================================
# 4. Fit functions
# ==============================================================================

#' Simulate one bell-curve dataset
#' @return tibble(subject, t, y)
simulate_dataset <- function() {
  b_mat <- mvrnorm(n = n, mu = rep(0, 2), Sigma = Sigma_b)
  colnames(b_mat) <- c("b1", "b2")

  dat_list <- vector("list", n)
  for (i in 1:n) {
    t     <- seq(-2, 2, length.out = m)
    mu_ij <- alpha + b_mat[i, "b1"] + exp(-0.5 * (t - b_mat[i, "b2"])^2)
    dat_list[[i]] <- tibble(subject = i,
                            t = t,
                            y = mu_ij + rnorm(m, mean = 0, sd = sigma))
  }
  dat <- bind_rows(dat_list)
  dat$subject <- as.integer(dat$subject)
  dat
}

#' Fit snmmTMB
#' @param dat tibble(subject, t, y)
#' @return the EDF, or NULL
fit_snmmTMB <- function(dat) {

  on.exit(if (exists("spline_start", envir = .GlobalEnv, inherits = FALSE)) {
    rm("spline_start", envir = .GlobalEnv)
  }, add = TRUE)

  tryCatch({

    t_min <- min(dat$t)
    t_max <- max(dat$t)

    Data0 <- list(
      y = dat$y, x = (dat$t - t_min) / (t_max - t_min),
      knots = as.numeric(knot_vec), degree = as.integer(degree), K = as.integer(K),
      Upos = as.matrix(Upos), U0 = as.matrix(U0), dpos = as.numeric(dpos),
      spline_ci = as.integer(1), x_grid = seq(0, 1, length.out = ngrid)
    )

    Params0 <- list(beta1 = 0, log_sigma = log(1),
                    vpos = rep(0, length(dpos)), gamma0 = rep(0, ncol(U0)),
                    log_lambda = log(1))

    obj0 <- MakeADFun(Data0, Params0, random = "vpos",
                      DLL = "starting_points", silent = TRUE)
    opt0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
    obj0$par <- opt0$par
    rep0 <- sdreport(obj0)

    c_hat <- rep0$value[grep("c", names(rep0$value))]
    m_var <- as.numeric(obj0$report()$m)

    # must be visible from .GlobalEnv for nlme() to find it from inside a function
    spline_start <- function(t) {
      sweep(as.matrix(splineDesign(knot_vec, (t - t_min) / (t_max - t_min), 4, outer.ok = TRUE))[, 1:Kint, drop = FALSE],
            2, m_var, FUN = "-") %*% c_hat
    }
    assign("spline_start", spline_start, envir = .GlobalEnv)  # removed by the on.exit() above

    # general (correlated) random-effect covariance
    nlmeobj <- nlme(y ~ b1 + spline_start(t - b2),
                    fixed   = list(b1 ~ 1),
                    random  = list(b1 + b2 ~ 1),
                    data    = dat,
                    groups  = ~subject,
                    start   = mean(dat$y),
                    control = list(returnObject = TRUE, tolerance = .01))

    # snmmTMB fit
    Data <- list(
      y = dat$y, t = dat$t,
      group = as.integer(dat$subject), nGroup = as.integer(n),
      subj_flag = as.integer(c(1, 1, rep(0, n - 2))),
      knots = as.numeric(knot_vec), degree = as.integer(degree), K = as.integer(K),
      Bmean = as.numeric(Bmean),
      Upos = as.matrix(Upos), U0 = as.matrix(U0), dpos = as.numeric(dpos),
      spline_ci = as.integer(1), t_grid = as.numeric(t_grid)
    )

    parList0 <- obj0$env$parList()
    Params <- list(
      alpha = nlmeobj$coefficients$fixed["b1"],
      b1 = as.numeric(nlmeobj$coefficients$random$subject[, "b1"]),
      b2 = as.numeric(nlmeobj$coefficients$random$subject[, "b2"]),
      log_sd_b1 = log(as.numeric(VarCorr(nlmeobj)["b1", "StdDev"])),
      log_sd_b2 = log(as.numeric(VarCorr(nlmeobj)["b2", "StdDev"])),
      log_sigma = log(as.numeric(VarCorr(nlmeobj)["Residual", "StdDev"])),
      transf_rho = 0,
      vpos = as.numeric(parList0$vpos), gamma0 = as.numeric(parList0$gamma0),
      log_lambda = as.numeric(parList0$log_lambda)
    )

    obj <- MakeADFun(Data, Params, random = c("b1", "b2", "vpos"),
                     DLL = "snmmTMB_likelihood_bellcurve_simulation", silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 1e4, iter.max = 1e4))
    obj$par <- opt$par
    obj$fn(opt$par)

    compute_edf_tmb(obj, exp(opt$par["log_lambda"]))

  }, error = function(e) {
    message("  snmmTMB failed: ", conditionMessage(e))
    NULL
  })
}

#' Fit assist
#' @param dat tibble(subject, t, y)
#' @return the EDF, or NULL
fit_assist <- function(dat) {

  on.exit(if (exists("start_val", envir = .GlobalEnv, inherits = FALSE)) {
    rm("start_val", envir = .GlobalEnv)
  }, add = TRUE)

  on.exit({
    if (exists("f", envir = .GlobalEnv, inherits = FALSE) &&
        is.function(get("f", envir = .GlobalEnv))) {
      rm("f", envir = .GlobalEnv)
    }
  }, add = TRUE)


  dat_fit <- as.data.frame(dat[, c("subject", "t", "y")])
  dat_fit$subject <- as.factor(dat_fit$subject)
  dat_fit$t_scaled <- dat_fit$t / 10 + 1 / 2

  start_val <- mean(dat_fit$y)
  assign("start_val", start_val, envir = .GlobalEnv)  # removed by the on.exit() above

  tryCatch({

    snm_fit <- snm(y ~ b1 + f(t_scaled - b2),
                   func    = f(u) ~ list(~u - 1, cubic(u)),
                   fixed   = list(b1 ~ 1),
                   random  = pdSymm(b1 + b2 ~ 1, value = Sigma_b),
                   data    = dat_fit,
                   groups  = ~subject,
                   start   = start_val,
                   verbose = FALSE)

    as.numeric(snm_fit$forCI$rkpk.obj$df)

  }, error = function(e) {
    message("  assist failed: ", conditionMessage(e))
    NULL
  })
}

#' Fit snmmAGQ
#' @param dat tibble(subject, t, y)
#' @return list(edf, nk) with the EDF and the AIC-selected knot count, or NULL
fit_snmmAGQ <- function(dat) {

  subj_idx <- split(seq_len(nrow(dat)), dat$subject)
  y_list <- lapply(subj_idx, function(idx) dat$y[idx])
  t_list <- lapply(subj_idx, function(idx) dat$t[idx])

  tryCatch({
    fit <- fit_select_nk(dat, y_list, t_list)
    list(edf = compute_edf_agq(fit, dat), nk = fit$nk_selected)
  }, error = function(e) {
    message("  snmmAGQ failed: ", conditionMessage(e))
    NULL
  })
}


# ==============================================================================
# 5. Monte Carlo loop
# ==============================================================================

if (run_simulation) {

  # compile once, then load
  for (cpp in c("../../src/starting_points.cpp",
                "../../src/snmmTMB_likelihood_bellcurve_simulation.cpp")) {
    dll <- dynlib(sub("\\.cpp$", "", cpp))
    if (!file.exists(dll)) TMB::compile(cpp)
    dyn.load(dll)
  }

  dir.create(dirname(results_file), showWarnings = FALSE)

  edf <- tibble(iter = seq_len(iterations),
                snmmTMB = NA_real_, assist = NA_real_, snmmAGQ = NA_real_,
                nk_agq = NA_integer_)

  set.seed(123)

  for (iter in seq_len(iterations)) {

    dat <- simulate_dataset()

    edf_tmb <- fit_snmmTMB(dat)
    edf_assist <- fit_assist(dat)
    agq <- fit_snmmAGQ(dat)

    if (!is.null(edf_tmb))    edf$snmmTMB[iter] <- edf_tmb
    if (!is.null(edf_assist)) edf$assist[iter] <- edf_assist
    if (!is.null(agq)) {
      edf$snmmAGQ[iter] <- agq$edf
      edf$nk_agq[iter] <- agq$nk
    }

    message(sprintf("[%s] iter %2d/%d  snmmTMB %s  assist %s  snmmAGQ %s",
                    format(Sys.time(), "%H:%M:%S"), iter, iterations,
                    format(edf$snmmTMB[iter], digits = 3),
                    format(edf$assist[iter], digits = 3),
                    format(edf$snmmAGQ[iter], digits = 3)))

    if (iter %% 10 == 0 || iter == iterations) {
      saveRDS(list(setting = list(n = n, m = m, sigma = sigma, D_diag = D_diag,
                                  rho = rho, K = K, nk_grid = nk_grid,
                                  iterations = iterations),
                   edf = edf),
              results_file)
    }
  }
}


# ==============================================================================
# 6. Plot
# ==============================================================================

results <- readRDS(results_file)

edf_long <- results$edf %>%
  filter(!is.na(snmmTMB), !is.na(assist), !is.na(snmmAGQ)) %>%
  pivot_longer(c(snmmTMB, assist, snmmAGQ), names_to = "method", values_to = "edf") %>%
  mutate(method = factor(method, levels = c("snmmTMB", "assist", "snmmAGQ")))

message(sprintf("%d of %d replications retained (all three methods succeeded)",
                n_distinct(edf_long$iter), results$setting$iterations))

method_colors <- c(snmmTMB = "#0046FA", assist = "#FF6242", snmmAGQ = "#A1E600")

supplementary_figure3 <- ggplot(edf_long, aes(x = method, y = edf, fill = method)) +
  geom_boxplot(width = 0.35, linewidth = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 1, alpha = 0.35) +
  scale_fill_manual(values = method_colors) +
  scale_y_log10(limits = c(2.5, NA), minor_breaks = NULL,
                breaks = c(5, 10, 20, 30, 40, 50, 60, 70, 80),
                labels = c("5", "10", "20", "", "40", "", "60", "", "80")) +
  labs(x = NULL, y = "EDF (log scale)") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

ggsave("Supplementary_Figure3.pdf", supplementary_figure3, width = 7, height = 4.5, units = "in")
