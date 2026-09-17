# ------------------------------------------------------------------------------
# snmm_monotonesimulation.R
#
# This script runs the illustration of the monotonicity penalty of Online
# Appendix D: a penalized B-spline is fitted to data simulated from the
# non-monotone function y = x^3 - x, under increasing values of the penalty
# weight lambda_c, and the fitted curves with their simultaneous confidence
# bands are plotted against the truth.
#
# Run from this directory:
#
#   Rscript snmm_monotonesimulation.R
#
# Writes Supplementary_Figure4.pdf   fitted curve and simultaneous band, one panel per lambda_c
#
# RUNTIME: A few seconds on a laptop, plus a one-off TMB compilation of about a
# minute the first time.
#
# Structure:
#   1. Settings
#   2. Simulate one dataset
#   3. Basis, penalty and centering constants
#   4. Fit snmmTMB for each lambda_c
#   5. Plot
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm()
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(TMB)      # snmmTMB
  library(splines)  # splineDesign()
})

set.seed(0)


# ==============================================================================
# 1. Settings
#
# Data-generating model (Online Appendix D):
#
#   y_ij = x_j^3 - x_j + eps_ij
#   x_j  = xmin + (j-1)/(m-1) * (xmax - xmin) in [-2,2]
#   eps_ij ~ N(0, sigma^2)
#
# There are no subject-level random effects: the n subjects share the same
# curve, so the model reduces to a penalized spline regression with a
# monotonicity penalty. The curve is decreasing on (-1/sqrt(3), 1/sqrt(3)), so
# an unpenalized fit (lambda_c = 0) is correctly non-monotone, and the fitted
# curve becomes less decreasing as lambda_c grows.
# ==============================================================================

n      <- 20                  # subjects
m      <- 20                  # observations per subject
sigma  <- 0.4                 # residual sd
xmin   <- -2                  # covariate range
xmax   <-  2

K      <- 10                  # snmmTMB basis functions
degree <- 3                   # cubic B-splines

lambda_c_values <- c(0, 0.1, 1, 10)   # monotonicity penalty weights, one panel each
Nsims  <- 10000               # draws for the simultaneous-band critical value

x_grid <- seq(xmin, xmax, length.out = 50)   # grid the curve is reported on, and
                                             # at which f' is penalized

# true curve
mu_true <- x_grid^3 - x_grid


# ==============================================================================
# 2. Simulate one dataset
# ==============================================================================

dat_list <- vector("list", n)
for (i in 1:n) {
  x     <- seq(xmin, xmax, length.out = m)
  mu_ij <- x^3 - x
  dat_list[[i]] <- tibble(subject = i,
                          x = x,
                          y = mu_ij + rnorm(m, mean = 0, sd = sigma))
}
dat <- bind_rows(dat_list)
dat$subject <- as.integer(dat$subject)


# ==============================================================================
# 3. Basis, penalty and centering constants (snmmTMB)
#
# Cubic B-spline basis with K functions. The covariate is scaled to [0,1] and
# the knots are fixed on that interval.
# ==============================================================================

scale_x <- function(x) (x - xmin) / (xmax - xmin)

Kint <- K - 1   # the last basis column is dropped so the smooth cannot absorb beta1
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

# Centering constants: column means of the basis over a fixed grid spanning the knot range
edges  <- seq(0, 1, length.out = 1001)
u_grid <- (head(edges, -1) + tail(edges, -1)) / 2
Bmean  <- colMeans(splineDesign(knot_vec, u_grid, ord = degree + 1,
                                outer.ok = TRUE)[, 1:Kint, drop = FALSE])

# Basis derivatives on x_grid, for the monotonicity penalty. The derivative is
# taken with respect to the scaled covariate, which only rescales lambda_c by
# the constant factor xmax - xmin.
X_deriv <- splineDesign(knot_vec, scale_x(x_grid), ord = degree + 1,
                        outer.ok = TRUE, derivs = 1)[, 1:Kint, drop = FALSE]

# compile once, then load
cpp <- "../../src/snmmTMB_derivpenalty.cpp"
dll <- dynlib(sub("\\.cpp$", "", cpp))
if (!file.exists(dll)) TMB::compile(cpp)
dyn.load(dll)


# ==============================================================================
# 4. Fit snmmTMB for each lambda_c
# ==============================================================================

Data <- list(
  y = dat$y, x = scale_x(dat$x),
  knots = as.numeric(knot_vec), degree = as.integer(degree), K = as.integer(K),
  Bmean = as.numeric(Bmean),
  Upos = as.matrix(Upos), U0 = as.matrix(U0), dpos = as.numeric(dpos),
  X_deriv = as.matrix(X_deriv), lambda_c = 0,
  spline_ci = as.integer(1), x_grid = scale_x(x_grid)
)

Params <- list(beta1 = 0, log_sigma = log(1),
               vpos = rep(0, length(dpos)), gamma0 = rep(0, ncol(U0)),
               log_lambda = log(1))

fit_list <- vector("list", length(lambda_c_values))
for (k in seq_along(lambda_c_values)) {
  Data$lambda_c <- lambda_c_values[k]

  obj <- MakeADFun(Data, Params, random = "vpos",
                   DLL = "snmmTMB_derivpenalty", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(eval.max = 1e4, iter.max = 1e4))
  obj$par <- opt$par
  rep <- sdreport(obj)

  # fitted curve with pointwise and simultaneous bands
  idx <- grep("h_grid", names(rep$value))
  tmb_mu <- rep$value[idx]
  Cov_mu <- rep$cov[idx, idx]
  Cov_mu <- (Cov_mu + t(Cov_mu)) / 2
  tmb_se <- sqrt(diag(Cov_mu))

  mu_sims  <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu)), Sigma = Cov_mu)
  tmb_crit <- quantile(apply(abs(sweep(mu_sims, 1, tmb_se, FUN = "/")), 1, max), 0.95)

  fit_list[[k]] <- tibble(lambda_c = lambda_c_values[k], x = x_grid,
                          mu = tmb_mu, truth = mu_true,
                          lwr_p = tmb_mu - 1.96 * tmb_se,
                          upr_p = tmb_mu + 1.96 * tmb_se,
                          lwr_s = tmb_mu - tmb_crit * tmb_se,
                          upr_s = tmb_mu + tmb_crit * tmb_se)

  cat("lambda_c =", lambda_c_values[k],
      ": convergence", opt$convergence, ", nll", round(opt$objective, 3), "\n")
}
fit <- bind_rows(fit_list)


# ==============================================================================
# 5. Plot
# ==============================================================================

fit$lambda_lab <- factor(paste0("lambda[c] == ", fit$lambda_c),
                         levels = paste0("lambda[c] == ", lambda_c_values))

supplementary_figure4 <- ggplot(fit, aes(x)) +
  geom_point(data = dat, aes(x, y), colour = "grey70", size = 1, alpha = 0.35) +
  geom_line(aes(y = truth), colour = "black", linetype = "dashed", linewidth = 0.9) +
  geom_ribbon(aes(ymin = lwr_s, ymax = upr_s), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = mu), colour = "blue", linewidth = 1.1) +
  facet_wrap(~lambda_lab, nrow = 1, labeller = label_parsed) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "x", y = "y") +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.8),
        panel.spacing = unit(1.2, "lines"),
        strip.background = element_rect(fill = "grey95", colour = "grey40"),
        strip.text = element_text(size = 16),
        panel.grid.minor = element_blank())

ggsave("Supplementary_Figure4.pdf", supplementary_figure4, width = 12, height = 6, units = "in")
