# ------------------------------------------------------------------------------
# sinecurve_example.R
#
# This script runs one replication of the sine-curve simulation of manuscript Section 3,
# fitted with all three methods, with the estimated population and subject curves and
# their confidence bands plotted against the truth.
#
# The other "server" scripts in this folder are meant to run this same simulation over 16 settings
# and 500 replications, on a HPC cluster.
#
# RUNTIME: A few seconds per method on a laptop (changes with sample size),
# plus a one-off TMB compilation of about a minute the first time.
#
# Structure:
#   1. Settings
#   2. Simulate one dataset
#   3. Basis, penalty and centering constants
#   4. Fit snmmTMB
#   5. Fit assist
#   6. Fit snmmAGQ
#   7. Plots
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm()
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(TMB)      # snmmTMB
  library(nlme)     # starting values
  library(splines)  # splineDesign()
  library(assist)   # snm()
})

# startingvalues(), bspline.snm(), neglogLik(), fit_select_nk(), nk_grid,
# subj_curve_fun(), make_pd()
source("../src/snmmAGQ_functions_sinecurve.R")

set.seed(1)


# ==============================================================================
# 1. Settings
#
# Data-generating model (manuscript Section 3):
#
#   y_ij = alpha + b1_i + 2*exp(b2_i) * sin(2*pi*(t_j - expit(b3_i))) + eps_ij
#   t_j  = (j-1)/(m-1) in [0,1]
#   b_i  = (b1_i, b2_i, b3_i) ~ N(0, sigma^2 * diag(D_diag))
#   eps_ij ~ N(0, sigma^2)
#
# b1_i shifts the curve vertically, b2_i scales its amplitude, b3_i shifts its
# phase. The values below are setting 1 of the manuscript grid.
# ==============================================================================

n      <- 10                  # subjects
m      <- 10                  # observations per subject
sigma  <- 1                   # residual sd
D_diag <- c(1, 0.25, 0.16)    # random-effect variance multipliers
alpha  <- 5                   # population intercept

K      <- 15                  # snmmTMB basis functions
degree <- 3                   # cubic B-splines
                              # snmmAGQ: interior knots selected by AIC over
                              # nk_grid (0:6, set in the sourced functions file)

subjects_to_plot <- c(1, 2)   # subjects to report subject-level curves for
Nsims  <- 10000               # draws for the simultaneous-band critical value

t_grid <- seq(0, 1, length.out = 50)   # grid the curves are reported on
ngrid  <- length(t_grid)

# true population curve: at b = 0 the amplitude is 2*exp(0) = 2 and the phase
# expit(0) = 0.5
mu_true <- alpha + 2 * sin(2 * pi * (t_grid - 0.5))


# ==============================================================================
# 2. Simulate one dataset
# ==============================================================================

Sigma_b <- sigma^2 * diag(D_diag)
b_mat <- mvrnorm(n = n, mu = rep(0, 3), Sigma = Sigma_b)
colnames(b_mat) <- c("b1", "b2", "b3")

dat_list <- vector("list", n)
for (i in 1:n) {
  t     <- seq(0, 1, length.out = m)
  amp   <- 2 * exp(b_mat[i, "b2"])
  phase <- exp(b_mat[i, "b3"]) / (1 + exp(b_mat[i, "b3"]))
  mu_ij <- alpha + b_mat[i, "b1"] + amp * sin(2 * pi * (t - phase))
  dat_list[[i]] <- tibble(subject = i,
                          t = t,
                          y = mu_ij + rnorm(m, mean = 0, sd = sigma))
}
dat <- bind_rows(dat_list)
dat$subject <- as.integer(dat$subject)

# the simulated data with the true population curve on top
ggplot(dat, aes(t, y)) +
  geom_line(aes(group = subject), colour = "grey60", alpha = 0.7) +
  geom_point(colour = "grey40", size = 1) +
  geom_line(data = tibble(t = t_grid, y = mu_true), linewidth = 1) +
  labs(title = "Simulated data", subtitle = "black = true population curve") +
  theme_minimal()


# ==============================================================================
# 3. Basis, penalty and centering constants (snmmTMB)
#
# Cubic B-spline basis with K functions. The spline argument is
# t_ij - expit(b3_i), which lies in [-1,1] whatever b3 is, so the knots are
# fixed on that interval (Section 2.3.1, first case).
# ==============================================================================

Kint <- K - 1   # the last basis column is dropped so the smooth cannot absorb alpha
n_internal <- K - (degree + 1)

internal_knots <- seq(-1, 1, length.out = n_internal + 2)
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
edges  <- seq(-1, 1, length.out = 1001)
u_grid <- (head(edges, -1) + tail(edges, -1)) / 2
Bmean  <- colMeans(splineDesign(knot_vec, u_grid, ord = degree + 1,
                                outer.ok = TRUE)[, 1:Kint, drop = FALSE])

# compile once, then load
for (cpp in c("../src/starting_points.cpp",
              "../src/snmmTMB_likelihood_sinecurve_simulation.cpp")) {
  dll <- dynlib(sub("\\.cpp$", "", cpp))
  if (!file.exists(dll)) TMB::compile(cpp)
  dyn.load(dll)
}


# ==============================================================================
# 4. Fit snmmTMB
# ==============================================================================

Data0 <- list(
  y = dat$y, x = dat$t - 0.5,
  knots = as.numeric(knot_vec), degree = as.integer(degree), K = as.integer(K),
  Upos = as.matrix(Upos), U0 = as.matrix(U0), dpos = as.numeric(dpos),
  spline_ci = as.integer(1), x_grid = t_grid - 0.5
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

spline_start <- function(t) {
  sweep(as.matrix(splineDesign(knot_vec, t, 4, outer.ok = TRUE))[, 1:Kint, drop = FALSE],
        2, m_var, FUN = "-") %*% c_hat
}

nlmeobj <- nlme(y ~ b1 + exp(b2) * spline_start(t - exp(b3) / (1 + exp(b3))),
                fixed   = list(b1 ~ 1),
                random  = pdDiag(b1 + b2 + b3 ~ 1),
                data    = dat,
                groups  = ~subject,
                start   = mean(dat$y),
                control = list(returnObject = TRUE, tolerance = .01))

# snmmTMB fit
Data <- list(
  y = dat$y, t = dat$t,
  group = as.integer(dat$subject), nGroup = as.integer(n),
  subj_flag = as.integer(seq_len(n) %in% subjects_to_plot),
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
  b3 = as.numeric(nlmeobj$coefficients$random$subject[, "b3"]),
  log_sd_b1 = log(as.numeric(VarCorr(nlmeobj)["b1", "StdDev"])),
  log_sd_b2 = log(as.numeric(VarCorr(nlmeobj)["b2", "StdDev"])),
  log_sd_b3 = log(as.numeric(VarCorr(nlmeobj)["b3", "StdDev"])),
  log_sigma = log(nlmeobj$sigma),
  vpos = as.numeric(parList0$vpos), gamma0 = as.numeric(parList0$gamma0),
  log_lambda = as.numeric(obj0$par["log_lambda"])
)

obj <- MakeADFun(Data, Params, random = c("b1", "b2", "b3", "vpos"),
                 DLL = "snmmTMB_likelihood_sinecurve_simulation", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 1e4, iter.max = 1e4))
obj$par <- opt$par
rep <- sdreport(obj)

# population curve with pointwise and simultaneous bands
idx <- grep("h_grid", names(rep$value))
tmb_mu <- rep$value[idx]
Cov_mu <- rep$cov[idx, idx]
Cov_mu <- (Cov_mu + t(Cov_mu)) / 2
tmb_se <- sqrt(diag(Cov_mu))

mu_sims  <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu)), Sigma = Cov_mu)
tmb_crit <- quantile(apply(abs(sweep(mu_sims, 1, tmb_se, FUN = "/")), 1, max), 0.95)


# ==============================================================================
# 5. Fit assist
# ==============================================================================

dat_fit <- as.data.frame(dat[, c("subject", "t", "y")])
dat_fit$subject <- as.factor(dat_fit$subject)

start_val <- mean(dat_fit$y)   # snm() resolves `start` on the search path

snm_fit <- snm(y ~ b1 + exp(b2) * f(t - alogit(b3)),
               func    = f(u) ~ list(~sin(2*pi*u) + cos(2*pi*u) - 1,
                                     lspline(u, type = "sine0")),
               fixed   = list(b1 ~ 1),
               random  = pdDiag(b1 + b2 + b3 ~ 1),
               data    = dat_fit,
               groups  = ~subject,
               start   = start_val,
               verbose = FALSE)

# assist pointwise band
bci <- intervals(snm_fit, newdata = data.frame(u = t_grid - 0.5))
assist_mu <- as.numeric(as.matrix(bci$fit)[, 1]) + snm_fit$coefficients$fixed
assist_se_point <- as.numeric(as.matrix(bci$pstd))

# --- reconstructed simultaneous band ------------------------------------------
#
# We obtain the critical value by simulation as the other two methods do.
#
# All required quantities are taken from snm_fit$forCI object, which is used for assist's own pstd

u_new <- t_grid - 0.5
fc <- snm_fit$forCI
ro <- fc$rkpk.obj

u_obs <- fc$data[[1]]$u        # t_ij - alogit(b3_hat_i): shifted observations
d1    <- fc$delta1[, 1]        # exp(b2_hat_i): each row's amplitude multiplier
Vinv  <- crossprod(fc$weight)  # weight W satisfies W V W' = I, so this is V^-1
lam   <- 10^ro$nlaht           # penalty coefficient
sig2  <- ro$varht              # residual variance

# forCI$q[[1]] is the kernel already multiplied by the amplitude d1, so we recover the plain kernel K(u_i, u_j).
K_obs <- fc$q[[1]] / tcrossprod(d1)
K_new <- lspline(u_new, u_obs, type = "sine0")   # K(u_new, u_obs), 50 x nobs

# K is rank deficient, so we take the eigendecomposition K = U L U', drop the numerically
# null directions, and reparameterise c = U L^(-1/2) beta_c.

eg    <- eigen((K_obs + t(K_obs)) / 2, symmetric = TRUE)
keep  <- eg$values > 1e-12 * eg$values[1]
Uk    <- eg$vectors[, keep, drop = FALSE]
Lk    <- eg$values[keep]
nkeep <- length(Lk)

Xt     <- cbind(1, d1 * cbind(sin(2 * pi * u_obs), cos(2 * pi * u_obs)))
Xt_new <- cbind(1, cbind(sin(2 * pi * u_new), cos(2 * pi * u_new)))
Xc     <- d1 * (Uk %*% diag(sqrt(Lk), nrow = nkeep))
Xc_new <- K_new %*% (Uk %*% diag(1 / sqrt(Lk), nrow = nkeep))

X     <- cbind(Xt, Xc)        # nobs x (3 + nkeep)
X_new <- cbind(Xt_new, Xc_new)  # 50 x (3 + nkeep)
Pen   <- diag(c(rep(0, ncol(Xt)), rep(1, nkeep)))   # no penalty on the null space

XtVX      <- t(X) %*% Vinv %*% X
Ainv      <- solve(XtVX + lam * Pen)
Cov_beta  <- sig2 * Ainv

Cov_curve <- X_new %*% Cov_beta %*% t(X_new)
Cov_curve <- (Cov_curve + t(Cov_curve)) / 2
assist_se_sim <- sqrt(pmax(diag(Cov_curve), 0))

sims <- mvrnorm(Nsims, mu = rep(0, length(assist_se_sim)), Sigma = Cov_curve)
assist_crit <- quantile(apply(abs(sweep(sims, 2, assist_se_sim, FUN = "/")), 1, max), 0.95)


# ==============================================================================
# 6. Fit snmmAGQ
# ==============================================================================

subj_idx <- split(seq_len(nrow(dat)), dat$subject)
y_list <- lapply(subj_idx, function(idx) dat$y[idx])
t_list <- lapply(subj_idx, function(idx) dat$t[idx])

fit_agq <- fit_select_nk(dat, y_list, t_list)
cat("snmmAGQ: AIC over nk =", nk_grid, "->", round(fit_agq$nk_aic, 1),
    "; selected nk =", fit_agq$nk_selected,
    "(", length(fit_agq$a), "basis functions )\n")

a         <- fit_agq$a
alpha_agq <- fit_agq$alpha
knotseq   <- fit_agq$knotseq
col_means <- fit_agq$col_means
modes     <- fit_agq$modes
hessians  <- fit_agq$hessians

# Hessian of the log-likelihood in the fixed effects, by numerical differentiation
G <- numDeriv::hessian(function(par) neglogLik(fit_agq$cache, par[-1], par[1], fit_agq$theta),
                       c(alpha_agq, a))
Ginv <- solve(G)

X_pop <- cbind(1, sweep(splineDesign(knotseq, t_grid - 0.5, 4, outer.ok = TRUE),
                        2, col_means, FUN = "-"))
agq_mu  <- as.numeric(X_pop %*% c(alpha_agq, a))
Cov_agq <- make_pd(X_pop %*% Ginv %*% t(X_pop))
agq_se  <- sqrt(diag(Cov_agq))

mu_sims  <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_agq)), Sigma = Cov_agq)
agq_crit <- quantile(apply(abs(sweep(mu_sims, 2, agq_se, FUN = "/")), 1, max), 0.95)


# ==============================================================================
# 7. Plots
# ==============================================================================

pop <- bind_rows(
  tibble(method = "snmmTMB", t = t_grid, mu = tmb_mu,
         lwr_p = tmb_mu - 1.96 * tmb_se,
         upr_p = tmb_mu + 1.96 * tmb_se,
         lwr_s = tmb_mu - tmb_crit * tmb_se,
         upr_s = tmb_mu + tmb_crit * tmb_se),
  tibble(method = "assist", t = t_grid, mu = assist_mu,
         lwr_p = assist_mu - 1.96 * assist_se_point,
         upr_p = assist_mu + 1.96 * assist_se_point,
         lwr_s = assist_mu - assist_crit * assist_se_sim,
         upr_s = assist_mu + assist_crit * assist_se_sim),
  tibble(method = "snmmAGQ", t = t_grid, mu = agq_mu,
         lwr_p = agq_mu - 1.96 * agq_se,
         upr_p = agq_mu + 1.96 * agq_se,
         lwr_s = agq_mu - agq_crit * agq_se,
         upr_s = agq_mu + agq_crit * agq_se)
)
pop$method <- factor(pop$method, levels = c("snmmTMB", "assist", "snmmAGQ"))
truth <- tibble(t = t_grid, mu = mu_true)

# one panel per method
ggplot(pop, aes(t)) +
  geom_ribbon(aes(ymin = lwr_s, ymax = upr_s), fill = "grey75", alpha = 0.7) +
  geom_ribbon(aes(ymin = lwr_p, ymax = upr_p), fill = "grey45", alpha = 0.5) +
  geom_point(data = dat, aes(t, y), colour = "grey50", size = 0.7, alpha = 0.5) +
  geom_line(aes(y = mu), linewidth = 0.9) +
  geom_line(data = truth, aes(y = mu), linetype = "dashed") +
  facet_wrap(~method) +
  labs(title = "Population curve",
       y = "y") +
  theme_minimal()

# three methods together, simultaneous bands only
ggplot(pop, aes(t, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lwr_s, ymax = upr_s), alpha = 0.15, colour = NA) +
  geom_line(aes(y = mu, linetype = method), linewidth = 0.9) +
  geom_line(data = truth, aes(t, mu), colour = "black", linetype = "dashed",
            inherit.aes = FALSE) +
  labs(title = "Population curve",
       y = "y") +
  theme_minimal()

# subject curves. snmmTMB and snmmAGQ only

subj_list <- vector("list", 0)

for (k in seq_along(subjects_to_plot)) {
  s <- subjects_to_plot[k]

  b <- b_mat[s, ]
  mu_true_s <- alpha + b["b1"] +
    2 * exp(b["b2"]) * sin(2 * pi * (t_grid - exp(b["b3"]) / (1 + exp(b["b3"]))))

  # snmmTMB
  idx_s  <- grep("mu_sel", names(rep$value))[(ngrid * (k - 1) + 1):(ngrid * k)]
  mu_s   <- rep$value[idx_s]
  Cov_s  <- rep$cov[idx_s, idx_s]
  Cov_s  <- (Cov_s + t(Cov_s)) / 2
  se_s   <- sqrt(diag(Cov_s))
  sims_s <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_s)), Sigma = Cov_s)
  crit_s <- quantile(apply(abs(sweep(sims_s, 1, se_s, FUN = "/")), 1, max), 0.95)

  # snmmAGQ: X_i carries the fixed-effect uncertainty, Z_i the Jacobian of the
  # subject curve in b_i, so the two sources add
  b_hat_i <- modes[s, ]
  X_i <- cbind(1, exp(b_hat_i[2]) *
                 sweep(splineDesign(knotseq, t_grid - expit(b_hat_i[3]), 4,
                                    outer.ok = TRUE),
                       2, col_means, FUN = "-"))
  Z_i <- numDeriv::jacobian(function(b) subj_curve_fun(b, alpha_agq, a, knotseq, col_means),
                            b_hat_i)
  Var_i  <- make_pd(X_i %*% Ginv %*% t(X_i) + Z_i %*% solve(hessians[[s]]) %*% t(Z_i))
  se_i   <- sqrt(diag(Var_i))
  sims_i <- mvrnorm(Nsims, mu = rep(0, ncol(Var_i)), Sigma = Var_i)
  crit_i <- quantile(apply(abs(sweep(sims_i, 2, se_i, FUN = "/")), 1, max), 0.95)
  mu_i   <- subj_curve_fun(b_hat_i, alpha_agq, a, knotseq, col_means)

  subj_list[[k]] <- bind_rows(
    tibble(subject = s, method = "snmmTMB", t = t_grid,
           mu = mu_s, truth = as.numeric(mu_true_s),
           lwr_s = mu_s - crit_s * se_s, upr_s = mu_s + crit_s * se_s),
    tibble(subject = s, method = "snmmAGQ", t = t_grid,
           mu = mu_i, truth = as.numeric(mu_true_s),
           lwr_s = mu_i - crit_i * se_i, upr_s = mu_i + crit_i * se_i)
  )
}

subj <- bind_rows(subj_list)
subj$method <- factor(subj$method, levels = c("snmmTMB", "snmmAGQ"))

ggplot(subj, aes(t)) +
  geom_ribbon(aes(ymin = lwr_s, ymax = upr_s), fill = "grey75", alpha = 0.7) +
  geom_point(data = dat[dat$subject %in% subjects_to_plot, ], aes(t, y),
             colour = "grey40", size = 1.2) +
  geom_line(aes(y = mu), linewidth = 0.9) +
  geom_line(aes(y = truth), linetype = "dashed") +
  facet_grid(method ~ paste("Subject", subject)) +
  labs(title = "Subject-specific curves",
       y = "y") +
  theme_minimal()
