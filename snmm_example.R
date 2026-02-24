# ------------------------------------------------------------------------------
# snmmTMB_example.R
#
# This script illustrates the use of snmmTMB on simulated longitudinal data of
# a sine curve with subject-specific amplitude and phase.
#
# 1. Generate the data
# 2. Construct the knot vector and penalty matrices
# 3. Fit the SNMM via TMB
# 4. Visualize the results
# ------------------------------------------------------------------------------

set.seed(1)

library(MASS)
library(dplyr)
library(ggplot2)
library(TMB)




# 1. Generate the data

n <- 10        # number of subjects
m <- 10       # number of measurements per subject
N <- n*m
sigma <- 0.4  # observation noise sd
D_diag <- c(1, 0.25, 0.16)  # diagonal multiplier for random-effect covariance

alpha <- 5 # intercept
Sigma_b <- sigma^2 * diag(D_diag)
b_mat <- mvrnorm(n = n, mu = rep(0,3), Sigma = Sigma_b)
colnames(b_mat) <- c("b1","b2","b3")

t_grid <- seq(0, 1, length.out = 50) # grid to evaluate the function on




dat_list <- vector("list", n)
for(i in 1:n){
  b1 <- b_mat[i,"b1"]
  b2 <- b_mat[i,"b2"]
  b3 <- b_mat[i,"b3"]
  
  t <- seq(0, 1, length.out = m)
  
  amp <- 2 * exp(b2) 
  phase <- exp(b3) / (1 + exp(b3))
  mu_ij <- alpha + b1 + amp * sin(2 * pi * (t - phase)) # subject specific curve
  y_ij  <- mu_ij + rnorm(m, mean = 0, sd = sigma)
  dat_list[[i]] <- data.frame(
    subject = rep(i,m),
    j = 1:m,
    t = t,
    y = y_ij,
    mu = mu_ij,
    b1 = rep(b1,m),
    b2 = rep(b2,m),
    b3 = rep(b3,m)
  )
}


dat <- bind_rows(dat_list)

mu_pop <- alpha + 2 * sin(2 * pi * (t_grid - 0.5)) # true population curve
pop_df <- data.frame(t = t_grid, mu = mu_pop)

ggplot(dat, aes(x = t, y = y, colour = factor(subject))) +
  geom_point(alpha = 0.45, size = 1.1, show.legend = FALSE) +
  geom_line(aes(y = mu, group = subject), data = dat[,c('subject','t','mu')] %>% distinct(),
            colour = "grey50", linewidth = 0.9, alpha = 0.3) +
  geom_line(data = pop_df, aes(x = t, y = mu), colour = "black", linewidth = 1.2) +
  labs(title = "Simulated data",
       subtitle = "black = population mean; grey = subject-specific curves",
       x = "t = j/n", y = "y_ij") +
  theme_minimal()












# 2. Construct the knot vector and penalty matrices

degree <- 3
K <- 15 # number of basis functions
n_internal <- K - (degree + 1)

# knot range
knot_lowlim <- -1
knot_uplim <- 1
internal_knots <- seq(knot_lowlim, knot_uplim, length.out = n_internal+2)

if (length(internal_knots) > 1) {
  dx_left  <- diff(internal_knots)[1]
  dx_right <- diff(internal_knots)[length(internal_knots)-1]
} else {
  dx_left <- dx_right <- 1
}

knot_vec <- c(internal_knots[1] - dx_left * (degree:1), internal_knots, internal_knots[length(internal_knots)] + dx_right * (1:degree))

D <- diff(diag(K), differences=2)[, 1:(K-1)] # drop one column to enforce the sum-to-zero constraint on the smooth
P <- t(D) %*% D
P <- P * qr(P)$rank / sum(diag(P))

# decompose into null and active space of the penalty
ev <- eigen(P, symmetric = TRUE)
pos_idx <- which(ev$values > 1e-12)
zero_idx <- which(ev$values <= 1e-12)

Upos <- if(length(pos_idx) > 0) ev$vectors[, pos_idx, drop = FALSE] else matrix(0, nrow = K, ncol = 0)
dpos <- if(length(pos_idx) > 0) ev$values[pos_idx] else numeric(0)
U0 <- if(length(zero_idx) > 0) ev$vectors[, zero_idx, drop = FALSE] else matrix(0, nrow = K, ncol = 0)












# 3. Fit the SNMM via TMB

cpp_file <- "snmmTMB_likelihood_sinecurve_simulation.cpp"
# the first time the model is fitted, the C++ file needs to be compiled
TMB::compile(cpp_file) 
dyn.load(dynlib(sub("\\.cpp$", "", cpp_file)))


# data to give TMB
Data <- list(
  y = dat$y,
  t = dat$t,
  group = as.integer(dat$subject), # integer 1 to n
  nGroup = as.integer(n),
  subj_flag = as.integer(c(1,1,rep(0,n-2))), # boolean: which subject curves to get CI for
  knots = as.numeric(knot_vec),
  degree = as.integer(degree),
  K = as.integer(K),
  Upos = as.matrix(Upos),
  U0 = as.matrix(U0),
  dpos = as.numeric(dpos),
  spline_ci = as.integer(1), # boolean: get the CI for the population curve
  t_grid = as.numeric(t_grid) # time values to evaluate the curve on
)


# parameters initialization
# in this example we initialize at standard values, otherwise the procedure to find suitable starting points can be used
Params <- list(
  alpha = 0,
  b1 = rep(0, n),
  b2 = rep(0, n),
  b3 = rep(0, n),
  log_sd_b1 = log(0.5),
  log_sd_b2 = log(0.5),
  log_sd_b3 = log(0.5),
  log_sigma = log(0.5),
  vpos = rep(0, length(dpos)),
  gamma0 = rep(0, ncol(U0)),
  log_lambda = log(1)
)

# declare random parameters to integrate out
randoms <- c("b1", "b2", "b3", "vpos")

# make function to evaluate log marginal likelihood
obj <- MakeADFun(data = Data, parameters = Params, random = randoms, DLL = sub("\\.cpp$", "", cpp_file), silent = TRUE)

# optimize
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 1e4, iter.max = 1e4))
obj$par <- opt$par

# compute inference for reported quantities with delta method
rep <- sdreport(obj, getJointPrecision = TRUE)


















# 4. Visualize the results

# summary including all reported variables (fixed and random effects, evaluated curves)
print(summary(rep))





# plot true and predicted random effects

parList <- obj$env$parList()
b1_hat <- if("b1" %in% names(parList)) parList$b1 else rep(0,n)
b2_hat <- if("b2" %in% names(parList)) parList$b2 else rep(0,n)
b3_hat <- if("b3" %in% names(parList)) parList$b3 else rep(0,n)

df <- data.frame(
  true = c(b_mat[,'b1'], b_mat[,'b2'], b_mat[,'b3']),
  predicted = c(b1_hat, b2_hat, b3_hat),
  b_element = rep(c("b1", "b2", "b3"), each = nrow(b_mat))
)

# Plot
ggplot(df, aes(x = true, y = predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~b_element, scales = "free") +
  labs(x = "True b", y = "Predicted b", title = "True vs Predicted b") +
  theme_bw()







# plot true and predicted population curve

idx <- grep("h_grid", names(rep$value))
mu_hat <- rep$value[idx]  
Cov_mu <- rep$cov[idx, idx]

Cov_mu <- (Cov_mu + t(Cov_mu)) / 2  

# obtain simultaneous critical value through simulation
Nsims <- 10000
mu_sims <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu)), Sigma = Cov_mu)

se_point <- sqrt(diag(Cov_mu))
absDev <- abs(sweep(mu_sims, 1, se_point, FUN = "/"))


max_dev <- apply(abs(absDev), 1, max)
crit <- quantile(max_dev, 0.95)

df <- data.frame(x = t_grid, y = mu_hat,
                 lwrS = mu_hat - crit*se_point,
                 uprS = mu_hat + crit*se_point,
                 lwrP = mu_hat - 1.96*se_point,
                 uprP = mu_hat + 1.96*se_point,
                 truth = alpha + 2 * sin(2 * pi * (t_grid-0.5)))


ggplot(df) +
  geom_ribbon(aes(x = x, ymin = lwrS, ymax = uprS, fill = "Simultaneous"), alpha = 0.2) +
  geom_ribbon(aes(x = x, ymin = lwrP, ymax = uprP, fill = "Pointwise"), alpha = 0.2) +
  geom_path(aes(x = x, y = y, colour = "Estimated"), lwd = 2) +
  geom_line(aes(x = x, y = truth, colour = "Truth"), linetype = "dashed", size = 1.1) +
  scale_fill_manual(name = "Confidence Intervals", 
                    values = c('Simultaneous' = 'red', 'Pointwise' = 'blue')) +
  scale_color_manual(name = "Model Fit", 
                     values = c("Estimated" = "blue", "Truth" = "black")) +
  labs(y = "y",
       x = "u",
       title = "Point-wise & Simultaneous 95% confidence intervals")









# plot estimated curves for individual subjects

for(subject_index in 1:2){
  
  b1 <- b_mat[subject_index,"b1"]
  b2 <- b_mat[subject_index,"b2"]
  b3 <- b_mat[subject_index,"b3"]
  
  mu_true_subj <- alpha + b1 + 2 * exp(b2) * sin(2 * pi * (t_grid - exp(b3) / (1 + exp(b3))))

  
  idx <- grep("mu_sel", names(rep$value))[(50*(subject_index-1)+1):(50*subject_index)]
  mu_hat_subj <- rep$value[idx]  
  Cov_mu_subj <- rep$cov[idx, idx]
  
  
  Cov_mu_subj <- (Cov_mu_subj + t(Cov_mu_subj)) / 2  
  
  
  Nsims <- 10000
  mu_sims_subj <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu_subj)), Sigma = Cov_mu_subj)
  
  
  se_point_subj <- sqrt(diag(Cov_mu_subj)) 
  absDev_subj <- abs(sweep(mu_sims_subj, 1, se_point_subj, FUN = "/"))
  
  max_dev_subj <- apply(abs(absDev_subj), 1, max)
  crit_subj <- quantile(max_dev_subj, 0.95)  
  
  
  # Plotting df
  df_sub <- data.frame(
    x_true = t_grid,
    y_true = mu_true_subj,
    x = t_grid,
    mu_hat = mu_hat_subj,
    se = se_point_subj,
    lwr_point = mu_hat_subj - 1.96 * se_point_subj,
    upr_point = mu_hat_subj + 1.96 * se_point_subj,
    lwr_sim = mu_hat_subj - as.numeric(crit_subj) * se_point_subj,
    upr_sim = mu_hat_subj + as.numeric(crit_subj) * se_point_subj,
    subject = factor(i)
  )
  
  p <- ggplot(data = df_sub) + 
    geom_ribbon(aes(x = x, ymin = lwr_sim, ymax = upr_sim, fill = 'Simultaneous'), 
                alpha = 0.2) +
    geom_ribbon(aes(x = x, ymin = lwr_point, ymax = upr_point, fill = 'Pointwise'), 
                alpha = 0.2) +
    geom_line(aes(x = x, y = mu_hat, color = 'Estimated'), size = 1.1) +
    geom_line(aes(x = x_true, y = y_true, color = 'Truth'), lty = 'longdash', size = 1.1) +
    geom_point(data = dat[dat$subject == subject_index,], 
               aes(x = t, y = y, color = 'Observed'), size = 3) +
    scale_fill_manual(name = "Confidence Bands", 
                      values = c('Simultaneous' = 'red', 'Pointwise' = 'blue')) +
    scale_color_manual(name = "Model Fit", 
                       values = c('Estimated' = 'blue', 'Truth' = 'black', 'Observed' = 'green')) +
    labs(x = "u", y = "y",
         title = paste0("Subject ", subject_index)) +
    theme_minimal()
  
  print(p)
  

}

