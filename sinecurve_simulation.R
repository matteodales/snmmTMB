set.seed(123)

library(MASS)
library(dplyr)
library(ggplot2)
library(TMB)
library(nlme)
library(splines)




# construct knot vector and penalty matrix

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



# functions to compute coverage

inCI <- function(x, upr, lwr) {
  all(x >= lwr & x <= upr)
}

meanCI <- function(x, upr, lwr) {
  mean(x >= lwr & x <= upr)
}







# Define parameter grid for simulation settings

iterations <- 300

n_vals <- c(10, 20)
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

# Store results in a nested list
results <- vector("list", nrow(param_grid))


for (setting_idx in 1:16) {
  
  
  n <- param_grid$n[setting_idx]
  m <- param_grid$m[setting_idx]
  N <- n*m
  sigma <- param_grid$sigma[setting_idx]
  D_diag <- D_diag_vals[[param_grid$D_idx[setting_idx]]]
  subj_flag_vec <- c(1,1,rep(0,n-2))
  
  # Pre-allocate result storage for this setting
  setting_results <- list(
    pointwise_coverage = rep(NA,iterations),
    simultaneous_coverage = rep(NA,iterations),
    subject_pointwise_coverage = rep(NA,iterations),
    subject_simultaneous_coverage = rep(NA,iterations),
    pointwise_CI_length = rep(NA,iterations),
    simultaneous_CI_length = rep(NA,iterations),
    subject_pointwise_CI_length = rep(NA,iterations),
    subject_simultaneous_CI_length = rep(NA,iterations),
    b1_mse = rep(NA,iterations),
    b2_mse = rep(NA,iterations),
    b3_mse = rep(NA,iterations),
    b1_coverage = rep(NA,iterations),
    b2_coverage = rep(NA,iterations),
    b3_coverage = rep(NA,iterations),
    bias_vec = vector("list", iterations),
    b1 = vector("list", iterations),
    b2 = vector("list", iterations),
    b3 = vector("list", iterations),
    time = rep(NA,iterations)
  )
  
  

  
  # Simulation loop
  for (iter in 1:iterations) {
    
  print(iter)
    
    
  # Generate data
  Sigma_b <- sigma^2 * diag(D_diag)
  b_mat <- mvrnorm(n = n, mu = rep(0,3), Sigma = Sigma_b)
  colnames(b_mat) <- c("b1","b2","b3")
  
  t_grid <- seq(0, 1, length.out = 50) # grid to evaluate the function on
  alpha = 5
  
  
  
  dat_list <- vector("list", n)
  for(i in 1:n){
    b1 <- b_mat[i,"b1"]
    b2 <- b_mat[i,"b2"]
    b3 <- b_mat[i,"b3"]
    
    t <- seq(0, 1, length.out = m)
    
    amp <- 2 * exp(b2) 
    phase <- exp(b3) / (1 + exp(b3))
    mu_ij <- alpha + b1 + amp * sin(2 * pi * (t - phase))
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
  
  
  
  
  
  
  
  
  ### get starting points for spline coefficients, fixed and random effects
  
  # fit a simple GAM with the same knot vector and penalty matrix as the final model
  
  
  Params <- list(
    beta1 = as.numeric(0), #intercept
    log_sigma = log(1),
    vpos = rep(0, length(dpos)),
    gamma0 = rep(0, ncol(U0)),
    log_lambda = log(1)
  )
  
  

  cpp_file <- "starting_points.cpp"
  # TMB::compile(cpp_file)
  dyn.load(dynlib(sub("\\.cpp$", "", cpp_file)))
  
  Data <- list(
    y = dat$y,
    x = dat$t-0.5,
    knots = as.numeric(knot_vec),
    degree = as.integer(degree),
    K = as.integer(K),
    Upos = as.matrix(Upos),
    U0 = as.matrix(U0),
    dpos = as.numeric(dpos),
    spline_ci = as.integer(1),
    x_grid = t_grid-0.5
  )
  
  obj <- MakeADFun(data = Data, parameters = Params, random = c("vpos"), DLL = sub("\\.cpp$", "", cpp_file), silent = TRUE)
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  obj$par <- opt$par
  rep <- sdreport(obj, getJointPrecision = TRUE)
  parList <- obj$env$parList()
  
  idx <- grep("c", names(rep$value))
  c <- rep$value[idx]
  m_var <- as.numeric(obj$report()$m)
  
  spline_start <-function(t){
    result <- sweep(as.matrix(splineDesign(knot_vec,t,4,outer.ok=T))[,1:(length(knot_vec)-5), drop = FALSE], 2, m_var, FUN="-")%*%c
    result
  }
  
  .env=.GlobalEnv
  assign("spline_start",spline_start,env=.env)
  
  ## run nlme model with spline_start
  
  nlmeobj <- tryCatch(nlme(y ~ b1+exp(b2)*spline_start(t-exp(b3)/(1 + exp(b3))),
                           fixed=list(b1~1),
                           random=pdDiag(b1+b2+b3~1),
                           data=dat[,1:4],
                           groups=~subject,
                           verbose=F,
                           start=mean(dat$y),
                           control=list(returnObject=T,tolerance=.01)),
                      error = function(e) {
                        message("Skipping iteration due to error: ", e$message)
                        return(NULL)
                        }
                      )
  
  if(is.null(nlmeobj)){next}

  
  # fitting model with TMB
  
  cpp_file <- "snmmTMB_likelihood_sinecurve_simulation.cpp"
  # TMB::compile(cpp_file)
  dyn.load(dynlib(sub("\\.cpp$", "", cpp_file)))
  
  # Data for TMB
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
  
  # Parameters initialization
  Params <- list(
    alpha = nlmeobj$coefficients$fixed["b1"],
    b1 = as.numeric(nlmeobj$coefficients$random$subject[,"b1"]),
    b2 = as.numeric(nlmeobj$coefficients$random$subject[,"b2"]),
    b3 = as.numeric(nlmeobj$coefficients$random$subject[,"b3"]),
    log_sd_b1 = log(as.numeric(VarCorr(nlmeobj)["b1","StdDev"])),
    log_sd_b2 = log(as.numeric(VarCorr(nlmeobj)["b2","StdDev"])),
    log_sd_b3 = log(as.numeric(VarCorr(nlmeobj)["b3","StdDev"])),
    log_sigma = log(nlmeobj$sigma),
    vpos = as.numeric(parList$vpos),
    gamma0 = as.numeric(parList$gamma0),
    log_lambda = as.numeric(obj$par["log_lambda"])
  )
  
  
  randoms <- c("b1", "b2", "b3", "vpos")
  
  start_time <- Sys.time()
  
  obj <- MakeADFun(data = Data, parameters = Params, random = randoms, DLL = sub("\\.cpp$", "", cpp_file), silent = TRUE)
  
  # Optimization
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(eval.max = 1e4, iter.max = 1e4))
  
  obj$par <- opt$par
  rep <- sdreport(obj, getJointPrecision = TRUE)
  
  end_time  <- Sys.time()
  time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))
  

  
  
  
  # fixed parameters
  true_par <- c(alpha, log(sqrt(diag(Sigma_b))), log(sigma))
  bias_vec <- true_par - opt$par[1:5]
  
  
  parList <- obj$env$parList()
  # b1_hat, b2_hat, b3_hat
  b1_hat <- if("b1" %in% names(parList)) parList$b1 else rep(0,n)
  b2_hat <- if("b2" %in% names(parList)) parList$b2 else rep(0,n)
  b3_hat <- if("b3" %in% names(parList)) parList$b3 else rep(0,n)
  
  se <- summary(rep)
  se_b1 <- se[rownames(se) == "b1","Std. Error"]
  se_b2 <- se[rownames(se) == "b2","Std. Error"]
  se_b3 <- se[rownames(se) == "b3","Std. Error"]
  
  b1_mse <- sqrt(sum((b_mat[,'b1']-b1_hat)**2)/length(b1_hat))
  b2_mse <- sqrt(sum((b_mat[,'b2']-b2_hat)**2)/length(b2_hat))
  b3_mse <- sqrt(sum((b_mat[,'b3']-b3_hat)**2)/length(b3_hat))
  
  b1_bias <- sum((b_mat[,'b1']-b1_hat))/length(b1_hat)
  b2_bias <- sum((b_mat[,'b2']-b2_hat))/length(b2_hat)
  b3_bias <- sum((b_mat[,'b3']-b3_hat))/length(b3_hat)
  
  cover_b1 <- mean(b_mat[,"b1"] >= (b1_hat - 1.96*se_b1) &
                     b_mat[,"b1"] <= (b1_hat + 1.96*se_b1))
  cover_b2 <- mean(b_mat[,"b2"] >= (b2_hat - 1.96*se_b2) &
                     b_mat[,"b2"] <= (b2_hat + 1.96*se_b2))
  cover_b3 <- mean(b_mat[,"b3"] >= (b3_hat - 1.96*se_b3) &
                     b_mat[,"b3"] <= (b3_hat + 1.96*se_b3))
  
  bias_vec <- c(bias_vec, b1_bias, b2_bias, b3_bias)
  
  
  
  
  
  
  
  
  # population curve coverage
  
  mu_true <- alpha + 2 * sin(2 * pi * (t_grid - 0.5)) # true population curve
  
  idx <- grep("h_grid", names(rep$value))
  mu_hat <- rep$value[idx]  
  Cov_mu <- rep$cov[idx, idx]
  
  # symmetric
  Cov_mu <- (Cov_mu + t(Cov_mu)) / 2  

  # simulation for simultaneous coverage critical value
  Nsims <- 10000
  mu_sims <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu)), Sigma = Cov_mu)
  
  se_point <- sqrt(diag(Cov_mu))  # pointwise SE
  absDev <- abs(sweep(mu_sims, 1, se_point, FUN = "/"))
  
  # take maximum absolute deviation over the grid for each simulation
  max_dev <- apply(abs(absDev), 1, max)
  crit <- quantile(max_dev, 0.95) # 95% simultaneous critical value
  
  # save
  pointwise_coverage <- meanCI(mu_true, upr = mu_hat + 1.96*se_point, lwr = mu_hat - 1.96*se_point)
  simultaneous_coverage <-  inCI(df$truth, upr = mu_hat + crit*se_point, lwr = mu_hat - crit*se_point)
  pointwise_CI_length <- mean(2*1.96*se_point)
  simultaneous_CI_length <- mean(2*crit*se_point)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # subject curve coverage
  
  sim_cov_subj <- c(0,0)
  point_cov_subj <- c(0,0)
  
  CI_length_point_subj <- c(0,0)
  CI_length_sim_subj <- c(0,0)
  
  
  for(subject_index in 1:2){
    
    # true individual curve
    
    b1 <- b_mat[subject_index,"b1"]
    b2 <- b_mat[subject_index,"b2"]
    b3 <- b_mat[subject_index,"b3"]
    
    amp <- 2 * exp(b2) 
    phase <- exp(b3) / (1 + exp(b3))   # logistic(b3)
    mu_true_subj <- b1 + amp * sin(2 * pi * (t_grid - phase))


    # obtain estimated curve and CI
    
    idx <- grep("mu_sel", names(rep$value))[(50*(subject_index-1)+1):(50*subject_index)]
    mu_hat_subj <- rep$value[idx]  
    Cov_mu_subj <- rep$cov[idx, idx]
    
    
    Cov_mu_subj <- (Cov_mu_subj + t(Cov_mu_subj)) / 2  
    
    
    Nsims <- 10000
    mu_sims_subj <- mvrnorm(Nsims, mu = rep(0, ncol(Cov_mu_subj)), Sigma = Cov_mu_subj)
    
    
    se_point_subj <- sqrt(diag(Cov_mu_subj))           # pointwise SE
    absDev_subj <- abs(sweep(mu_sims_subj, 1, se_point_subj, FUN = "/"))
    
    # Take maximum absolute deviation over the grid for each simulation
    max_dev_subj <- apply(abs(absDev_subj), 1, max)
    crit_subj <- quantile(max_dev_subj, 0.95)         # 95% simultaneous critical value

    sim_cov_subj[subject_index] <- inCI(mu_true_subj, upr = mu_hat_subj + as.numeric(crit_subj) * se_point_subj, lwr = mu_hat_subj - as.numeric(crit_subj) * se_point_subj)
    point_cov_subj[subject_index] <- meanCI(mu_true_subj, upr = mu_hat_subj + 1.96 * se_point_subj, lwr = mu_hat_subj - 1.96 * se_point_subj)
    CI_length_sim_subj[subject_index] <- mean(2*as.numeric(crit_subj) * se_point_subj)
    CI_length_point_subj[subject_index] <- mean(2 * 1.96 * se_point_subj)
    
  }
  
  
  setting_results$pointwise_coverage[iter] <- pointwise_coverage
  setting_results$simultaneous_coverage[iter] <- simultaneous_coverage
  
  
  setting_results$subject_pointwise_coverage[iter] <- mean(point_cov_subj)
  setting_results$subject_simultaneous_coverage[iter] <- mean(sim_cov_subj)
  
  
  setting_results$pointwise_CI_length[iter] <- pointwise_CI_length
  setting_results$simultaneous_CI_length[iter] <- simultaneous_CI_length
  
  setting_results$subject_pointwise_CI_length[iter] <- mean(CI_length_point_subj)
  setting_results$subject_simultaneous_CI_length[iter] <- mean(CI_length_sim_subj)
  
  
  setting_results$b1_mse[iter] <- b1_mse
  setting_results$b2_mse[iter] <- b2_mse
  setting_results$b3_mse[iter] <- b3_mse
  setting_results$b1_coverage[iter] <- cover_b1
  setting_results$b2_coverage[iter] <- cover_b2
  setting_results$b3_coverage[iter] <- cover_b3
  
  setting_results$bias_vec[[iter]] <- bias_vec
  setting_results$b1[[iter]] <- b1_hat
  setting_results$b2[[iter]] <- b2_hat
  setting_results$b3[[iter]] <- b3_hat
  
  setting_results$time[[iter]] <- time_diff
  
  }
  
  
  
  # Save this setting's results
  results[[setting_idx]] <- list(
    params = param_grid[setting_idx, ],
    metrics = setting_results
  )
  
  saveRDS(results, file = "check_sinesimulation_coverage_results_with_time.RDS")
}

saveRDS(results, file = "check_sinesimulation_coverage_results_with_time.RDS")
