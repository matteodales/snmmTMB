# ------------------------------------------------------------------------------
# snmmAGQ_functions_bellcurve.R
#
# Estimation routines for snmmAGQ, the adaptive Gauss-Hermite quadrature
# procedure of Elmi et al. (2011), adapted from the supplementary code published
# with that paper, for the BELL-CURVE simulation of manuscript Section 3:
#
#   y_ij = alpha + b1_i + f(t_ij - b2_i) + eps_ij
#   b_i = (b1_i, b2_i) ~ N(0, D), D unstructured with correlation rho
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(splines)   # splineDesign()
  library(nlme)      # nlme() inside startingvalues()
  library(statmod)   # gauss.quad() for the quadrature nodes
  library(numDeriv)  # hessian(), jacobian()
})

# knot placing functions from Elmi
knots.equispaced <- function(nk, xl, xr){
  #purpose: to compute equispaced knot locations for a cubic B-spline
  #nk=the number of knots
  #xl=left boundary of the covariate domain
  #xr=right boundary of the covariate domain
  dx <- (xr-xl)/(nk+1)
  inner.knots <- seq(xl,xr,dx)
  outer.knots <- c(xl-3*dx,xl-2*dx,xl-dx,xr+dx,xr+2*dx,xr+3*dx)
  sort(c(inner.knots,outer.knots))
}
knots.equiquantile <- function(nk, t, xl, xr){
  #purpose: to compute equiquantile knot locations for a cubic B-spline
  #nk=the number of knots
  #t=data values the interior knots are placed at the quantiles of
  #xl=left boundary of the covariate domain
  #xr=right boundary of the covariate domain
  dx <- 1/(nk+1)
  probs <- seq(0,1,dx)
  quantiles <- quantile(t, probs)
  inner.knots <- quantiles
  outer.knots <- c(xl-3*dx,xl-2*dx,xl-dx,xr+dx,xr+2*dx,xr+3*dx)
  sort(c(inner.knots,outer.knots))
}



agq_cache <- function(y_list, t_list, knotseq, col_means, modes, hessians){
  # Per-subject cache of everything needed for the AGQ likelihood
  #
  # Within one iteration of bspline.snm() below, the posterior modes, the
  # posterior Hessians, the knot sequence and the centering constants are all
  # held fixed while a, alpha and theta are optimized
  #
  # Everything else is constant across neglogLik() and is built once here:
  #   - SVD of the posterior Hessian, the resulting quadrature scaling
  #     matrix and its log-determinant,
  #   - adaptive quadrature nodes Bcand,
  #   - centered B-spline design matrix evaluated at those nodes
  #   - random-effect columns.
  #
  # y_list, t_list: per-subject list of response vectors / observation times
  # knotseq, col_means: spline knot sequence / basis column means (centering)
  # modes: matrix of posterior modes, one row per subject
  # hessians: list of posterior Hessians, one per subject
  
  lapply(seq_along(y_list), function(i){
    y <- y_list[[i]]; t <- t_list[[i]]

    sv <- svd(hessians[[i]])
    scale_mat <- sv$v %*% diag(1/sqrt(sv$d)) %*% t(sv$u)

    # the original script has a nested loop for the two random effects, here it's vectorized for speed
    Bcand <- sweep(sqrt(2) * (AGQ_NODES %*% t(scale_mat)), 2, modes[i,], FUN = "+")

    nq <- nrow(Bcand); mrows <- length(y)

    u_all <- rep(t, times = nq) - rep(Bcand[,2], each = mrows)  # covariate shifted by b2
    Bspl  <- sweep(splineDesign(knotseq, u_all, 4, outer.ok = TRUE), 2, col_means, FUN = "-") # (nq*mrows) x K, centered

    list(nq = nq, mrows = mrows, Bcand = Bcand, Bspl = Bspl,
         y_rep = rep(y, times = nq),
         b1_rep = rep(Bcand[,1], each = mrows),
         logdet_scale = log(det(scale_mat)))
  })
}


JointLogDensity <- function(cache_i, a, alpha, theta){
  # Log joint density of (y_i, b_i) for one subject, as in Eq. (2.9) of Elmi et al. (2011)
  #
  # cache_i: this subject's entry from agq_cache()
  # a: spline coefficients; alpha: fixed intercept; theta: (sd_b1,sd_b2,sigma,rho)

  s1 <- theta[1]; s2 <- theta[2]; rr <- theta[4]
  omr2 <- 1 - rr^2
  detD <- s1^2 * s2^2 * omr2
  Dinv <- matrix(c(1/(s1^2), -rr/(s1*s2),
                   -rr/(s1*s2), 1/(s2^2)), nrow = 2) / omr2
  Bmat <- cache_i$Bcand
  nq <- cache_i$nq; mrows <- cache_i$mrows

  h_all <- as.numeric(cache_i$Bspl %*% a)

  mu    <- alpha + cache_i$b1_rep + h_all                       # subject prediction
  resid <- cache_i$y_rep - mu
  rss   <- colSums(matrix(resid^2, nrow = mrows, ncol = nq))    # RSS per quadrature node

  quad_b <- rowSums((Bmat %*% Dinv) * Bmat)

  log_dens_y <- -0.5*mrows*log(2*pi*theta[3]^2) - 0.5*rss/theta[3]^2
  log_dens_b <- -1.0*log(2*pi) - 0.5*log(detD) - 0.5*quad_b     # -(d/2)*log(2*pi), d=2
  log_dens_y + log_dens_b
}


postmode <- function(b, y, t, a, alpha, knotseq, col_means, theta){
  # negative log posterior for one subject, optimized by nlm in initialest.re() to find the posterior mode of b_i (Eq. 2.10)
  #
  # b: candidate random-effects vector (b1,b2)
  # y, t: subject's response vector / covariate values
  # a, alpha: current spline coefficients / fixed intercept
  # knotseq, col_means: spline knot sequence / basis column means (centering)
  # theta: variance components (sd_b1,sd_b2,sigma,rho)

  # same closed-form correlated Dinv as in JointLogDensity()
  s1 <- theta[1]; s2 <- theta[2]; rr <- theta[4]
  omr2 <- 1 - rr^2
  Dinv <- matrix(c(1/(s1^2), -rr/(s1*s2),
                   -rr/(s1*s2), 1/(s2^2)), nrow = 2) / omr2
  u <- t - b[2]
  Bspl <- sweep(splineDesign(knotseq, u, 4, outer.ok=TRUE), 2, col_means, FUN = "-")
  resid <- y - alpha - b[1] - as.numeric(Bspl %*% a) #model equation
  sum(resid^2)/theta[3]^2 + as.numeric(t(b) %*% Dinv %*% b)
}


initialest.re <- function(y, t, a, alpha, knotseq, col_means, theta, b.start){
  # finds the posterior mode of b_i and the Hessian of the negative log posterior at the mode (Eq. 2.11)
  #
  # y, t: subject's response vector / covariate values
  # a, alpha, knotseq, col_means, theta: current parameter estimates, passed through to postmode()
  # b.start: starting value for b_i in the nlm optimization

  q2 <- nlm(postmode, b.start, y=y, t=t, a=a, alpha=alpha, knotseq=knotseq, col_means=col_means, theta=theta, hessian=TRUE)
  mode <- q2$estimate
  hess <- q2$hessian
  list(mode=mode, hessian=hess)
}


MarginalLogLik <- function(cache_i, a, alpha, theta){
  # approximate a subject's marginal log-likelihood via adaptive Gauss-Hermite quadrature (Eq. 2.12)
  #
  # cache_i: this subject's entry from agq_cache()
  # a, alpha, theta: current parameter estimates

  log_dens <- JointLogDensity(cache_i, a, alpha, theta)

  log_terms <- AGQ_LOGW + log_dens + AGQ_ZSQ
  mx <- max(log_terms)
  logsum <- mx + log(sum(exp(log_terms - mx)))

  log(2) + cache_i$logdet_scale + logsum # (d/2)*log(2), d=2
}

neglogLik <- function(cache, a, alpha, theta){
  # sum over marginal log-likelihood for all subjects
  #
  # cache: list of per-subject caches from agq_cache(), one element per subject
  # a, alpha, theta: current parameter estimates

  ll <- 0
  for(i in seq_along(cache)){
    ll <- ll + MarginalLogLik(cache[[i]], a, alpha, theta)
  }
  -ll
}


# wrapper functions for optimizing negloglik w.r.t. the different parameters
# the original code just calls nlm to negloglik with the different parameters as first argument
nll_a <- function(a, cache, alpha, theta)
  # a: spline coefficients
  neglogLik(cache, a, alpha, theta)

nll_alpha <- function(alpha, cache, a, theta)
  # alpha: fixed intercept
  neglogLik(cache, a, alpha, theta)

# theta = (sd_b1, sd_b2, sigma, rho) cannot be optimized on the log scale as a
# whole since rho lives in (-1, 1), so the optimizer works on the unconstrained
#     par = (log sd_b1, log sd_b2, log sigma, atanh(rho))
# and these two helpers map back and forth (same parameterisation as the TMB likelihood)

# unconstrained optimizer scale -> natural (sd_b1, sd_b2, sigma, rho)
theta_from_par <- function(par) c(exp(par[1:3]), tanh(par[4]))

# natural -> unconstrained, pulling rho off the boundary so atanh stays finite
par_from_theta <- function(theta) {
  r <- theta[4]
  if (!is.finite(r)) r <- 0
  c(log(theta[1:3]), atanh(max(min(r, 0.99), -0.99)))
}

nll_logtheta <- function(par, cache, a, alpha)
  # par: unconstrained variance components + atanh(rho)
  neglogLik(cache, a, alpha, theta_from_par(par))



startingvalues <- function(dat, knots){
  # starting values for a, alpha, knotseq, col_means, theta, and per-subject modes,
  # via a population B-spline shape fit + an approximating nlme() fit
  #
  # dat: tibble with this dataset's y, t, subject columns
  # knots: list(equispaced TRUE/FALSE, number of interior knots)

  y <- dat$y; t <- dat$t; g <- dat$subject

  # the reference horizontal shift is b2 = 0, so the shape fit is built directly on t
  if(knots[[1]]){
    knotseq <- knots.equispaced(knots[[2]], min(t), max(t))
  }else{
    knotseq <- knots.equiquantile(knots[[2]], t, min(t), max(t))
  }

  alpha <- mean(y)
  col_means <- colMeans(splineDesign(knotseq, t, 4, outer.ok = TRUE))

  B <- sweep(splineDesign(knotseq, t, 4, outer.ok=TRUE), 2, col_means, FUN = "-")
  a <- as.numeric(MASS::ginv(t(B) %*% B) %*% t(B) %*% (y - alpha)) # initial guess of spline shape

  f <- function(u) as.numeric(sweep(splineDesign(knotseq, u, 4, outer.ok=TRUE), 2, col_means, FUN = "-") %*% a)

  # `f` has to be made visible via .GlobalEnv but it's removed on exit from the function
  assign("f", f, envir = .GlobalEnv)
  on.exit(if (exists("f", envir = .GlobalEnv, inherits = FALSE)) rm("f", envir = .GlobalEnv), add = TRUE)


  nlmeobj <- tryCatch(
    nlme(y ~ b1 + f(t - b2),
         fixed  = list(b1~1),
         random = pdSymm(b1+b2~1), # pdSymm, not pdDiag: the starting values must estimate rho too
         data = data.frame(y,t,g),
         groups = ~g,
         start = alpha,
         verbose = FALSE,
         control = list(returnObject=TRUE, tolerance=.01, maxIter=25)),
    error = function(e){
      message("nlme starting-value search failed (", conditionMessage(e), "); falling back to simple starting values.")
      NULL
    })

  if(!is.null(nlmeobj)){
    alpha <- as.numeric(nlme::fixef(nlmeobj)["b1"])
    vcmat <- VarCorr(nlmeobj)
    # getVarCov() is not implemented for nlme objects, so the starting correlation is read off the pdSymm structure directly
    rho_start <- tryCatch({
      vc <- nlme::pdMatrix(nlmeobj$modelStruct$reStruct[[1]]) * nlmeobj$sigma^2
      r <- vc[1, 2] / sqrt(vc[1, 1] * vc[2, 2])
      if (is.finite(r)) max(min(r, 0.95), -0.95) else 0
    }, error = function(e) 0)
    theta <- c(as.numeric(vcmat[c("b1","b2","Residual"), "StdDev"]), rho_start)
    re <- ranef(nlmeobj)
    re <- re[order(as.numeric(rownames(re))), c("b1","b2")]
    modes <- as.matrix(re)
  }else{
    subj_mean_resid <- as.numeric(tapply(y - alpha, g, mean))
    b1_start <- subj_mean_resid - mean(subj_mean_resid)
    theta <- c(max(sd(b1_start), 0.2), 0.5, max(sd(y - alpha), 0.2), 0)
    modes <- cbind(b1 = b1_start, b2 = rep(0, length(unique(g))))
  }

  list(a=a, alpha=alpha, knotseq=knotseq, col_means=col_means, theta=theta, modes=modes)
}



bspline.snm <- function(dat, y_list, t_list, start, knots, control){
  # main algorithm from Elmi, corresponding to Equations (2.10)-(2.15) in the paper
  #
  # fits the model by alternating per-subject posterior modes/Hessians via
  # initialest.re() and nlm updates of a, then alpha, then theta on neglogLik()
  # until the log-likelihood converges (or max.iter is hit)
  #
  # dat: full dataset (tibble); y_list, t_list: per-subject response/covariate lists
  # start: starting values list from startingvalues() (a, alpha, theta, modes)
  # knots: knot-placement spec, same as startingvalues()
  # control: list(tol, max.iter, inner.iterlim)

  a <- start$a
  alpha <- start$alpha
  theta <- start$theta
  modes <- start$modes

  # the knot domain is recentered every iteration based on the shift estimates
  shift_hat0 <- modes[,2]
  u_all0 <- dat$t - shift_hat0[dat$subject]
  if(knots[[1]]){
    knotseq <- knots.equispaced(knots[[2]], min(u_all0), max(u_all0))
  }else{
    knotseq <- knots.equiquantile(knots[[2]], u_all0, min(u_all0), max(u_all0))
  }
  col_means <- colMeans(splineDesign(knotseq, dat$t, 4, outer.ok = TRUE))

  nsubj <- length(y_list)
  hessians <- vector("list", nsubj)
  loglik_prev <- -Inf
  iteration <- 0

  repeat{
    iteration <- iteration + 1

    for(i in 1:nsubj){
      pm <- initialest.re(y_list[[i]], t_list[[i]], a, alpha, knotseq, col_means, theta, modes[i,])
      modes[i,] <- pm$mode
      hessians[[i]] <- pm$hessian
    }

    # re-derive the knot domain (and column-centering constants)
    shift_hat <- modes[,2]
    u_all <- dat$t - shift_hat[dat$subject]
    if(knots[[1]]){
      knotseq <- knots.equispaced(knots[[2]], min(u_all), max(u_all))
    }else{
      knotseq <- knots.equiquantile(knots[[2]], u_all, min(u_all), max(u_all))
    }
    col_means <- colMeans(splineDesign(knotseq, dat$t, 4, outer.ok = TRUE))

    # modes, hessians, knotseq and col_means are now all fixed for the rest of
    # this iteration, so everything derived from them is built once here instead
    # of being recomputed inside every neglogLik() call
    cache <- agq_cache(y_list, t_list, knotseq, col_means, modes, hessians)

    # optimize the parameters
    a.obj <- nlm(nll_a, a, cache=cache,
                 alpha=alpha, theta=theta, iterlim=control$inner.iterlim)
    a <- a.obj$estimate

    alpha.obj <- nlm(nll_alpha, alpha, cache=cache,
                      a=a, theta=theta, iterlim=control$inner.iterlim)
    alpha <- alpha.obj$estimate

    theta.obj <- nlm(nll_logtheta, par_from_theta(theta), cache=cache,
                      a=a, alpha=alpha, iterlim=control$inner.iterlim)
    theta <- theta_from_par(theta.obj$estimate)

    logLik.new <- -theta.obj$minimum
    convcrit <- if(is.finite(loglik_prev)) (logLik.new-loglik_prev)/loglik_prev else Inf
    loglik_prev <- logLik.new

    # check for convergence
    if(abs(convcrit) < control$tol) break
    if(iteration >= control$max.iter) break
  }


  boundary <- theta[1:2] < 1e-4 * max(theta[1:2])
  if(any(boundary)){
    warning("variance component(s) ", paste(c("sd_b1","sd_b2")[boundary], collapse=", "),
            " collapsed to the zero boundary; this fit is degenerate", call. = FALSE)
  }

  list(a=a, alpha=alpha, knotseq=knotseq, col_means=col_means, theta=theta, modes=modes, hessians=hessians,
       cache=cache, logLik=logLik.new, iterations=iteration, boundary=boundary)
}

# Gauss-Hermite quadrature grid
qp <- 5 # 25 total nodes
gh <- statmod::gauss.quad(qp, "hermite")
AGQ_NODES <- as.matrix(expand.grid(z1=gh$nodes, z2=gh$nodes))
AGQ_LOGW  <- with(expand.grid(w1=gh$weights, w2=gh$weights), log(w1)+log(w2))
AGQ_ZSQ   <- rowSums(AGQ_NODES^2)   # constant across MarginalLogLik() calls

subj_curve_fun <- function(b, alpha_, a_, knotseq, col_means){
  # function computing the subject curve given the fixed and random parameters from fitting
  #
  # b: subject's random-effects vector (b1,b2)
  # alpha_, a_: fixed intercept / spline coefficients
  # knotseq, col_means: spline knot sequence / basis column means (centering)
  #
  # `t_grid` is resolved from the global environment, so it must be defined
  # before this is called
  u <- t_grid - b[2]
  alpha_ + b[1] + as.numeric(sweep(splines::splineDesign(knotseq, u, 4, outer.ok=TRUE), 2, col_means, FUN = "-") %*% a_)
}

# numerically estimated covariance matrices can end up just barely non-positive-definite
# make the smallest eigenvalue a small positive eps
#
# Sigma: covariance matrix to fix; eps: minimum eigenvalue floor
make_pd <- function(Sigma, eps = 1e-8){
  Sigma <- (Sigma + t(Sigma))/2
  ev <- eigen(Sigma, symmetric = TRUE)
  if(min(ev$values) < eps){
    Sigma <- Sigma + diag(eps - min(ev$values), nrow(Sigma))
  }
  Sigma
}


nk_grid <- 0:6   # interior knots considered for every fit; K = nk + 4

#' Fit snmmAGQ at every knot count in the grid and keep the AIC-best fit
#' @param dat full dataset (tibble)
#' @param y_list,t_list per-subject response / covariate lists
#' @param nk_values interior-knot counts to try
#' @param control EM controls, passed straight to bspline.snm()
#' @return the selected bspline.snm() fit

fit_select_nk <- function(dat, y_list, t_list, nk_values = nk_grid,
                          control = list(tol=.001, max.iter=15, inner.iterlim=100)){

  fits <- vector("list", length(nk_values))
  aics <- rep(NA_real_, length(nk_values))
  t_start <- rep(NA_real_, length(nk_values))
  t_fit <- rep(NA_real_, length(nk_values))
  last_error <- NULL

  for(k in seq_along(nk_values)){

    nk_k <- nk_values[k]

    f <- tryCatch({
      t0 <- Sys.time()
      start <- startingvalues(dat, knots=list(TRUE, nk_k))
      t1 <- Sys.time()
      ff <- bspline.snm(dat, y_list, t_list, start=start,
                        knots=list(TRUE, nk_k), control=control)
      ff$time_start <- as.numeric(difftime(t1, t0, units = "secs"))
      ff$time_fit <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
      ff
    }, error = function(e){
      last_error <<- conditionMessage(e)
      message("    nk = ", nk_k, " failed: ", conditionMessage(e))
      NULL
    })

    if(is.null(f)) next

    fits[[k]] <- f
    t_start[k] <- f$time_start
    t_fit[k] <- f$time_fit
    aics[k] <- 2 * (-f$logLik + length(f$a))
  }

  if(all(is.na(aics))){
    stop("every knot count in the grid failed (last error: ",
         if(is.null(last_error)) "none recorded" else last_error, ")")
  }

  best <- which.min(aics)
  fit <- fits[[best]]

  fit$nk_selected <- nk_values[best]
  fit$nk_grid <- nk_values
  fit$nk_aic <- aics
  fit$nk_ok <- sum(!is.na(aics))
  fit$nk_time_start <- t_start
  fit$nk_time_fit <- t_fit
  fit
}
