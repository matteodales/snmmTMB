# snmmTMB
**Semi-parametric Non-linear Mixed Models via Template Model Builder (TMB)**

Code and data relative to the article "A semiparametric nonlinear mixed effects model with penalized splines using automatic differentiation".

The code in this repository can be used to fit semi-parametric non-linear mixed models (SNMMs) to longitudinal data. The population-level shape is captured by a penalized B-spline, subject-specific random effects allow the curve to shift, scale, and warp per individual, and all random effects (both subject-level and the spline's mixed-model random component) are integrated out analytically using the Laplace approximation implemented in [TMB](https://github.com/kaskr/adcomp).

To fit the model, the user needs to first write a C++ file that defines the joint log-likelihood for their specific model structure; TMB compiles it, applies the Laplace approximation and obtains derivatives through automatic differentiation (AD), returning a function that can be optimized in R.

## Table of Contents

1. [Repository structure](#1-repository-structure)
2. [Writing the C++ likelihood file](#2-writing-the-c-likelihood-file)
3. [Fitting the model in R](#3-fitting-the-model-in-r)
---


## 1. Repository structure

- **`simulations/`** — code to reproduce the simulation studies in the paper.
- **`application/`** — code to reproduce real-data application in the paper.
- **`snmmTMB_example.R`** — a self-contained example script showing how to set up, fit, and visualise results from the model.

---

## 2. Writing the C++ likelihood file

The C++ file is the core of the TMB workflow. It defines the joint log-likelihood that TMB marginalizes with the Laplace approximation and computes derivatives of. Below is a description of how to set it up, illustrated by the sine-curve example.

### 2.1 Headers and helper functions

```cpp
#include <TMB.hpp>
```

The example implements B-splines through Cox–de Boor recursion using `CppAD::CondExp*` comparison primitives (required for AD compatibility):

```cpp
template<class Type>
vector<Type> bspline_basis(Type x, vector<Type> knots, int degree) { ... }
```

If you want to use a different basis, you can define the function to evaluate it, making sure its AD-safe.

### 2.2 Declaring data inputs

Inside `objective_function<Type>::operator()`, declare every object passed from R using the appropriate TMB macro (e.g. `DATA_VECTOR`, `DATA_INTEGER`, `DATA_MATRIX`)

The example passes the response `y`, the time covariate `t`, subject group indices `group`, the knot vector, the penalty eigendecomposition (`Upos`, `U0`, `dpos`), and grid vectors for posterior curve evaluation.

### 2.3 Declaring parameters

Fixed effects and variance components use `PARAMETER`; for vectors of subject random effects use `PARAMETER_VECTOR`.

```cpp
PARAMETER(alpha);             // population intercept
PARAMETER_VECTOR(b1);         // random shifts
PARAMETER_VECTOR(b2);         // random log-amplitudes
PARAMETER_VECTOR(b3);         // random phases
PARAMETER(log_sd_b1);         // log standard deviations
PARAMETER(log_sd_b2);
PARAMETER(log_sd_b3);
PARAMETER(log_sigma);         // log residual SD
PARAMETER_VECTOR(vpos);       // penalized spline coordinates  ← declared random in R
PARAMETER_VECTOR(gamma0);     // unpenalized (null-space) coordinates
PARAMETER(log_lambda);        // log smoothing parameter
```

Add or remove random-effect vectors and variance parameters to match your specification.

### 2.4 Computing the negative log-likelihood

The `nll` includes three contributions:

1. **Observation likelihood** for all $N$ observations

2. **Random-effect priors** for all subjects

3. **Penalized spline prior** by placing a zero-mean multivariate normal on `vpos` with covariance $\text{diag}(1 / (\lambda \cdot d_k))$, where $d_k$ are the positive eigenvalues of the penalty matrix.

`nll` is the returned quantity from the C++ script.

### 2.5 Evaluating the population and subject curves



### 2.6 Reporting derived quantities

Use `ADREPORT` to flag quantities for which you want posterior means and delta-method standard errors:

```cpp
ADREPORT(h_grid);    // population curve evaluated on t_grid
ADREPORT(mu_sel);    // subject-specific curves for flagged subjects
```

Any expression computed from parameters can be `ADREPORT`-ed. The example only reports the population curve and the subject curves for the flagged subjects, both evaluated on the provided grid.

---

## 3. Fitting the model in R

### 3.1 — Compile and load

```r
library(TMB)
TMB::compile("your_model.cpp")
dyn.load(dynlib("your_model"))
```

### 3.2 — Construct the knot vector and penalty matrices

Choose the number of basis functions `K` and spline degree (typically `degree = 3` for cubic splines). Generate internal knots and extend them with repeated boundary knots:

```r
degree <- 3
K      <- 15
n_internal <- K - (degree + 1)
internal_knots <- seq(knot_lowlim, knot_uplim, length.out = n_internal + 2)
# ... extend with boundary knots (see snmmTMB_example.R for full code)
```

Build the second-difference penalty matrix, normalise it, and eigendecompose it to separate the penalized from the null-space directions:

```r
D    <- diff(diag(K), differences = 2)[, 1:(K-1)]  # sum-to-zero reduction
P    <- t(D) %*% D
P    <- P * qr(P)$rank / sum(diag(P))

ev      <- eigen(P, symmetric = TRUE)
Upos    <- ev$vectors[, ev$values > 1e-12]   # penalized directions
dpos    <- ev$values[ev$values > 1e-12]
U0      <- ev$vectors[, ev$values <= 1e-12]  # null-space directions
```

### 3.3 — Assemble the `Data` and `Params` lists

Pass all data objects to TMB exactly as declared in the C++ file:

```r
Data <- list(
  y = dat$y, t = dat$t,
  group = as.integer(dat$subject), nGroup = as.integer(n),
  subj_flag = as.integer(c(1, 1, rep(0, n-2))),
  knots = as.numeric(knot_vec), degree = as.integer(degree), K = as.integer(K),
  Upos = as.matrix(Upos), U0 = as.matrix(U0), dpos = as.numeric(dpos),
  spline_ci = as.integer(1), t_grid = as.numeric(t_grid)
)
```

Provide initial values for all parameters:

```r
Params <- list(
  alpha = 0, b1 = rep(0, n), b2 = rep(0, n), b3 = rep(0, n),
  log_sd_b1 = log(0.5), log_sd_b2 = log(0.5), log_sd_b3 = log(0.5),
  log_sigma = log(0.5),
  vpos = rep(0, length(dpos)), gamma0 = rep(0, ncol(U0)),
  log_lambda = log(1)
)
```

### 3.4 — Build the TMB objective and optimise

Declare which parameters are random (to be integrated out by Laplace):

```r
randoms <- c("b1", "b2", "b3", "vpos")

obj <- MakeADFun(data = Data, parameters = Params,
                 random = randoms, DLL = "your_model", silent = TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 1e4, iter.max = 1e4))
obj$par <- opt$par
```

### 3.5 — Extract inference

```r
rep <- sdreport(obj, getJointPrecision = TRUE)
print(summary(rep))
```

`sdreport` returns posterior means and delta-method standard errors for all fixed parameters and all `ADREPORT`-ed quantities.

---
