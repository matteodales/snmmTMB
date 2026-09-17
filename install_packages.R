# ------------------------------------------------------------------------------
# install_packages.R
# ------------------------------------------------------------------------------

pkgs <- c(
  # model fitting
  "TMB",          # snmmTMB: automatic differentiation and Laplace approximation
  "nlme",         # starting values (nlme()), used inside assist and snmmAGQ
  "assist",       # snm(): Ke and Wang 2001
  "statmod",      # gauss.quad(): quadrature nodes for snmmAGQ (Elmi et al. 2011)
  "numDeriv",     # numerical Hessians for snmmAGQ
  "sitar",        # growth-curve application
  "MASS",         # mvrnorm()
  "Matrix",       # used by TMB
  # data handling and figures
  "dplyr", "tibble", "tidyr", "purrr", "forcats",
  "ggplot2", "patchwork", "scales", "ggrepel",
  # optional: single-threaded BLAS in the timing study (Figure 3)
  "RhpcBLASctl"
)

installed <- rownames(installed.packages())
to_install <- setdiff(pkgs, installed)

if (length(to_install) == 0) {
  message("all packages already installed")
} else {
  message("installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

still_missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(still_missing)) {
  stop("could not install: ", paste(still_missing, collapse = ", "))
}

message("done")
