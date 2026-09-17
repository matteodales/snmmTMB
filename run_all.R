# ------------------------------------------------------------------------------
# run_all.R
#
# Reproduces the figures and tables of the manuscript from the scripts in this
# folder. Before running, set the working directory to the one of this file.
#
# Every script is run in a new R session from its own folder. The run stops at the
# first script that fails.
# ------------------------------------------------------------------------------

# ==============================================================================
# 1. What to run
# ==============================================================================

figures     <- TRUE    # every figure and table, redrawn from the shipped results   ~5 min
examples    <- TRUE    # one replication of each simulation, all three methods     ~1 min
application <- FALSE   # refit the SMOCC data with snmmTMB and sitar, then Table 1 ~7 min
                       # and Figure 5 (rerun `figures` afterwards for Figure 6)
simulations <- FALSE   # the simulation studies behind Figures 1-4 and Supplementary
                       # Figures 1-2 from scratch, sequentially on this machine:
                       # about 2400 CPU-hours in total, see README.md for the
                       # cluster launchers and for a reduced run

# Two figures have their own recompute switch inside the script instead:
#   Supplementary Figure 3   run_simulation <- TRUE  in appendices/appendixC/bellcurve_edf_diagnostic.R   (~25 min)
#   Figure 6                 run_bootstrap  <- TRUE  in application_smocc/smocc_parambootstrap.R         (~7 h)

# ==============================================================================
# 2. Helper
# ==============================================================================

if (!file.exists("run_all.R")) {
  stop("set the working directory to the folder that contains run_all.R first")
}

rscript <- file.path(R.home("bin"), "Rscript")

#' Run one script in a fresh R session from its own folder
#' @param dir    folder of the script, relative to this one
#' @param script file name
#' @param args   command-line arguments passed to the script
run_script <- function(dir, script, args = character(0)) {
  old <- setwd(dir)
  on.exit(setwd(old))
  cat("\n---", format(Sys.time(), "%H:%M:%S"), file.path(dir, script), args, "\n")
  t0 <- Sys.time()
  status <- system2(rscript, c(script, args))
  if (status != 0) stop("script failed: ", file.path(dir, script))
  cat("--- done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
}

# ==============================================================================
# 3. Scripts
# ==============================================================================

if (examples) {
  # the plots go to Rplots.pdf in each folder
  run_script("simulation_sinecurve", "sinecurve_example.R")
  run_script("simulation_bellcurve", "bellcurve_example.R")
}

if (application) {
  run_script("application_smocc", "smocc_model_snmmTMB.R")   # results/smocc_snmmTMB_results.RDS
  run_script("application_smocc", "smocc_model_sitar.R")     # results/smocc_sitar_results.RDS
}

if (simulations) {
  # Each study is three job scripts (one per method), run here as a single job
  # each: Rscript <script> 1 1. On a cluster use the *_server_run_*.sh launchers.
  run_script("simulation_sinecurve", "sinecurve_server_compile.R")
  run_script("simulation_bellcurve", "bellcurve_server_compile.R")
  for (method in c("snmmTMB", "assist", "snmmAGQ")) {
    run_script("simulation_sinecurve", paste0("sinecurve_server_simulation_", method, ".R"), c(1, 1))
    run_script("simulation_bellcurve", paste0("bellcurve_server_simulation_", method, ".R"), c(1, 1))
    run_script("appendices/appendixB/t_residuals", paste0("sinecurve_server_simulation_t_residuals_", method, ".R"), c(1, 1))
    run_script("appendices/appendixB/t_randomeffects", paste0("sinecurve_server_simulation_t_randomeffects_", method, ".R"), c(1, 1))
  }
  run_script("simulation_sinecurve/scalability", "sinecurve_scalability.R", c(1, 1))
}

if (figures) {
  run_script("simulation_sinecurve", "sinecurve_server_merge_and_plots.R")                        # Figures 1, 2
  run_script("simulation_sinecurve/scalability", "sinecurve_scalability_merge_and_plots.R")       # Figure 3
  run_script("simulation_bellcurve", "bellcurve_server_merge_and_plots.R")                        # Figure 4
  run_script("application_smocc", "smocc_comparison_snmmTMB_vs_sitar.R")                          # Table 1, Figure 5
  run_script("application_smocc", "smocc_parambootstrap.R")                                       # Figure 6
  run_script("appendices/appendixB/t_residuals", "sinecurve_server_merge_and_plots_t_residuals.R")           # Supplementary Figure 1
  run_script("appendices/appendixB/t_randomeffects", "sinecurve_server_merge_and_plots_t_randomeffects.R")   # Supplementary Figure 2
  run_script("appendices/appendixC", "bellcurve_edf_diagnostic.R")                                # Supplementary Figure 3
  run_script("appendices/appendixD", "snmm_monotonesimulation.R")                                 # Supplementary Figure 4
}

cat("\nrun_all.R finished\n")
