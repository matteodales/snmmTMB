# ------------------------------------------------------------------------------
# bellcurve_server_compile.R
#
# Compiles the two TMB models used by the bell-curve simulation.
#
# Run this once on server before launching the job scripts.
#
# The first run takes a few minutes. After that, a model is recompiled only when
# its .cpp is newer than its shared library, since TMB::compile() has no
# up-to-date check of its own.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(TMB))

src_dir <- "../src"

cpp_files <- file.path(src_dir, c(
  "starting_points.cpp",
  "snmmTMB_likelihood_bellcurve_simulation.cpp"
))

# Compile a TMB model only when its .cpp is newer than its shared library
compile_if_stale <- function(cpp_file) {
  if (!file.exists(cpp_file)) {
    stop("Missing source file: ", cpp_file)
  }
  dll <- dynlib(sub("\\.cpp$", "", cpp_file))
  stale <- !file.exists(dll) ||
    file.info(dll)$mtime < file.info(cpp_file)$mtime
  if (stale) {
    message("Compiling ", cpp_file)
    status <- TMB::compile(cpp_file)
    if (!identical(as.integer(status), 0L) || !file.exists(dll)) {
      stop("Compilation of ", cpp_file, " failed (status ", status, ").")
    }
  } else {
    message("Up to date, not recompiling: ", dll)
  }
  invisible(dll)
}

for (f in cpp_files) {
  dll <- compile_if_stale(f)
  dyn.load(dll)
  message("OK: ", normalizePath(dll))
  dyn.unload(dll)
}

message("All models compiled.")
