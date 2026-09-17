# ------------------------------------------------------------------------------
# sinecurve_scalability_merge_and_plots.R
#
# Merges the per-job files written by sinecurve_scalability.R and draws Figure 3
# of manuscript Section 3.
#
# Run this after all jobs have finished, from this directory:
#
#   Rscript sinecurve_scalability_merge_and_plots.R
#
# Reads  results/<prefix>_job<j>_of<n>.RDS
# Writes results/sinecurve_scalability_results.RDS   (merged, one row per cell x replication x method)
#        Figure3.pdf   median wall-clock time against n, m and K, all three methods
#
# ------------------------------------------------------------------------------

library(dplyr)
library(tibble)
library(ggplot2)

# ==============================================================================
# 1. Configuration
# ==============================================================================

jobs_dir <- "results"     # where the per-job files are, and where the merged
                          # RDS file is written

prefix <- "sinecurve_scalability"                 # output_prefix of sinecurve_scalability.R
merged_file <- paste0(prefix, "_results.RDS")

do_merge <- FALSE         # FALSE: skip the merge and read the already merged
                          # RDS file from jobs_dir instead

# ==============================================================================
# 2. Merge helpers
# ==============================================================================

#' @param prefix   output_prefix used by the job script
#' @param jobs_dir directory holding the *_job<j>_of<n>.RDS files
#' @return list(meta = the shared configuration, rows = one tibble of
#'   per-(cell, replication, method) rows across all jobs)
read_jobs <- function(prefix, jobs_dir) {

  files <- list.files(jobs_dir,
                      pattern = paste0("^", prefix, "_job[0-9]+_of[0-9]+\\.RDS$"),
                      full.names = TRUE)

  if (length(files) == 0) {
    stop("No job files matching ", prefix, "_job*_of*.RDS in ", jobs_dir)
  }

  jobs <- lapply(files, readRDS)

  # every file must agree on the split it belongs to
  key <- function(j) paste(j$njobs, j$config$iterations, j$config$seed_base, sep = "|")
  keys <- vapply(jobs, key, character(1))
  if (length(unique(keys)) > 1L) {
    stop("Job files for ", prefix, " disagree on njobs/iterations/seed_base:\n  ",
         paste(unique(keys), collapse = "\n  "),
         "\nThese come from different runs - merging them would corrupt the ",
         "results. Move the stale ones out of ", jobs_dir, ".")
  }
  for (k in seq_along(jobs)) {
    if (!identical(jobs[[k]]$cell_grid, jobs[[1]]$cell_grid)) {
      stop("Job file ", files[k], " was written with a different grid than ",
           files[1], ". Refusing to merge them.")
    }
  }

  njobs <- as.integer(jobs[[1]]$njobs)
  ids <- vapply(jobs, function(j) as.integer(j$job), integer(1))

  if (anyDuplicated(ids)) {
    stop("Duplicate job ids for ", prefix, ": ",
         paste(sort(ids[duplicated(ids)]), collapse = ", "))
  }
  if (length(ids) < njobs) {
    warning(prefix, ": only ", length(ids), " of ", njobs,
            " job files found - missing job(s) ",
            paste(setdiff(seq_len(njobs), ids), collapse = ", "),
            ". Their replications are absent (every job holds a share of every ",
            "cell, so the loss is spread evenly over the cells).",
            call. = FALSE)
  }

  # one row per (cell, replication, method), keeping the columns the figure needs
  rows <- bind_rows(lapply(jobs, function(j) {
    bind_rows(Filter(Negate(is.null), j$records))
  })) %>%
    select(sweep, cell, n, m, K, iter, method, ok, time_total)

  if (nrow(rows) == 0L) stop("The job files for ", prefix, " contain no records.")

  dup <- rows %>% count(cell, iter, method) %>% filter(n > 1L)
  if (nrow(dup) > 0L) {
    stop("Replication (cell ", dup$cell[1], ", iter ", dup$iter[1], ", ",
         dup$method[1], ") appears in more than one job file for ", prefix, ".")
  }

  list(
    meta = list(
      njobs = njobs,
      iterations = as.integer(jobs[[1]]$config$iterations),
      seed_base = jobs[[1]]$config$seed_base,
      cell_grid = jobs[[1]]$cell_grid,
      jobs_found = sort(ids)
    ),
    rows = rows
  )
}

# ==============================================================================
# 3. Merge
# ==============================================================================

if (!do_merge) {

  merged_path <- file.path(jobs_dir, merged_file)
  if (!file.exists(merged_path)) {
    stop("do_merge = FALSE but ", merged_path, " does not exist. ",
         "Set do_merge <- TRUE to build it from the per-job files.")
  }
  rows <- readRDS(merged_path)
  message("read: ", merged_path)

} else {

  jb <- read_jobs(prefix, jobs_dir)
  rows <- jb$rows

  message("njobs = ", jb$meta$njobs,
          " (files found: ", length(jb$meta$jobs_found), ")",
          ", iterations = ", jb$meta$iterations,
          ", seed_base = ", jb$meta$seed_base)

  report <- rows %>%
    group_by(sweep, cell, n, m, K, method) %>%
    summarize(
      attempts = n(),
      n_ok = sum(ok),
      n_failed = sum(!ok),
      median_time = median(time_total[ok]),
      .groups = "drop"
    )
  message("")
  print(as.data.frame(report), row.names = FALSE, digits = 3)

  saveRDS(rows, file = file.path(jobs_dir, merged_file))
  message("written: ", file.path(jobs_dir, merged_file))
}

# ==============================================================================
# 4. Figure 3: median wall-clock time against n, m and K, all three methods
# ==============================================================================

method_levels <- c("snmmTMB", "assist", "snmmAGQ")
method_colors <- c(snmmTMB = "#0046FA", assist = "#FF6242", snmmAGQ = "#A1E600")

sweep_labels <- c(
  subjects     = "Number of subjects",
  observations = "Observations per subject",
  basis        = "Basis functions"
)

plot_df <- rows %>%
  filter(ok) %>%
  mutate(
    x_value = case_when(
      sweep == "subjects"     ~ as.numeric(n),
      sweep == "observations" ~ as.numeric(m),
      sweep == "basis"        ~ as.numeric(K)
    )
  ) %>%
  group_by(sweep, cell, x_value, method) %>%
  summarize(
    time_median = median(time_total),
    time_q25 = quantile(time_total, 0.25),
    time_q75 = quantile(time_total, 0.75),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = method_levels),
    sweep = factor(sweep_labels[sweep], levels = unname(sweep_labels))
  )

figure3 <- ggplot(plot_df, aes(x = x_value, y = time_median, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = time_q25, ymax = time_q75), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ sweep, nrow = 1, scales = "free_x") +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = method_colors, name = "Method", drop = FALSE) +
  scale_fill_manual(values = method_colors, name = "Method", drop = FALSE) +
  labs(x = "Swept quantity (log scale)", y = "Time (seconds, log scale)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 0, b = 0),
    legend.box.spacing = unit(0.25, "lines"),
    legend.key.height = unit(0.9, "lines"),
    strip.text = element_text(size = 9.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 2, r = 6, b = 2, l = 2)
  )

ggsave("Figure3.pdf", figure3, width = 6.5, height = 3.1, units = "in")

message("\nwritten: Figure3.pdf")
