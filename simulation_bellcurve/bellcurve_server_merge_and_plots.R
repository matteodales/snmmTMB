# ------------------------------------------------------------------------------
# bellcurve_server_merge_and_plots.R
#
# Merges the per-job files written by the three bell-curve simulation scripts and
# draws Figure 4 of manuscript Section 3.
#
# Run this after all jobs have finished, from this directory:
#
#   Rscript bellcurve_server_merge_and_plots.R
#
# Reads  results/bellcurve_server_<method>_job<j>_of<n>.RDS
# Writes results/bellcurve_server_<method>_results.RDS   (merged, per method)
#        Figure4.pdf   population curve, all three methods
#
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(forcats)
library(ggplot2)
library(patchwork)
library(scales)   # percent() labels on the failure-rate panel

# ==============================================================================
# 1. Configuration
# ==============================================================================

jobs_dir <- "results"     # where the per-job files are, and where the merged
                          # per-method RDS files are written

do_merge <- FALSE         # if FALSE: skip the merge and read the already merged
                          # per-method RDS files from jobs_dir instead

methods <- tibble(method = c("snmmTMB", "assist", "snmmAGQ")) %>%
  mutate(prefix   = paste0("bellcurve_server_", method),
         out_file = paste0(prefix, "_results.RDS"))

# ==============================================================================
# 2. Merge helpers
# ==============================================================================

#' @param prefix   output_prefix used by the job script
#' @param jobs_dir directory holding the *_job<j>_of<n>.RDS files
#' @return list(meta = the shared configuration, records = one list of per-replication records across all jobs)
read_jobs <- function(prefix, jobs_dir) {

  files <- list.files(jobs_dir,
                      pattern = paste0("^", prefix, "_job[0-9]+_of[0-9]+\\.RDS$"),
                      full.names = TRUE)

  if (length(files) == 0) {
    stop("No job files matching ", prefix, "_job*_of*.RDS in ", jobs_dir)
  }

  jobs <- lapply(files, readRDS)

  # every file must agree on the split it belongs to
  key <- function(j) paste(j$method, j$njobs, j$iterations, j$seed_base,
                           paste(j$settings_to_run, collapse = ","), sep = "|")
  keys <- vapply(jobs, key, character(1))
  if (length(unique(keys)) > 1L) {
    stop("Job files for ", prefix, " disagree on method/njobs/iterations/",
         "seed_base/settings:\n  ", paste(unique(keys), collapse = "\n  "),
         "\nThese come from different runs - merging them would corrupt the ",
         "results. Move the stale ones out of ", jobs_dir, ".")
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
            ". Their replications will be counted as MISSING, not as failures.",
            call. = FALSE)
  }

  list(
    meta = list(
      method = jobs[[1]]$method,
      njobs = njobs,
      iterations = as.integer(jobs[[1]]$iterations),
      seed_base = jobs[[1]]$seed_base,
      settings_to_run = jobs[[1]]$settings_to_run,
      param_grid = jobs[[1]]$param_grid,
      jobs_found = sort(ids)
    ),
    records = unlist(lapply(jobs, `[[`, "records"), recursive = FALSE)
  )
}


#' Create per-setting results list
#' @param records list of list(setting_idx, iter, ok, metrics)
#' @param meta    shared configuration from read_jobs()
build_results <- function(records, meta) {

  iterations <- meta$iterations
  param_grid <- meta$param_grid
  settings <- seq_len(nrow(param_grid))

  ok_records <- Filter(function(r) isTRUE(r$ok), records)
  if (length(ok_records) == 0L) {
    stop("Not a single successful replication for method ", meta$method)
  }

  proto <- ok_records[[1]]$metrics
  fields <- names(proto)
  is_scalar <- vapply(proto, function(v) is.atomic(v) && length(v) == 1L, logical(1))

  empty_metrics <- function() {
    out <- vector("list", 0)
    for (f in fields) {
      out[[f]] <- if (is_scalar[[f]]) rep(NA, iterations) else vector("list", iterations)
    }
    out
  }

  metrics_by_setting <- lapply(settings, function(i) empty_metrics())
  present <- matrix(FALSE, nrow = length(settings), ncol = iterations)
  succeeded <- matrix(FALSE, nrow = length(settings), ncol = iterations)

  for (r in records) {
    if (is.null(r)) next
    s <- as.integer(r$setting_idx)
    it <- as.integer(r$iter)
    if (is.na(s) || is.na(it) || s > length(settings) || it > iterations) {
      warning("Record outside the grid (setting ", s, ", iter ", it, ") ignored.",
              call. = FALSE)
      next
    }
    if (present[s, it]) {
      stop("Replication (setting ", s, ", iter ", it, ") appears in more than ",
           "one job file for method ", meta$method, ".")
    }
    present[s, it] <- TRUE
    if (!isTRUE(r$ok)) next
    succeeded[s, it] <- TRUE
    m <- r$metrics
    for (f in fields) {
      if (is_scalar[[f]]) {
        metrics_by_setting[[s]][[f]][it] <- m[[f]]
      } else {
        metrics_by_setting[[s]][[f]][[it]] <- m[[f]]
      }
    }
  }

  results <- lapply(settings, function(i) {
    list(params = param_grid[i, ], metrics = metrics_by_setting[[i]])
  })

  report <- tibble(
    setting = settings,
    n = param_grid$n[settings],
    m = param_grid$m[settings],
    sigma = param_grid$sigma[settings],
    rho = param_grid$rho[settings],
    n_present = rowSums(present),
    n_ok = rowSums(succeeded),
    n_failed = rowSums(present) - rowSums(succeeded),
    n_missing = iterations - rowSums(present)
  )

  list(results = results, report = report)
}


# ==============================================================================
# 3. Merge every method
# ==============================================================================

merged <- list()

for (i in seq_len(nrow(methods))) {

  meth <- methods$method[i]
  message("\n=== ", meth, " ===")

  if (!do_merge) {
    merged_file <- file.path(jobs_dir, methods$out_file[i])
    if (!file.exists(merged_file)) {
      stop("do_merge = FALSE but ", merged_file, " does not exist. ",
           "Set do_merge <- TRUE to build it from the per-job files.")
    }
    merged[[meth]] <- list(results = readRDS(merged_file))
    message("read: ", merged_file)
    next
  }

  jb <- read_jobs(methods$prefix[i], jobs_dir)
  br <- build_results(jb$records, jb$meta)

  message("njobs = ", jb$meta$njobs,
          " (files found: ", length(jb$meta$jobs_found), ")",
          ", iterations = ", jb$meta$iterations,
          ", seed_base = ", jb$meta$seed_base)
  print(as.data.frame(br$report), row.names = FALSE)

  saveRDS(br$results, file = file.path(jobs_dir, methods$out_file[i]))
  message("written: ", file.path(jobs_dir, methods$out_file[i]))

  merged[[meth]] <- br
}



# ==============================================================================
# 4. Layout constants and palette
# ==============================================================================

group_spacing <- 12
setting_spacing <- 2.5
method_spacing <- 0.7

method_levels <- c("snmmTMB", "assist", "snmmAGQ")
method_colors <- c(snmmTMB = "#0046FA", assist = "#FF6242", snmmAGQ = "#A1E600")

# display order of the 8 settings, as used in the manuscript figure: the four
# n = 10 settings first, then the four n = 50 ones
new_order <- c(1, 3, 5, 7, 2, 4, 6, 8)

results_snmmTMB <- merged$snmmTMB$results[new_order]
results_assist  <- merged$assist$results[new_order]
results_snmmAGQ <- merged$snmmAGQ$results[new_order]

# ==============================================================================
# 5. Reshaping and plotting helpers
# ==============================================================================

#' Convert one method's `results` list into a long-format tibble
results_to_long <- function(results, method) {
  scalar_fields <- c(
    "pointwise_coverage", "simultaneous_coverage",
    "pointwise_CI_length", "simultaneous_CI_length",
    "subject_pointwise_coverage", "subject_simultaneous_coverage",
    "subject_pointwise_CI_length", "subject_simultaneous_CI_length",
    "b1_mse", "b2_mse",
    "b1_coverage", "b2_coverage",
    "rho_true", "rho_hat", "rho_se", "rho_bias", "rho_coverage",
    "time"
  )
  map_dfr(seq_along(results), function(i) {
    metrics <- results[[i]]$metrics
    niter <- length(metrics$pointwise_coverage)
    if (niter == 0) return(tibble())
    present <- intersect(scalar_fields, names(metrics))
    as_tibble(metrics[present]) %>%
      mutate(setting = factor(i), iter = seq_len(niter), method = method, .before = 1)
  })
}

#' Add a dodge x-position for grouped-by-setting, grouped-by-method plots
add_xpos <- function(df, group_size = 4) {
  df %>%
    mutate(
      setting_num = as.integer(as.character(setting)),
      group_id = ceiling(setting_num / group_size),
      setting_in_group = setting_num - (group_id - 1) * group_size,
      group_base = (group_id - 1) * group_spacing,
      setting_offset = (setting_in_group - 1) * setting_spacing,
      method_rank = as.integer(droplevels(method)) - mean(as.integer(droplevels(method))),
      xpos = group_base + setting_offset + method_rank * method_spacing
    )
}

#' Per-setting/method coverage, CI width and failure rate
summarize_coverage <- function(df, pointwise_cov, simultaneous_cov, pointwise_len, simultaneous_len) {
  df %>%
    group_by(setting, method) %>%
    summarize(
      n_total = n(),
      n_sim = sum(!is.na(.data[[simultaneous_cov]])),
      failure_rate = 1 - n_sim / n_total,
      x_sim = sum(.data[[simultaneous_cov]], na.rm = TRUE),
      mean_simultaneous = x_sim / n_sim,
      mean_sim_CI_len = mean(.data[[simultaneous_len]], na.rm = TRUE),
      sim_ci_lower = ifelse(x_sim == 0, 0, qbeta(0.025, x_sim, n_sim - x_sim + 1)),
      sim_ci_upper = ifelse(x_sim == n_sim, 1, qbeta(0.975, x_sim + 1, n_sim - x_sim)),
      sim_mcse = sqrt(mean_simultaneous * (1 - mean_simultaneous) / n_sim),
      .groups = "drop"
    )
}

#' The y scale shared by the coverage panels
coverage_y_scale <- function(...) {
  lo <- suppressWarnings(min(c(...), na.rm = TRUE))
  lower <- if (is.finite(lo)) max(0, floor(lo * 10) / 10) else 0.9
  scale_y_continuous(limits = c(lower, 1), breaks = seq(lower, 1, 0.1))
}

#' Coverage against setting, with binomial error bars and the nominal 0.95 line
plot_coverage <- function(summary_df,
                          y_col = "mean_simultaneous",
                          lo_col = "sim_ci_lower", hi_col = "sim_ci_upper",
                          ylab = "Sim. coverage") {
  summary_df <- add_xpos(summary_df)
  label_df <- summary_df %>%
    group_by(setting) %>%
    summarize(label_pos = mean(xpos), .groups = "drop") %>%
    arrange(as.numeric(as.character(setting)))

  ggplot(summary_df, aes(x = xpos, y = .data[[y_col]], fill = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
    geom_errorbar(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]), width = 0.5, colour = "black", linewidth = 0.6) +
    geom_line(aes(group = interaction(method, group_id)), colour = "black", linetype = "dotted", alpha = 0.6) +
    geom_point(shape = 21, size = 3.5, colour = "black", stroke = 0.4) +
    scale_fill_manual(values = method_colors, drop = TRUE) +
    scale_x_continuous(breaks = label_df$label_pos, labels = label_df$setting, expand = expansion(mult = c(0.02, 0.02))) +
    coverage_y_scale(summary_df[[y_col]], summary_df[[lo_col]], summary_df[[hi_col]]) +
    labs(x = "Setting", y = ylab, fill = "Method") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
}

#' Boxplot of a per-iteration metric against setting
plot_boxplot_by_setting <- function(df, value_col, ylab, ylim = NULL, log_y = FALSE,
                                    breaks = NULL) {
  df <- add_xpos(df)
  label_df <- df %>%
    group_by(setting) %>%
    summarize(label_pos = mean(xpos), .groups = "drop") %>%
    arrange(as.numeric(as.character(setting)))

  p <- ggplot(df, aes(x = xpos, y = .data[[value_col]], fill = method, group = interaction(xpos, method))) +
    # outlier.shape = NA: far points dropped, box and whiskers unchanged - with
    # 8 settings x 3 methods the dots swamped the boxes
    geom_boxplot(width = 0.95, linewidth = 0.3, outlier.shape = NA, alpha = 0.95) +
    scale_fill_manual(values = method_colors, drop = TRUE) +
    scale_x_continuous(breaks = label_df$label_pos, labels = label_df$setting, expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = "Setting", y = ylab, fill = "Method") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none", panel.grid.major.x = element_blank())

  if (log_y) {
    p <- p + if (is.null(breaks)) {
      scale_y_log10()
    } else {
      scale_y_log10(breaks = breaks, labels = breaks,
                    minor_breaks = sort(unique(c(breaks, breaks * 1.5))))
    }
  } else if (!is.null(breaks)) {
    p <- p + scale_y_continuous(breaks = breaks)
  }

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

#' Failure rate against setting
plot_failure_rate <- function(summary_df, rate_col = "failure_rate") {
  summary_df <- add_xpos(summary_df)
  label_df <- summary_df %>%
    group_by(setting) %>%
    summarize(label_pos = mean(xpos), .groups = "drop") %>%
    arrange(as.numeric(as.character(setting)))
  ylim <- c(0, max(0.05, max(summary_df[[rate_col]], na.rm = TRUE) * 1.15))

  ggplot(summary_df, aes(x = xpos, y = .data[[rate_col]], fill = method)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = method_colors, drop = TRUE) +
    scale_x_continuous(breaks = label_df$label_pos, labels = label_df$setting, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = ylim, labels = percent) +
    labs(x = "Setting", y = "Failure rate", fill = "Method") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
}

df_all <- bind_rows(
  results_to_long(results_snmmTMB, "snmmTMB"),
  results_to_long(results_assist, "assist"),
  results_to_long(results_snmmAGQ, "snmmAGQ")
) %>%
  mutate(method = factor(method, levels = method_levels))


# ==============================================================================
# 6. Figure 4: population curve, all three methods
# ==============================================================================

pop_summary <- summarize_coverage(df_all, "pointwise_coverage", "simultaneous_coverage",
                                  "pointwise_CI_length", "simultaneous_CI_length")

p_top <- plot_coverage(pop_summary)
p_bottom <- plot_boxplot_by_setting(df_all, "simultaneous_CI_length", "CI width",
                                    ylim = c(0, 2)) +
  guides(fill = "none")
p_failure <- plot_failure_rate(pop_summary) +
  guides(fill = "none")

figure4 <- (p_top / p_bottom / p_failure) +
  plot_layout(heights = c(0.4, 0.4, 0.2), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Figure4.pdf", figure4, width = 11, height = 8, units = "in")

message("\nwritten: Figure4.pdf")
