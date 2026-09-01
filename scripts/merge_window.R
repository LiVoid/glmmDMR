#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stats)
})

#========================= Options =============================================
option_list <- list(

  # --- Input / Output ---
  make_option(c("--windows"), type = "character",
              help = "GLMM window-level results (tsv/tsv.gz). Required."),
  make_option(c("--out-prefix"), type = "character", default = "results/dmr",
              help = "Output prefix for results [default %default]"),

  # --- Mode ---
  make_option(c("--merge-mode"), type = "character", default = "hybrid_seed",
              help = "DMR merge strategy: single_seed, multi_seed, hybrid_seed [default %default]"),

  # --- Seed / Extension (core detection parameters) ---
  make_option(c("--p-seed"), type = "numeric", default = 0.05,
              help = "Seed p-value threshold for starting a DMR [default %default]"),
  make_option(c("--p-extend"), type = "numeric", default = 0.05,
              help = "Extend p-value threshold for growing a DMR [default %default]"),
  make_option(c("--max-gap-bp"), type = "integer", default = 200,
              help = "Max gap in bp to continue merging [default %default]"),
  make_option(c("--min-windows"), type = "integer", default = 1,
              help = "Minimum windows for a DMR [default %default]"),
  make_option(c("--min-delta"), type = "numeric", default = 0,
              help = "Minimum absolute delta for extension (effect size filter) [default %default]"),
  make_option(c("--max-p-degradation"), type = "numeric", default = 1.2,
              help = "Maximum allowed p-value degradation during extension (1.0=no degradation) [default %default]"),
  make_option(c("--max-final-p"), type = "numeric", default = 1.0,
              help = "Maximum final combined p-value for DMRs [default %default]"),
  make_option(c("--min-strong-windows"), type = "numeric", default = 0.5,
              help = "Minimum fraction of windows with p <= p-seed in final DMR [default %default]"),

  # --- Adaptive delta threshold ---
  make_option(c("--adaptive-delta"), action = "store_true", default = FALSE,
              help = "Use adaptive delta threshold based on candidate DMR distribution [default %default]"),
  make_option(c("--adaptive-delta-method"), type = "character", default = "median_ratio",
              help = "Adaptive delta method: median_ratio, q50, q25, q10, mad [default %default]"),
  make_option(c("--adaptive-delta-ratio"), type = "numeric", default = 0.6,
              help = "Ratio for median_ratio method (0.6 = 60%% of median) [default %default]"),

  # --- Multi-seed specific ---
  make_option(c("--seed-min-windows"), type = "integer", default = 1,
              help = "Minimum windows for multi-seed Stouffer method (1=single+multi, 2=multi only) [default %default]"),

  # --- Post-filter (consistency check after detection) ---
  make_option(c("--post-filter"), action = "store_true", default = FALSE,
              help = "Enable post-filtering to reduce false positives [default %default]"),
  make_option(c("--min-median-p"), type = "numeric", default = 0.01,
              help = "Maximum median p-value of windows in a DMR (post-filter) [default %default]"),
  make_option(c("--min-consistent-frac"), type = "numeric", default = 0.5,
              help = "Minimum fraction of windows with p <= p-seed (post-filter) [default %default]"),

  # --- Length / median-p filters ---
  make_option(c("--min-dmr-length"), type = "integer", default = 0,
              help = "Minimum DMR length in bp (0=no filter) [default %default]"),
  make_option(c("--max-median-p"), type = "numeric", default = 1.0,
              help = "Maximum median p-value for final DMRs [default %default]"),

  # --- Overlap merge (post-detection re-merge) ---
  make_option(c("--merge-overlaps"), action = "store_true", default = FALSE,
              help = "Merge overlapping/nearby same-direction DMRs after detection [default %default]"),
  make_option(c("--merge-overlaps-gap"), type = "integer", default = 0,
              help = "Max gap (bp) between DMRs to merge [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$windows)) stop("--windows argument is required.")

#========================= Helper Functions ====================================
safe_rbindlist <- function(lst) {
  if (length(lst) == 0) return(data.table())
  data.table::rbindlist(lst)
}

load_data <- function(path) {
  message("[INFO] Loading data: ", path)
  dt <- fread(path)
  stopifnot(c("chr", "start", "end", "p", "delta") %in% names(dt))
  dt <- dt[order(chr, start)]
  dt[, direction := ifelse(delta >= 0, "hyper", "hypo")]
  dt[is.finite(p) & p > 0 & p <= 1 & is.finite(delta)]
}

write_bed <- function(dmrs, out_path) {
  bed_data <- dmrs[, .(
    chr,
    start = as.integer(start - 1L),
    end = as.integer(end),
    score = as.integer(n_windows),
    direction
  )]
  fwrite(bed_data, out_path, sep = "\t", col.names = FALSE)
}

p_stouffer <- function(p_values) {
  p_values <- p_values[is.finite(p_values)]
  k <- length(p_values)
  if (k == 0) return(1)
  # Clip away from 0/1 so qnorm() never returns +-Inf: a p=0 (LRT statistic
  # underflow) combined with a p=1 (null-effect window) in the same run would
  # otherwise sum to Inf + (-Inf) = NaN and crash downstream if() checks.
  p_values <- pmin(pmax(p_values, 1e-15), 1 - 1e-15)
  z_scores <- qnorm(1 - p_values)
  combined_z <- sum(z_scores) / sqrt(k)
  2 * pnorm(-abs(combined_z))
}

add_delta_summary <- function(dmrs, win) {
  if (nrow(dmrs) == 0) return(dmrs)
  rows <- lapply(seq_len(nrow(dmrs)), function(i) {
    row <- dmrs[i, ]
    overlapping_windows <- win[chr == row$chr & start >= row$start & end <= row$end]
    row[, `:=`(
      delta_mean = mean(overlapping_windows$delta, na.rm = TRUE),
      delta_max = max(overlapping_windows$delta, na.rm = TRUE)
    )]
    row
  })
  safe_rbindlist(rows)
}

calculate_adaptive_delta_threshold <- function(dmrs, win, method = "median_ratio", ratio = 0.6) {
  if (nrow(dmrs) == 0) {
    message("[WARN] No DMRs found for adaptive threshold calculation. Using default min-delta.")
    return(0)
  }

  all_deltas <- unlist(lapply(seq_len(nrow(dmrs)), function(i) {
    row <- dmrs[i, ]
    overlapping_windows <- win[chr == row$chr & start >= row$start & end <= row$end]
    abs(overlapping_windows$delta)
  }))
  all_deltas <- all_deltas[is.finite(all_deltas) & all_deltas > 0]

  if (length(all_deltas) == 0) {
    message("[WARN] No valid delta values found. Using default min-delta.")
    return(0)
  }

  median_delta <- median(all_deltas, na.rm = TRUE)
  q25_delta <- as.numeric(quantile(all_deltas, 0.25, na.rm = TRUE))
  q10_delta <- as.numeric(quantile(all_deltas, 0.10, na.rm = TRUE))
  mad_delta <- mad(all_deltas, na.rm = TRUE)

  threshold <- switch(method,
    "median_ratio" = median_delta * ratio,
    "q50" = median_delta,
    "q25" = q25_delta,
    "q10" = q10_delta,
    "mad" = median_delta - mad_delta,
    median_delta * ratio
  )

  message(sprintf("[INFO] Adaptive delta threshold: %.4f (method=%s, median=%.4f, Q25=%.4f, Q10=%.4f)",
                  threshold, method, median_delta, q25_delta, q10_delta))
  max(0, threshold)
}

post_filter_dmrs <- function(dmrs, win, min_median_p, min_consistent_frac, p_seed) {
  if (nrow(dmrs) == 0) return(dmrs)

  message("[INFO] Applying post-filter: median_p <= ", min_median_p,
          ", consistent_frac >= ", min_consistent_frac)

  keep <- sapply(seq_len(nrow(dmrs)), function(i) {
    row <- dmrs[i, ]
    overlapping_windows <- win[chr == row$chr & direction == row$direction &
                                 start >= row$start & end <= row$end]
    if (nrow(overlapping_windows) == 0) return(FALSE)

    median_p <- median(overlapping_windows$p, na.rm = TRUE)
    consistent_frac <- mean(overlapping_windows$p <= p_seed, na.rm = TRUE)
    median_p <= min_median_p && consistent_frac >= min_consistent_frac
  })

  filtered <- dmrs[keep, ]
  message("[INFO] Post-filter: ", nrow(dmrs), " -> ", nrow(filtered),
          " DMRs (removed ", nrow(dmrs) - nrow(filtered), ")")
  filtered
}

merge_overlapping_dmrs <- function(dmrs, win, max_gap_bp = 0, p_seed = 0.05) {
  if (nrow(dmrs) == 0) return(dmrs)

  dmrs <- dmrs[order(chr, start, end)]
  merged <- list()
  i <- 1L

  while (i <= nrow(dmrs)) {
    current <- dmrs[i, ]
    j <- i + 1L
    while (j <= nrow(dmrs) && dmrs$chr[j] == current$chr && dmrs$direction[j] == current$direction) {
      gap <- dmrs$start[j] - current$end
      if (gap <= max_gap_bp) {
        current$end <- max(current$end, dmrs$end[j])
        j <- j + 1L
      } else {
        break
      }
    }

    merged_windows <- win[chr == current$chr & direction == current$direction &
                           start >= current$start & end <= current$end]
    if (nrow(merged_windows) == 0) {
      i <- j
      next
    }

    out_row <- data.table(
      chr = current$chr,
      start = current$start,
      end = current$end,
      n_windows = nrow(merged_windows),
      direction = current$direction,
      combined_p = p_stouffer(merged_windows$p)
    )
    if ("strong_frac" %in% names(dmrs)) {
      out_row[, strong_frac := mean(merged_windows$p <= p_seed, na.rm = TRUE)]
    }

    merged[[length(merged) + 1L]] <- out_row
    i <- j
  }

  safe_rbindlist(merged)
}

#========================= DMR Detection =======================================
detect_dmrs_single_seed <- function(win, p_seed, p_extend, max_gap_bp, min_windows, min_delta = 0,
                                    max_p_degradation = 1.2, max_final_p = 1.0, min_strong_windows = 0.5) {
  res <- list()
  i <- 1L
  n <- nrow(win)

  while (i <= n) {
    if (win$p[i] > p_seed) {
      i <- i + 1L
      next
    }

    chr <- win$chr[i]
    direction <- win$direction[i]
    run_idx <- i
    current_stouffer_p <- win$p[i]

    j <- i + 1L
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[max(run_idx)]
      if (gap > max_gap_bp) break
      if (min_delta > 0 && abs(win$delta[j]) < min_delta) {
        j <- j + 1L
        next
      }

      test_idx <- c(run_idx, j)
      stouffer_p <- p_stouffer(win$p[test_idx])
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        j <- j + 1L
      } else {
        break
      }
    }

    k <- i - 1L
    while (k >= 1 && win$chr[k] == chr && win$direction[k] == direction) {
      gap <- win$start[min(run_idx)] - win$end[k]
      if (gap > max_gap_bp) break
      if (min_delta > 0 && abs(win$delta[k]) < min_delta) {
        k <- k - 1L
        next
      }

      test_idx <- c(k, run_idx)
      stouffer_p <- p_stouffer(win$p[test_idx])
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        k <- k - 1L
      } else {
        break
      }
    }

    final_stouffer_p <- p_stouffer(win$p[run_idx])
    strong_windows_frac <- sum(win$p[run_idx] <= p_seed) / length(run_idx)

    if (length(run_idx) >= min_windows &&
        final_stouffer_p <= max_final_p &&
        strong_windows_frac >= min_strong_windows) {
      res[[length(res) + 1L]] <- list(
        chr = chr,
        start = min(win$start[run_idx]),
        end = max(win$end[run_idx]),
        n_windows = length(run_idx),
        direction = direction,
        combined_p = final_stouffer_p,
        strong_frac = strong_windows_frac
      )
    }

    i <- max(run_idx) + 1L
  }

  safe_rbindlist(res)
}

detect_dmrs_stouffer_multi_seed <- function(win, p_seed, p_extend, max_gap_bp, min_windows, min_delta = 0,
                                            max_p_degradation = 1.2, max_final_p = 1.0, min_strong_windows = 0.5,
                                            seed_min_windows = 2) {
  res <- list()
  i <- 1L
  n <- nrow(win)

  while (i <= n) {
    chr <- win$chr[i]
    direction <- win$direction[i]

    seed_idx <- i
    j <- i + 1L
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[j - 1L]
      if (gap > max_gap_bp) break
      seed_idx <- c(seed_idx, j)
      j <- j + 1L

      if (length(seed_idx) >= seed_min_windows) {
        stouffer_seed_p <- p_stouffer(win$p[seed_idx])
        if (stouffer_seed_p <= p_seed) break
      }
    }

    if (length(seed_idx) < seed_min_windows) {
      i <- max(seed_idx) + 1L
      next
    }

    stouffer_seed_p <- p_stouffer(win$p[seed_idx])
    if (stouffer_seed_p > p_seed) {
      i <- max(seed_idx) + 1L
      next
    }

    run_idx <- seed_idx
    current_stouffer_p <- stouffer_seed_p

    j <- max(run_idx) + 1L
    while (j <= n && win$chr[j] == chr && win$direction[j] == direction) {
      gap <- win$start[j] - win$end[max(run_idx)]
      if (gap > max_gap_bp) break
      if (min_delta > 0 && abs(win$delta[j]) < min_delta) {
        j <- j + 1L
        next
      }

      test_idx <- c(run_idx, j)
      stouffer_p <- p_stouffer(win$p[test_idx])
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        j <- j + 1L
      } else {
        break
      }
    }

    k <- min(seed_idx) - 1L
    while (k >= 1 && win$chr[k] == chr && win$direction[k] == direction) {
      gap <- win$start[min(run_idx)] - win$end[k]
      if (gap > max_gap_bp) break
      if (min_delta > 0 && abs(win$delta[k]) < min_delta) {
        k <- k - 1L
        next
      }

      test_idx <- c(k, run_idx)
      stouffer_p <- p_stouffer(win$p[test_idx])
      if (stouffer_p <= p_extend && stouffer_p <= current_stouffer_p * max_p_degradation) {
        run_idx <- test_idx
        current_stouffer_p <- stouffer_p
        k <- k - 1L
      } else {
        break
      }
    }

    final_stouffer_p <- p_stouffer(win$p[run_idx])
    strong_windows_frac <- sum(win$p[run_idx] <= p_seed) / length(run_idx)

    if (length(run_idx) >= min_windows &&
        final_stouffer_p <= max_final_p &&
        strong_windows_frac >= min_strong_windows) {
      res[[length(res) + 1L]] <- list(
        chr = chr,
        start = min(win$start[run_idx]),
        end = max(win$end[run_idx]),
        n_windows = length(run_idx),
        direction = direction,
        combined_p = final_stouffer_p,
        strong_frac = strong_windows_frac
      )
    }

    i <- max(run_idx) + 1L
  }

  safe_rbindlist(res)
}

detect_dmrs_hybrid_seed <- function(win, p_seed, p_extend, max_gap_bp, min_windows, min_delta = 0,
                                    max_p_degradation = 1.2, max_final_p = 1.0, min_strong_windows = 0.5,
                                    seed_min_windows = 2) {
  message("[INFO] Hybrid seed detection: multi_seed (primary) + single_seed (complementary)")

  dmrs_multi <- detect_dmrs_stouffer_multi_seed(win, p_seed, p_extend, max_gap_bp, min_windows,
                                                 min_delta, max_p_degradation, max_final_p,
                                                 min_strong_windows, seed_min_windows)
  if (nrow(dmrs_multi) > 0) {
    dmrs_multi[, detection_method := "multi_seed"]
  } else {
    dmrs_multi <- data.table(chr = character(), start = integer(), end = integer(),
                             n_windows = integer(), direction = character(),
                             combined_p = numeric(), strong_frac = numeric(),
                             detection_method = character())
  }

  covered_idx <- integer(0)
  if (nrow(dmrs_multi) > 0) {
    for (i in seq_len(nrow(dmrs_multi))) {
      dmr <- dmrs_multi[i, ]
      overlapping <- which(win$chr == dmr$chr &
                           win$direction == dmr$direction &
                           win$start >= dmr$start &
                           win$end <= dmr$end)
      covered_idx <- c(covered_idx, overlapping)
    }
    covered_idx <- unique(covered_idx)
  }

  remaining_win <- win[!seq_len(nrow(win)) %in% covered_idx]
  if (nrow(remaining_win) > 0) {
    dmrs_single <- detect_dmrs_single_seed(remaining_win, p_seed, p_extend, max_gap_bp, min_windows,
                                           min_delta, max_p_degradation, max_final_p, min_strong_windows)
    if (nrow(dmrs_single) > 0) {
      dmrs_single[, detection_method := "single_seed"]
    } else {
      dmrs_single <- data.table(chr = character(), start = integer(), end = integer(),
                                n_windows = integer(), direction = character(),
                                combined_p = numeric(), strong_frac = numeric(),
                                detection_method = character())
    }
  } else {
    dmrs_single <- data.table(chr = character(), start = integer(), end = integer(),
                              n_windows = integer(), direction = character(),
                              combined_p = numeric(), strong_frac = numeric(),
                              detection_method = character())
  }

  dmrs_hybrid <- rbind(dmrs_multi, dmrs_single)
  if (nrow(dmrs_hybrid) > 0) dmrs_hybrid <- dmrs_hybrid[order(chr, start)]
  dmrs_hybrid
}

#========================= Main Script =========================================
win <- load_data(opt$windows)

valid_modes <- c("single_seed", "multi_seed", "hybrid_seed")
if (!opt$`merge-mode` %in% valid_modes) {
  stop(sprintf("Invalid --merge-mode '%s'. Must be one of: %s",
               opt$`merge-mode`, paste(valid_modes, collapse = ", ")))
}

if (opt$`merge-mode` == "multi_seed") {
  if (opt$`adaptive-delta`) {
    message("[INFO] Using adaptive delta threshold...")
    dmrs_initial <- detect_dmrs_stouffer_multi_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                                     min_delta = 0, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                                     seed_min_windows = opt$`seed-min-windows`)
    adaptive_min_delta <- calculate_adaptive_delta_threshold(dmrs_initial, win,
                                                              method = opt$`adaptive-delta-method`,
                                                              ratio = opt$`adaptive-delta-ratio`)
    dmrs <- detect_dmrs_stouffer_multi_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                            adaptive_min_delta, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                            seed_min_windows = opt$`seed-min-windows`)
  } else {
    dmrs <- detect_dmrs_stouffer_multi_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                            opt$`min-delta`, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                            seed_min_windows = opt$`seed-min-windows`)
  }

  if (opt$`min-dmr-length` > 0) dmrs <- dmrs[end - start >= opt$`min-dmr-length`, ]

  if (opt$`max-median-p` < 1.0 && nrow(dmrs) > 0) {
    dmrs <- safe_rbindlist(lapply(seq_len(nrow(dmrs)), function(i) {
      row <- dmrs[i, ]
      overlapping_windows <- win[chr == row$chr & direction == row$direction &
                                   start >= row$start & end <= row$end]
      if (median(overlapping_windows$p, na.rm = TRUE) <= opt$`max-median-p`) row else NULL
    }))
  }

  if (opt$`post-filter`) dmrs <- post_filter_dmrs(dmrs, win, opt$`min-median-p`, opt$`min-consistent-frac`, opt$`p-seed`)
  if (opt$`merge-overlaps`) dmrs <- merge_overlapping_dmrs(dmrs, win, max_gap_bp = opt$`merge-overlaps-gap`, p_seed = opt$`p-seed`)

  dmrs <- add_delta_summary(dmrs, win)
  fwrite(dmrs, sprintf("%s_dmrs_multi_seed.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs, sprintf("%s_dmrs_multi_seed.bed", opt$`out-prefix`))

} else if (opt$`merge-mode` == "hybrid_seed") {
  if (opt$`adaptive-delta`) {
    message("[INFO] Hybrid seed with adaptive delta threshold...")
    dmrs_initial <- detect_dmrs_stouffer_multi_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                                     min_delta = 0, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                                     seed_min_windows = opt$`seed-min-windows`)
    adaptive_min_delta <- calculate_adaptive_delta_threshold(dmrs_initial, win,
                                                              method = opt$`adaptive-delta-method`,
                                                              ratio = opt$`adaptive-delta-ratio`)
    dmrs <- detect_dmrs_hybrid_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                    adaptive_min_delta, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                    seed_min_windows = opt$`seed-min-windows`)
  } else {
    dmrs <- detect_dmrs_hybrid_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                    opt$`min-delta`, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`,
                                    seed_min_windows = opt$`seed-min-windows`)
  }

  if (opt$`min-dmr-length` > 0) dmrs <- dmrs[end - start >= opt$`min-dmr-length`, ]

  if (opt$`max-median-p` < 1.0 && nrow(dmrs) > 0) {
    dmrs <- safe_rbindlist(lapply(seq_len(nrow(dmrs)), function(i) {
      row <- dmrs[i, ]
      overlapping_windows <- win[chr == row$chr & direction == row$direction &
                                   start >= row$start & end <= row$end]
      if (median(overlapping_windows$p, na.rm = TRUE) <= opt$`max-median-p`) row else NULL
    }))
  }

  if (opt$`post-filter`) dmrs <- post_filter_dmrs(dmrs, win, opt$`min-median-p`, opt$`min-consistent-frac`, opt$`p-seed`)
  if (opt$`merge-overlaps`) dmrs <- merge_overlapping_dmrs(dmrs, win, max_gap_bp = opt$`merge-overlaps-gap`, p_seed = opt$`p-seed`)

  dmrs <- add_delta_summary(dmrs, win)
  fwrite(dmrs, sprintf("%s_dmrs_hybrid_seed.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs, sprintf("%s_dmrs_hybrid_seed.bed", opt$`out-prefix`))

} else if (opt$`merge-mode` == "single_seed") {
  if (opt$`adaptive-delta`) {
    message("[INFO] Using adaptive delta threshold...")
    dmrs_initial <- detect_dmrs_single_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                            min_delta = 0, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`)
    adaptive_min_delta <- calculate_adaptive_delta_threshold(dmrs_initial, win,
                                                              method = opt$`adaptive-delta-method`,
                                                              ratio = opt$`adaptive-delta-ratio`)
    dmrs <- detect_dmrs_single_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                    adaptive_min_delta, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`)
  } else {
    dmrs <- detect_dmrs_single_seed(win, opt$`p-seed`, opt$`p-extend`, opt$`max-gap-bp`, opt$`min-windows`,
                                    opt$`min-delta`, opt$`max-p-degradation`, opt$`max-final-p`, opt$`min-strong-windows`)
  }

  if (opt$`min-dmr-length` > 0) dmrs <- dmrs[end - start >= opt$`min-dmr-length`, ]

  if (opt$`max-median-p` < 1.0 && nrow(dmrs) > 0) {
    dmrs <- safe_rbindlist(lapply(seq_len(nrow(dmrs)), function(i) {
      row <- dmrs[i, ]
      overlapping_windows <- win[chr == row$chr & direction == row$direction &
                                   start >= row$start & end <= row$end]
      if (median(overlapping_windows$p, na.rm = TRUE) <= opt$`max-median-p`) row else NULL
    }))
  }

  if (opt$`post-filter`) dmrs <- post_filter_dmrs(dmrs, win, opt$`min-median-p`, opt$`min-consistent-frac`, opt$`p-seed`)
  if (opt$`merge-overlaps`) dmrs <- merge_overlapping_dmrs(dmrs, win, max_gap_bp = opt$`merge-overlaps-gap`, p_seed = opt$`p-seed`)

  dmrs <- add_delta_summary(dmrs, win)
  fwrite(dmrs, sprintf("%s_dmrs_single_seed.tsv", opt$`out-prefix`), sep = "\t")
  write_bed(dmrs, sprintf("%s_dmrs_single_seed.bed", opt$`out-prefix`))
}

message("[INFO] DMR analysis complete.")
