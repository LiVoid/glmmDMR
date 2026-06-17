#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(foreach)
})

option_list <- list(
  make_option(c("--input-dir"), type = "character", default = "output_for_MACAU",
              help = "Directory containing MACAU2 input BED files [default %default]"),
  make_option(c("--output-prefix"), type = "character", default = "output_for_MACAU/MACAU2",
              help = "Output file prefix [default %default]"),
  make_option(c("--macau2-r-dir"), type = "character", default = "",
              help = "Path to MACAU2 R directory containing MACAU2.R and RcppExports.R"),
  make_option(c("--group1-name"), type = "character", default = "WT",
              help = "Group 1 label [default %default]"),
  make_option(c("--group2-name"), type = "character", default = "MT",
              help = "Group 2 label [default %default]"),
  make_option(c("--rep-per-group"), type = "integer", default = 4,
              help = "Replicates per group [default %default]"),
  make_option(c("--max-sites"), type = "integer", default = 0,
              help = "Maximum number of sites (0 = all) [default %default]"),
  make_option(c("--merge-gap"), type = "integer", default = 300,
              help = "Maximum gap in bp to merge significant sites [default %default]"),
  make_option(c("--p-threshold"), type = "double", default = 0.05,
              help = "Site significance threshold [default %default]"),
  make_option(c("--num-core"), type = "integer", default = 1,
              help = "numCore passed to macau2 [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$`input-dir`)) {
  stop("Input directory not found: ", opt$`input-dir`)
}

if (!requireNamespace("MACAU2", quietly = TRUE)) {
  if (nzchar(opt$`macau2-r-dir`)) {
    source(file.path(opt$`macau2-r-dir`, "RcppExports.R"))
    source(file.path(opt$`macau2-r-dir`, "MACAU2.R"))
    if (!exists("macau2", mode = "function")) {
      stop("Could not load macau2() from --macau2-r-dir")
    }
  } else {
    stop("MACAU2 package not found. Install MACAU2 or set --macau2-r-dir")
  }
}

samples_g1 <- sprintf("%s%02d", opt$`group1-name`, seq_len(opt$`rep-per-group`))
samples_g2 <- sprintf("%s%02d", opt$`group2-name`, seq_len(opt$`rep-per-group`))
samples <- c(samples_g1, samples_g2)

files <- file.path(opt$`input-dir`, sprintf("windows_CG_forMACAU_%s.bed", samples))
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing MACAU2 input files:\n", paste(missing_files, collapse = "\n"))
}

all_sites <- lapply(files, function(f) {
  dt <- fread(f, header = FALSE,
              col.names = c("chr", "start", "end", "meth_rate", "coverage"))
  dt[, .(chr, start, end)]
})
unique_sites <- unique(rbindlist(all_sites))
setorder(unique_sites, chr, start, end)

if (!is.null(opt$`max-sites`) && opt$`max-sites` > 0) {
  unique_sites <- unique_sites[seq_len(min(opt$`max-sites`, nrow(unique_sites)))]
  cat(sprintf("[INFO] Limiting to %d sites\n", nrow(unique_sites)))
}

n_sites <- nrow(unique_sites)
n_samples <- length(files)
count_matrix <- matrix(0, nrow = n_sites, ncol = n_samples)
coverage_matrix <- matrix(0, nrow = n_sites, ncol = n_samples)

for (i in seq_along(files)) {
  dt <- fread(files[i], header = FALSE,
              col.names = c("chr", "start", "end", "meth_rate", "coverage"))
  dt[, meth_count := round(meth_rate * coverage)]
  merged <- merge(unique_sites, dt, by = c("chr", "start", "end"), all.x = TRUE)
  count_matrix[, i] <- merged$meth_count
  coverage_matrix[, i] <- merged$coverage
}

count_matrix[is.na(count_matrix)] <- 0
coverage_matrix[is.na(coverage_matrix)] <- 0

phenotypes <- c(rep(0, opt$`rep-per-group`), rep(1, opt$`rep-per-group`))
relatedness <- diag(length(phenotypes))

foreach::registerDoSEQ()

cat("Running MACAU2 in BMM mode\n")
if (exists("macau2", mode = "function")) {
  result <- macau2(
    RawCountDataSet = count_matrix,
    Phenotypes = phenotypes,
    Covariates = NULL,
    RelatednessMatrix = relatedness,
    LibSize = coverage_matrix,
    fit.model = "BMM",
    numCore = opt$`num-core`,
    filtering = FALSE,
    verbose = TRUE
  )
} else {
  result <- MACAU2::macau2(
    RawCountDataSet = count_matrix,
    Phenotypes = phenotypes,
    Covariates = NULL,
    RelatednessMatrix = relatedness,
    LibSize = coverage_matrix,
    fit.model = "BMM",
    numCore = opt$`num-core`,
    filtering = FALSE,
    verbose = TRUE
  )
}

result_with_pos <- cbind(unique_sites, as.data.table(result))
p_col <- intersect(c("pvalue", "p.value", "P", "p"), names(result_with_pos))
if (length(p_col) == 0) {
  stop("No p-value column found in MACAU2 result")
}
p_col <- p_col[1]

sig_sites <- result_with_pos[get(p_col) < opt$`p-threshold` & !is.na(get(p_col)), ]

dir.create(dirname(opt$`output-prefix`), recursive = TRUE, showWarnings = FALSE)
fwrite(result_with_pos, paste0(opt$`output-prefix`, "_sites.tsv"), sep = "\t")
fwrite(sig_sites, paste0(opt$`output-prefix`, "_sig_sites.tsv"), sep = "\t")

if (nrow(sig_sites) > 0) {
  sig_sites <- sig_sites[order(chr, start, end)]
  merged_list <- list()
  cur <- sig_sites[1, .(chr, start, end)]
  if (nrow(sig_sites) > 1) {
    for (i in 2:nrow(sig_sites)) {
      row <- sig_sites[i]
      if (row$chr == cur$chr && row$start <= (cur$end + opt$`merge-gap`)) {
        cur$end <- max(cur$end, row$end)
      } else {
        merged_list[[length(merged_list) + 1]] <- cur
        cur <- row[, .(chr, start, end)]
      }
    }
  }
  merged_list[[length(merged_list) + 1]] <- cur
  merged_dmrs <- rbindlist(merged_list)
  fwrite(merged_dmrs, paste0(opt$`output-prefix`, "_sig_windows_merged.tsv"), sep = "\t")
}

cat("MACAU2 analysis completed.\n")
