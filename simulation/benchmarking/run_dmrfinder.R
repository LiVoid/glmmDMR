#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--input-dir"), type = "character", default = "output_for_DMRfinder",
              help = "Directory containing DMRfinder input files [default %default]"),
  make_option(c("--python-bin"), type = "character", default = "python",
              help = "Python executable [default %default]"),
  make_option(c("--combine-script"), type = "character", default = "combine_CpG_sites.py",
              help = "Path to combine_CpG_sites.py [default %default]"),
  make_option(c("--finddmrs-script"), type = "character", default = "findDMRs_fixed.r",
              help = "Path to findDMRs_fixed.r [default %default]"),
  make_option(c("--group1-name"), type = "character", default = "WT",
              help = "Group 1 label [default %default]"),
  make_option(c("--group2-name"), type = "character", default = "MT",
              help = "Group 2 label [default %default]"),
  make_option(c("--rep-per-group"), type = "integer", default = 4,
              help = "Replicates per group [default %default]"),
  make_option(c("--r"), type = "integer", default = 3,
              help = "combine_CpG_sites.py -r argument [default %default]"),
  make_option(c("--s"), type = "integer", default = 4,
              help = "combine_CpG_sites.py -s argument [default %default]"),
  make_option(c("--merge-dist"), type = "integer", default = 300,
              help = "combine_CpG_sites.py -d argument [default %default]"),
  make_option(c("--min-cpg"), type = "integer", default = 3,
              help = "findDMRs -c argument [default %default]"),
  make_option(c("--delta"), type = "double", default = 0,
              help = "findDMRs -d argument [default %default]"),
  make_option(c("--append-dummy-row"), action = "store_true", default = TRUE,
              help = "Append one dummy row to avoid chr1-only edge case [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$`input-dir`)) {
  stop("Input directory not found: ", opt$`input-dir`)
}

samples_g1 <- sprintf("%s%02d", opt$`group1-name`, seq_len(opt$`rep-per-group`))
samples_g2 <- sprintf("%s%02d", opt$`group2-name`, seq_len(opt$`rep-per-group`))

file_list <- c(
  file.path(opt$`input-dir`, sprintf("sites_CG_forDMRfinder_%s.txt", samples_g1)),
  file.path(opt$`input-dir`, sprintf("sites_CG_forDMRfinder_%s.txt", samples_g2))
)

missing_files <- file_list[!file.exists(file_list)]
if (length(missing_files) > 0) {
  stop("Missing DMRfinder input files:\n", paste(missing_files, collapse = "\n"))
}

combine_args <- c(
  opt$`combine-script`,
  "-o", file.path(opt$`input-dir`, "results.csv"),
  "-r", as.character(opt$r),
  "-s", as.character(opt$s),
  "-d", as.character(opt$`merge-dist`),
  file_list
)

status <- system2(opt$`python-bin`, args = combine_args)
if (status != 0) {
  stop("combine_CpG_sites.py failed with exit status: ", status)
}

results_csv <- file.path(opt$`input-dir`, "results.csv")
results_mod_csv <- file.path(opt$`input-dir`, "results.mod.csv")

status <- system(sprintf("cat %s | sed -e 's/sites_CG_forDMRfinder_//g' > %s",
                         shQuote(results_csv), shQuote(results_mod_csv)))
if (status != 0) {
  stop("Failed to create results.mod.csv")
}

if (isTRUE(opt$`append-dummy-row`)) {
  dummy_cmd <- sprintf("echo -e \"chr2\\t1\\t2\\t0\\t$(yes 1 | head -n %d | paste -s)\" >> %s",
                       opt$`rep-per-group` * 4, shQuote(results_mod_csv))
  status <- system(dummy_cmd)
  if (status != 0) {
    stop("Failed to append dummy row")
  }
}

finddmrs_args <- c(
  opt$`finddmrs-script`,
  "-i", results_mod_csv,
  "-o", file.path(opt$`input-dir`, "out_findDMRs.txt"),
  "-n", sprintf("%s,%s", opt$`group1-name`, opt$`group2-name`),
  "-c", as.character(opt$`min-cpg`),
  "-d", as.character(opt$delta),
  paste(samples_g1, collapse = ","),
  paste(samples_g2, collapse = ",")
)

status <- system2("Rscript", args = finddmrs_args)
if (status != 0) {
  stop("findDMRs_fixed.r failed with exit status: ", status)
}

cat("DMRfinder analysis completed.\n")
cat("Output:", file.path(opt$`input-dir`, "out_findDMRs.txt"), "\n")
