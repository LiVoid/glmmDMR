#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--input-file"), type = "character", default = "output_for_metilene/sites_CG_forMetilene.txt",
              help = "Input metilene matrix file [default %default]"),
  make_option(c("--output-file"), type = "character", default = "output_for_metilene/metilene_out.tsv",
              help = "Output metilene result file [default %default]"),
  make_option(c("--filter-file"), type = "character", default = "output_for_metilene/sites_CG_forMetilene.filter",
              help = "Output metilene filter file [default %default]"),
  make_option(c("--group1-name"), type = "character", default = "WT",
              help = "Group 1 label [default %default]"),
  make_option(c("--group2-name"), type = "character", default = "MT",
              help = "Group 2 label [default %default]"),
  make_option(c("--metilene-bin"), type = "character", default = "metilene",
              help = "Path to metilene binary [default %default]"),
  make_option(c("--metilene-output-pl"), type = "character", default = "metilene_output.pl",
              help = "Path to metilene_output.pl [default %default]"),
  make_option(c("--M"), type = "integer", default = 300,
              help = "metilene -M argument [default %default]"),
  make_option(c("--m"), type = "integer", default = 3,
              help = "metilene -m argument [default %default]"),
  make_option(c("--d"), type = "double", default = 0,
              help = "metilene -d argument [default %default]"),
  make_option(c("--t"), type = "integer", default = 36,
              help = "metilene -t argument [default %default]"),
  make_option(c("--c"), type = "integer", default = 2,
              help = "metilene -c argument [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$`input-file`)) {
  stop("Input file not found: ", opt$`input-file`)
}

dir.create(dirname(opt$`output-file`), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$`filter-file`), recursive = TRUE, showWarnings = FALSE)

metilene_cmd <- sprintf(
  "%s -a %s -b %s -M %d -m %d -d %s -t %d -c %d %s > %s",
  shQuote(opt$`metilene-bin`),
  shQuote(opt$`group1-name`),
  shQuote(opt$`group2-name`),
  as.integer(opt$M),
  as.integer(opt$m),
  as.character(opt$d),
  as.integer(opt$t),
  as.integer(opt$c),
  shQuote(opt$`input-file`),
  shQuote(opt$`output-file`)
)

status <- system(metilene_cmd)
if (status != 0) {
  stop("metilene failed with exit status: ", status)
}

filter_args <- c(
  opt$`metilene-output-pl`,
  "-q", opt$`input-file`,
  "-o", opt$`filter-file`,
  "-a", opt$`group1-name`,
  "-b", opt$`group2-name`
)

status <- system2("perl", args = filter_args)
if (status != 0) {
  stop("metilene_output.pl failed with exit status: ", status)
}

cat("metilene analysis completed.\n")
cat("Output:", opt$`output-file`, "\n")
cat("Filter:", opt$`filter-file`, "\n")
