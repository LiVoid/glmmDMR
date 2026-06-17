#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)   # コマンドライン引数処理
  library(data.table) # データ操作
  library(ggplot2)    # プロット作成
})

#========================= Options =============================================
option_list <- list(
  make_option(c("--simes"), type="character",
              help="Input file: DMRs detected using Simes method (tsv format). Required."),
  make_option(c("--stouffer"), type="character",
              help="Input file: DMRs detected using Stouffer method (tsv format). Required."),
  make_option(c("--combined"), type="character",
              help="Input file: DMRs from combined method (tsv format). Required."),
  make_option(c("--out-prefix"), type="character", default="dmr_eval",
              help="Output prefix for evaluation results [default %default].")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Uncomment the lines below to test with local files (for development only)
# opt$simes <- "results/dmrs_simes.tsv"
# opt$stouffer <- "results/dmrs_stouffer.tsv"
# opt$combined <- "results/dmrs_combined.tsv"
# opt$`out-prefix` <- "results/dmr_comparison"

if (is.null(opt$simes) || is.null(opt$stouffer) || is.null(opt$combined)) {
  stop("[ERROR] --simes, --stouffer, and --combined arguments are required.")
}

#========================= Helper Functions ====================================
load_dmrs <- function(path, method) {
  message("[INFO] Loading DMRs from: ", path)
  dt <- fread(path)
  dt[, method := method] # 手法名を追加
  return(dt)
}

plot_dmr_distributions <- function(dmrs, out_prefix) {
  cat("[INFO] Generating DMR distribution plots\n")
  
  # 長さ分布のヒストグラム
  dmrs[, length_bp := end - start + 1]
  p1 <- ggplot(dmrs, aes(x = length_bp, fill = method)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_x_log10() +
    labs(title = "DMR Length Distribution", x = "Log10 Length (bp)", y = "Count") +
    theme_minimal()
  
  # スタート位置の分布
  p2 <- ggplot(dmrs, aes(x = start, fill = method)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    labs(title = "DMR Start Position Distribution", x = "Start Position", y = "Count") +
    theme_minimal()
  
  # PDFで保存
  ggsave(sprintf("%s_dmr_length_distribution.pdf", out_prefix), p1, width = 8, height = 5)
  ggsave(sprintf("%s_dmr_start_distribution.pdf", out_prefix), p2, width = 8, height = 5)
}

summarize_dmrs <- function(dmrs) {
  message("[INFO] Summarizing DMR results")
  summary_table <- dmrs[, .(
    num_dmrs = .N,
    mean_length = mean(end - start + 1),
    median_length = median(end - start + 1),
    total_length = sum(end - start + 1)
  ), by = .(method)]
  return(summary_table)
}

compare_methods <- function(dmrs_simes, dmrs_stouffer, dmrs_combined) {
  cat("[INFO] Comparing detection methods\n")
  
  # 方法のマージ
  all_dmrs <- rbind(dmrs_simes, dmrs_stouffer, dmrs_combined, fill = TRUE)
  
  # メソッドごとのユニークDMR数と重複領域数
  all_dmrs[, overlap := .N > 1, by = .(chr, start, end)]
  
  comparison_table <- all_dmrs[, .(
    unique_dmrs = sum(!overlap),
    overlapping_dmrs = sum(overlap)
  ), by = .(method)]
  return(comparison_table)
}

#========================= Main Script ==========================================
# Load DMRs
dmrs_simes <- load_dmrs(opt$simes, "Simes")
dmrs_stouffer <- load_dmrs(opt$stouffer, "Stouffer")
dmrs_combined <- load_dmrs(opt$combined, "Combined")

# Summarize results
summary_table <- summarize_dmrs(rbind(dmrs_simes, dmrs_stouffer, dmrs_combined, fill = TRUE))
#fwrite(summary_table, sprintf("%s_summary_table.tsv", opt$`out-prefix`), sep = "\t")

# Compare methods
comparison_table <- compare_methods(dmrs_simes, dmrs_stouffer, dmrs_combined)
#fwrite(comparison_table, sprintf("%s_method_comparison.tsv", opt$`out-prefix`), sep = "\t")

# Generate plots
plot_dmr_distributions(rbind(dmrs_simes, dmrs_stouffer, dmrs_combined, fill = TRUE), opt$`out-prefix`)
