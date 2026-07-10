#!/usr/bin/env Rscript

Sys.setenv(FUTURE_AVAILABLECORES_METHODS = "system")

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(glmmTMB)
  library(future)
  library(future.apply)
})

# ---------- options ----------
opt_list <- list(
  make_option(c("-i","--infile"), type="character", help="input matrix (tsv[.gz]) with fixed columns"),
  make_option(c("-o","--out_prefix"), type="character", help="output prefix (dir/prefix)"),
  make_option(c("-c","--context"), type="character", default=NULL, help="filter by context (e.g., CG or CpG / CHG / CHH). NULL=use all"),
  make_option("--group1", type="character", help="name of group1 (e.g., WT)"),
  make_option("--group2", type="character", help="name of group2 (e.g., MT)"),
  make_option("--min_reps_g1", type="integer", default=2, help="minimum #replicates in group1 [default %default]"),
  make_option("--min_reps_g2", type="integer", default=2, help="minimum #replicates in group2 [default %default]"),
  make_option("--family", type="character", default="beta", help="model family: binom or beta [default %default]"),
  make_option("--mode", type="character", default="site", help="aggregate (window-sum per sample) or site (site-level) [default %default]"),
  make_option("--random_effect", action="store_true", default=TRUE, help="include (1|sample) [default %default]. Use --no-random_effect to disable"),
  make_option("--min_cov", type="integer", default=0, help="per-site minimum coverage (meth+unmeth) before window filters [default %default]"),
  make_option("--min_sites_win", type="integer", default=0, help="minimum #sites in a window to keep [default %default]"),
  make_option("--prefilter_delta", type="double", default=0.00, help="Quick screen by |Δ| >= threshold (0 disables) [default: %default]"),
  make_option("--workers", type="integer", default=4, help="#workers for future [default %default]"),
  make_option("--batches", type="integer", default=50, help="#batches (split windows) [default %default]"),
  make_option("--max_globals_mb", type="integer", default=1000, help="future.globals.maxSize (MiB) [default %default]"),
  make_option("--seed", type="integer", default=1, help="random seed [default %default]")
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$infile) || is.null(opt$out_prefix) || is.null(opt$group1) || is.null(opt$group2)) {
  stop("Required: --infile, --out_prefix, --group1, --group2")
}

set.seed(opt$seed)

# ---------- helpers ----------
COLS_FIXED <- c("chr","start","end",
                "cytosine_chr","cytosine_start","cytosine_end",
                "context","group","sample","strand","meth","unmeth","coverage")

fread_auto <- function(path) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("zcat", shQuote(path)))
  else fread(path)
}

mean_rate_by_group <- function(dt, gname) {
  dt[group==gname, mean(meth/(meth+unmeth), na.rm=TRUE)]
}

delta_group2_minus_group1 <- function(dt, g1, g2) {
  mean_rate_by_group(dt, g2) - mean_rate_by_group(dt, g1)
}

build_formulas <- function(family_tag = c("binom","beta"),
                           with_re = TRUE,
                           pooled = FALSE) {
  family_tag <- match.arg(family_tag)
  
  re_term <- if (with_re) "(1|sample)" else "1"
  add_re <- function(rhs) if (with_re) paste(rhs, "+ (1|sample)") else rhs
  
  if (family_tag == "binom") {
    list(
      family = "binom",
      alt  = as.formula( paste0("cbind(resp_m, resp_u) ~ ", add_re("group") ) ),
      null = as.formula( paste0("cbind(resp_m, resp_u) ~ ", add_re("1")    ) )
    )
  } else {
    list(
      family = "beta",
      alt  = as.formula( paste0("resp ~ ", add_re("group") ) ),
      null = as.formula( paste0("resp ~ ", add_re("1")    ) )
    )
  }
}

# single source of truth for the model label written to every output row
model_tag <- function(family, mode) {
  sprintf("ctx: GLMM %s (%s)",
          ifelse(family == "binom", "binomial", "beta"),
          mode)
}

# placeholder row returned when a window fails to fit
na_row <- function(chr, start, end, family, mode) {
  data.table(
    chr = chr, start = start, end = end,
    model = model_tag(family, mode),
    p = NA_real_, delta = NA_real_,
    mean_rate1 = NA_real_, mean_rate2 = NA_real_,
    aic_diff = NA_real_, bic_diff = NA_real_
  )
}

fit_one_window <- function(win_dt, fam=c("binom","beta"), mode=c("aggregate","site"), re=TRUE) {
  fam  <- match.arg(fam)
  mode <- match.arg(mode)
  # win_dt columns: chr, start, end, cytosine_*, group, sample, meth, unmeth, coverage
  # win_dt is a fresh, single-use subset (and workers get their own copy), so we
  # can modify it in place without copy()
  D <- win_dt
  D[, group := factor(group, levels=c(opt$group1, opt$group2))]

  pooled_flag <- (mode=="site")

  if (fam=="binom") {
    # use pseudocounts for model stability (does not affect filtering/output)
    D[, resp_m := meth + 0.1]
    D[, resp_u := unmeth + 0.1]
    fm <- build_formulas("binom", with_re=re, pooled=pooled_flag)
    fit_alt  <- glmmTMB(formula = fm$alt,  data = D, family = binomial)
    fit_null <- glmmTMB(formula = fm$null, data = D, family = binomial)
  } else if (fam == "beta") {
    # beta regression on rate; use meth/(meth+unmeth) to stay consistent with the
    # delta / mean_rate reported below (avoids divergence if an input `coverage`
    # column carries values other than meth+unmeth)
    D[, rate := meth / (meth + unmeth)]
    D[, resp := pmax(pmin(rate, 0.99), 0.01)]
    fm <- build_formulas("beta", with_re=re, pooled=pooled_flag)
    fit_alt  <- glmmTMB(formula=fm$alt,  family=beta_family(), data=D)
    fit_null <- glmmTMB(formula=fm$null, family=beta_family(), data=D)
  }

  # LRT
  llr <- suppressWarnings(anova(fit_null, fit_alt))
  pval <- tryCatch(llr$`Pr(>Chisq)`[2], error=function(e) NA_real_)

  # Δ = mean_rate(group2) − mean_rate(group1)
  dlt <- delta_group2_minus_group1(D, opt$group1, opt$group2)

  # info
  aic_diff <- AIC(fit_alt) - AIC(fit_null)
  bic_diff <- BIC(fit_alt) - BIC(fit_null)

  data.table(
#    chr = D$start[1] * 0 + D$chr[1],    # robust pull
    chr = D$chr[1],    # robust pull
    start = D$start[1], end = D$end[1],
    model = model_tag(fam, mode),
    p = pval, delta = dlt,
    mean_rate1 = mean_rate_by_group(D,opt$group1),
    mean_rate2 = mean_rate_by_group(D,opt$group2),
    aic_diff = aic_diff, bic_diff = bic_diff
  )
}

# ---------- load & filter ----------
DT <- fread_auto(opt$infile)
# attach header if missing
if (!all(COLS_FIXED %in% names(DT))) {
  # try to enforce fixed header if file was headerless in correct order
  if (ncol(DT) == length(COLS_FIXED)) {
    setnames(DT, COLS_FIXED)
  } else {
    stop("Input does not have required columns: ", paste(COLS_FIXED, collapse=", "))
  }
}

if (!is.null(opt$context)) {
  DT <- DT[context == opt$context]
}

# ensure numeric
num_cols <- c("start","end","cytosine_start","cytosine_end","meth","unmeth","coverage")
for (cc in num_cols) if (cc %in% names(DT)) DT[[cc]] <- as.numeric(DT[[cc]])

if (!("coverage" %in% names(DT))) DT[, coverage := meth + unmeth]
DT[is.na(coverage), coverage := meth + unmeth]

# site-level filter
if (opt$min_cov > 0L) {
  n_before <- data.table::uniqueN(DT[, .(chr, start, end)])
  DT <- DT[coverage >= opt$min_cov]
  n_after  <- data.table::uniqueN(DT[, .(chr, start, end)])
  message(sprintf("[site filter] kept %d / %d windows",
                  n_after, n_before))
}

# window site-count filter (after coverage filter)
if (opt$min_sites_win > 0L) {
  WIN_SITES <- DT[, .N, by = .(chr, start, end)]
  setnames(WIN_SITES, "N", "sites")
  KEEP_W <- WIN_SITES[sites >= opt$min_sites_win, .(chr, start, end)]
  setkey(KEEP_W, chr, start, end)
  setkey(DT, chr, start, end)
  n_before <- data.table::uniqueN(DT[, .(chr, start, end)])
  DT <- DT[KEEP_W, nomatch = 0L]
  n_after  <- data.table::uniqueN(DT[, .(chr, start, end)])
  message(sprintf("[coverage filter] kept %d / %d windows",
                  n_after, n_before))
  
}

# ---- replicate filter ----
if (opt$min_reps_g1 > 0L || opt$min_reps_g2 > 0L) {
  # 1) site count per window
  WIN_CHECK <- DT[, .(sites = .N), by = .(chr, start, end)]
  
  # 2) unique sample counts per window x group
  rep_by_win_grp <- DT[, .(
    reps = uniqueN(as.character(`sample`))
  ), by = .(chr, start, end, grp = as.character(`group`))]
  
  # 3) pivot to wide for group columns
  rep_wide <- data.table::dcast(
    rep_by_win_grp,
    chr + start + end ~ grp,
    value.var = "reps",
    fill = 0L
  )
  
  # 4) normalize group columns to reps_g1 / reps_g2
  if (opt$group1 %in% names(rep_wide)) data.table::setnames(rep_wide, opt$group1, "reps_g1") else rep_wide[, reps_g1 := 0L]
  if (opt$group2 %in% names(rep_wide)) data.table::setnames(rep_wide, opt$group2, "reps_g2") else rep_wide[, reps_g2 := 0L]
  
  # 5) join back to window table
  rep_wide <- rep_wide[, .(chr, start, end, reps_g1, reps_g2)]
  WIN_CHECK <- rep_wide[WIN_CHECK, on = .(chr, start, end)]
  
  KEEP_W <- WIN_CHECK[
    reps_g1 >= opt$min_reps_g1 &
    reps_g2 >= opt$min_reps_g2,
    .(chr,start,end)
  ]

  setkey(KEEP_W, chr,start,end)
  setkey(DT, chr,start,end)
  n_before <- data.table::uniqueN(DT[, .(chr, start, end)])
  DT <- DT[KEEP_W, nomatch = 0L]
  n_after  <- data.table::uniqueN(DT[, .(chr, start, end)])
  message(sprintf("[replicate filter] kept %d / %d windows",
                  n_after, n_before))
}

# ---- quick delta prefilter (optional; drop windows with small between-group difference) ----
if (!is.null(opt$prefilter_delta) && is.finite(opt$prefilter_delta) && opt$prefilter_delta >= 0) {
  
  # 1) sum meth/unmeth per window x group
  SAFE <- DT[, .(
    meth = sum(meth, na.rm = TRUE),
    unmeth = sum(unmeth, na.rm = TRUE)
  ), by = .(chr, start, end, grp = as.character(`group`))]
  
  # 2) coverage and rate
  SAFE[, cov := meth + unmeth]
  SAFE <- SAFE[cov > 0]
  SAFE[, rate := meth / cov]
  
  # 3) pivot to wide (mean rate per group)
  WIDE <- data.table::dcast(
    SAFE, chr + start + end ~ grp,
    value.var = "rate",
    fun.aggregate = mean,
    fill = NA_real_
  )
  
  # 4) ensure both group columns exist
  need_cols <- setdiff(c(opt$group1, opt$group2), names(WIDE))
  if (length(need_cols)) WIDE[, (need_cols) := NA_real_]
  
  # 5) approx delta = mean_rate(group2) - mean_rate(group1)
  WIDE[, delta_pref := get(opt$group2) - get(opt$group1)]
  
  # 6) keep windows
  if (opt$prefilter_delta == 0) {
    KEEP_W <- WIDE[is.finite(delta_pref) & delta_pref != 0, .(chr, start, end)]
  } else {
    KEEP_W <- WIDE[is.finite(delta_pref) & abs(delta_pref) >= opt$prefilter_delta, .(chr, start, end)]
  }
  
  # 7) filter DT by windows
  n_before <- data.table::uniqueN(DT[, .(chr, start, end)])
  data.table::setkey(KEEP_W, chr, start, end)
  data.table::setkey(DT,     chr, start, end)
  DT <- DT[KEEP_W, nomatch = 0L]
  n_after  <- data.table::uniqueN(DT[, .(chr, start, end)])
  message(sprintf("[delta filter] kept %d / %d windows (%s %.3g)",
                  n_after, n_before,
                  if (opt$prefilter_delta == 0) "Δ≠0" else "|Δ|≥",
                  if (opt$prefilter_delta == 0) 0 else opt$prefilter_delta))
}




# ---------- per-window preprocessing ----------
# aggregate mode: sum meth/unmeth/coverage per (window, sample) -> one row per sample
# site mode: pass site-level rows through unchanged
prep_window <- function(sub, mode = c("aggregate", "site")) {
  mode <- match.arg(mode)
  sub[, `:=`(sample = as.character(sample),
             group  = as.character(group))]
  if (mode == "aggregate") {
    sub <- sub[, .(
      meth     = sum(meth,     na.rm = TRUE),
      unmeth   = sum(unmeth,   na.rm = TRUE),
      coverage = sum(coverage, na.rm = TRUE),
      group    = unique(group)[1L]
    ), by = .(chr, start, end, sample)]
  }
  sub
}

# ---------- per-window runner ----------
# fit one window, returning an NA placeholder row on failure
fit_window_safe <- function(sub2) {
  tryCatch(
    fit_one_window(sub2, fam = opt$family, mode = opt$mode,
                   re = isTRUE(opt$random_effect)),
    error = function(e) na_row(sub2$chr[1L], sub2$start[1L], sub2$end[1L],
                               opt$family, opt$mode)
  )
}

# ---------- parallel ----------
plan(multisession, workers = opt$workers)
options(future.rng.onMisuse = "ignore",
        future.globals.maxSize = opt$max_globals_mb * 1024^2)

# ---------- build per-window subsets (single pass) ----------
# Preprocess the whole table once, then split it into per-window chunks. This
# replaces the previous O(nW x N) per-window logical scan of DT with a single
# O(N) grouped pass; for aggregate mode the per-sample summation is likewise done
# once for the whole table instead of repeatedly per window.
DT <- prep_window(DT, opt$mode)
SUB_ALL <- split(DT, by = c("chr", "start", "end"), sorted = FALSE)
rm(DT); invisible(gc(FALSE))

nW <- length(SUB_ALL)
if (nW == 0L) stop("No windows after filtering.")

# split window subsets into batches (roughly equal-sized chunks) to cap how many
# payloads are in flight to workers at once
split_ix <- parallel::splitIndices(nW, min(opt$batches, nW))

# process batches sequentially; fit windows within a batch in parallel
RES_LIST <- vector("list", length(split_ix))

for (bi in seq_along(split_ix)) {
  SUB_LIST <- SUB_ALL[ split_ix[[bi]] ]
  batch_res <- future.apply::future_lapply(
    SUB_LIST, fit_window_safe, future.seed = TRUE
  )
  RES_LIST[[bi]] <- data.table::rbindlist(batch_res, use.names = TRUE, fill = TRUE)
}
RES <- data.table::rbindlist(RES_LIST, use.names = TRUE, fill = TRUE)

# ---------- sort & write ----------
# chromosome order: numeric-like then others
setorderv(RES, c("chr","start","end"))

out_dir <- dirname(opt$out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
outfile <- paste0(opt$out_prefix, "_fit_", opt$family, "_", opt$mode, ".tsv.gz")
fwrite(RES, outfile, sep="\t")

message("Done: ", outfile)
