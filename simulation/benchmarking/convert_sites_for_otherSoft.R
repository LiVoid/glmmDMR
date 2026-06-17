#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(tidyr); library(stringr) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: convert_sites_for_dss_metilene_methylkit.R <sites.tsv.gz>")
infile <- args[1]
#infile <- "/home/epigenome/dev_tools/glmm_dmrs/simulation_2nd/results/site_window_sim_forPaper/tsv/sites_CG.tsv.gz"
DT <- fread(infile)

# --- DSS format ---
dss_out <- DT[, .(chr, pos, N = meth + unmeth, X = meth, sample)]
dir.create("output_for_DSS", showWarnings=FALSE)
fwrite(dss_out, file = "output_for_DSS/sites_CG_forDSS.tsv", sep = "\t")


# --- methylKit format ---
df <- DT[, .(chrBase=paste(chr, pos, sep = "_"), chr, base=pos, strand="R", coverage = meth + unmeth, 
             freqC=(meth/(meth + unmeth)) * 100, freqT=(unmeth/(meth + unmeth))*100, sample)]
dir.create("output_for_methylKit", showWarnings=FALSE)



##         chrBase   chr    base strand coverage freqC  freqT
## 1 chr21.9764539 chr21 9764539      R       12 25.00  75.00
## 2 chr21.9764513 chr21 9764513      R       12  0.00 100.00
## 3 chr21.9820622 chr21 9820622      F       13  0.00 100.00
## 4 chr21.9837545 chr21 9837545      F       11  0.00 100.00
## 5 chr21.9849022 chr21 9849022      F      124 72.58  27.42

WT1 <- df[df$sample == "WT01", .(chrBase,chr,base,strand,coverage,freqC,freqT)]
WT2 <- df[df$sample == "WT02",.(chrBase,chr,base,strand,coverage,freqC,freqT)]
WT3 <- df[df$sample == "WT03",.(chrBase,chr,base,strand,coverage,freqC,freqT)]
WT4 <- df[df$sample == "WT04",.(chrBase,chr,base,strand,coverage,freqC,freqT)]

MT1 <- df[df$sample == "MT01",.(chrBase,chr,base,strand,coverage,freqC,freqT)]
MT2 <- df[df$sample == "MT02",.(chrBase,chr,base,strand,coverage,freqC,freqT)]
MT3 <- df[df$sample == "MT03",.(chrBase,chr,base,strand,coverage,freqC,freqT)]
MT4 <- df[df$sample == "MT04",.(chrBase,chr,base,strand,coverage,freqC,freqT)]


fwrite(WT1, file = "output_for_methylKit/sites_CG_forMethylKit_WT1.txt", sep = "\t")
fwrite(WT2, file = "output_for_methylKit/sites_CG_forMethylKit_WT2.txt", sep = "\t")
fwrite(WT3, file = "output_for_methylKit/sites_CG_forMethylKit_WT3.txt", sep = "\t")
fwrite(WT4, file = "output_for_methylKit/sites_CG_forMethylKit_WT4.txt", sep = "\t")
fwrite(MT1, file = "output_for_methylKit/sites_CG_forMethylKit_MT1.txt", sep = "\t")
fwrite(MT2, file = "output_for_methylKit/sites_CG_forMethylKit_MT2.txt", sep = "\t")
fwrite(MT3, file = "output_for_methylKit/sites_CG_forMethylKit_MT3.txt", sep = "\t")
fwrite(MT4, file = "output_for_methylKit/sites_CG_forMethylKit_MT4.txt", sep = "\t")

# --- metilene format ---

dir.create("output_for_metilene", showWarnings=FALSE)
#fwrite(meta_wide, file = "output_for_metilene/sites_CG_forMetilene.txt", sep = "\t", na = "NA")

# Calculate methylation percentage
DT[, rate := meth / (meth + unmeth)]
wide <- dcast(DT, chr + pos ~ sample, value.var = "rate")
wide[, `:=`(start = pos, end = pos)]
setcolorder(wide, c("chr", "pos", setdiff(names(wide), c("chr","pos","start", "end"))))
setorder(wide, chr, start)
out <- wide[, c(1:10)]

fwrite(out, file = "output_for_metilene/sites_CG_forMetilene.txt", sep = "\t", quote = FALSE, na = "NA")
#message("Written metilene input to ", outfile)




# --- DMRfinder format ---
#dmrf_out <- DT[, .(chr, pos, pos+1, rate = meth /(meth + unmeth), meth, unmeth, sample)]
#df <- DT[, .(chr, pos, meth, unmeth, sample)]
df <- DT[, .(chr, pos, end=pos+1, pct=meth /(meth + unmeth), meth, unmeth, sample)]
WT1 <- df[df$sample == "WT01",.(chr, pos, end, pct, meth, unmeth)]
WT2 <- df[df$sample == "WT02",.(chr, pos, end, pct, meth, unmeth)]
WT3 <- df[df$sample == "WT03",.(chr, pos, end, pct, meth, unmeth)]
WT4 <- df[df$sample == "WT04",.(chr, pos, end, pct, meth, unmeth)]

MT1 <- df[df$sample == "MT01",.(chr, pos, end, pct, meth, unmeth)]
MT2 <- df[df$sample == "MT02",.(chr, pos, end, pct, meth, unmeth)]
MT3 <- df[df$sample == "MT03",.(chr, pos, end, pct, meth, unmeth)]
MT4 <- df[df$sample == "MT04",.(chr, pos, end, pct, meth, unmeth)]



dir.create("output_for_DMRfinder", showWarnings=FALSE)
fwrite(WT1, file = "output_for_DMRfinder/sites_CG_forDMRfinder_WT1.txt", sep = "\t", col.names = FALSE)
fwrite(WT2, file = "output_for_DMRfinder/sites_CG_forDMRfinder_WT2.txt", sep = "\t", col.names = FALSE)
fwrite(WT3, file = "output_for_DMRfinder/sites_CG_forDMRfinder_WT3.txt", sep = "\t", col.names = FALSE)
fwrite(WT4, file = "output_for_DMRfinder/sites_CG_forDMRfinder_WT4.txt", sep = "\t", col.names = FALSE)

fwrite(MT1, file = "output_for_DMRfinder/sites_CG_forDMRfinder_MT1.txt", sep = "\t", col.names = FALSE)
fwrite(MT2, file = "output_for_DMRfinder/sites_CG_forDMRfinder_MT2.txt", sep = "\t", col.names = FALSE)
fwrite(MT3, file = "output_for_DMRfinder/sites_CG_forDMRfinder_MT3.txt", sep = "\t", col.names = FALSE)
fwrite(MT4, file = "output_for_DMRfinder/sites_CG_forDMRfinder_MT4.txt", sep = "\t", col.names = FALSE)


# --- MACAU2 format ---
dir.create("output_for_MACAU", showWarnings = FALSE)
macau_df <- DT[, .(chr, start = pos, end = pos + 1,
             meth_rate = meth / (meth + unmeth),
             coverage = meth + unmeth, sample)]

fwrite(macau_df[sample == "WT01", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_WT01.bed", sep = "\t", col.names = FALSE)
fwrite(macau_df[sample == "WT02", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_WT02.bed", sep = "\t", col.names = FALSE)
fwrite(macau_df[sample == "WT03", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_WT03.bed", sep = "\t", col.names = FALSE)
fwrite(macau_df[sample == "WT04", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_WT04.bed", sep = "\t", col.names = FALSE)

fwrite(macau_df[sample == "MT01", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_MT01.bed", sep = "\t", col.names = FALSE)
fwrite(macau_df[sample == "MT02", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_MT02.bed", sep = "\t", col.names = FALSE)
fwrite(macau_df[sample == "MT03", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_MT03.bed", sep = "\t", col.names = FALSE)
fwrite(macau_df[sample == "MT04", .(chr, start, end, meth_rate, coverage)],
    file = "output_for_MACAU/windows_CG_forMACAU_MT04.bed", sep = "\t", col.names = FALSE)
