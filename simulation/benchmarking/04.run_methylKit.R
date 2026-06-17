#!/usr/bin/env Rscript
# run_methylKit.R
suppressPackageStartupMessages({
  library(methylKit); library(data.table)
})

file.list <- list(
  "output_for_methylKit/sites_CG_forMethylKit_MT1.txt",
  "output_for_methylKit/sites_CG_forMethylKit_MT2.txt",
  "output_for_methylKit/sites_CG_forMethylKit_MT3.txt",
  "output_for_methylKit/sites_CG_forMethylKit_MT4.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT1.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT2.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT3.txt",
  "output_for_methylKit/sites_CG_forMethylKit_WT4.txt"
)
myobjDB <- methRead(
  file.list,
  sample.id = list("MT1","MT2","MT3","MT4",
                              "WT1","WT2","WT3","WT4"),
  assembly = "sim",
  treatment = c(1,1,1,1,0,0,0,0),
  context = "CpG",
  dbtype = "tabix",
  dbdir = "methylKit_DB"
)

filtered.myobj=filterByCoverage(myobjDB,lo.count=5,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

meth=unite(filtered.myobj, destrand=FALSE)

# 差分メチル化
#myDiff=calculateDiffMeth(meth)

#myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)


tiles <- tileMethylCounts(meth, win.size=300, step.size=200)
tiles_diff <- calculateDiffMeth(tiles,qvalue=0.05)


write.table(getData(tiles_diff), file="output_for_methylKit/methylKit_diff.tsv", sep="\t", quote=FALSE, row.names=FALSE)

