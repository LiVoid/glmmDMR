# run_DSS.R
suppressPackageStartupMessages({
  library(DSS); library(data.table)
})

# 入力ファイル
# df <- fread("/home/epigenome/dev_tools/glmm_dmrs/simulation_2nd/output_for_DSS/sites_CG_forDSS.tsv")
df <- fread("output_for_DSS/sites_CG_forDSS.tsv")
# MT01   MT02   MT03   MT04   WT01   WT02   WT03   WT04 
WT1 <- df[df$sample == "WT01",]
WT2 <- df[df$sample == "WT02",]
WT3 <- df[df$sample == "WT03",]
WT4 <- df[df$sample == "WT04",]

MT1 <- df[df$sample == "MT01",]
MT2 <- df[df$sample == "MT02",]
MT3 <- df[df$sample == "MT03",]
MT4 <- df[df$sample == "MT04",]

BSobj = makeBSseqData(list(WT1, WT2, WT3, WT4, MT1, MT2, MT3, MT4),
                      c("WT01","WT02","WT03","WT04",
                        "MT01","MT02","MT03","MT04"))

# グループ指定
dmlTest = DMLtest(BSobj, 
                  group1 = c("WT01","WT02","WT03","WT04"),
                  group2 = c("MT01","MT02","MT03","MT04"),
                  smoothing=TRUE, smoothing.span=500)
dmrs <- callDMR(dmlTest, p.threshold=0.05, minlen=500, minCG=3, dis.merge=200)

# 出力
write.table(as.data.frame(dmlTest), "output_for_DSS/DSS_dmlTest.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(dmrs), "output_for_DSS/DSS_dmrs.tsv", sep="\t", quote=FALSE, row.names=FALSE)