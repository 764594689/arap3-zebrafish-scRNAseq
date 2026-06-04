# =====================================================================
# 03_coexpression_correlation.R — 数据驱动找 arap3 共表达基因(仅出表)
# 计算全体细胞中各基因与 arap3 的 Spearman 相关，输出 top 表。
# 不属于标准组合4图，仅作支撑数据/方法佐证。
# =====================================================================
## 运行前请将工作目录设为本仓库根目录(含 data/ figs/ tables/ 的文件夹)。
## RStudio: Session > Set Working Directory > To Source File Location;
## 或取消下一行注释、改成你的实际路径:
# setwd("path/to/Rscript")
library(Seurat)
dir.create("tables", showWarnings = FALSE)

obj  <- readRDS("data/zfdev.rds")
gene <- "arap3"
mat  <- GetAssayData(obj, slot = "data")       # genes x cells

## 只对在足够多细胞中表达的基因算相关，省时间 -----------------------
expr_frac <- Matrix::rowMeans(mat > 0)
keep <- names(expr_frac[expr_frac > 0.01])
v    <- mat[gene, ]

## Spearman 相关(对单细胞稀疏数据更稳健) ----------------------------
sub  <- as.matrix(mat[keep, ])
cors <- apply(sub, 1, function(x) suppressWarnings(cor(x, v, method = "spearman")))
res  <- data.frame(gene = names(cors), rho = cors)
res  <- res[order(-res$rho), ]
res  <- res[res$gene != gene, ]
write.csv(res, "tables/arap3_coexpression_correlation.csv", row.names = FALSE)

message("top correlated genes: ", paste(head(res$gene, 9), collapse = ", "))
message("done: tables/arap3_coexpression_correlation.csv")
