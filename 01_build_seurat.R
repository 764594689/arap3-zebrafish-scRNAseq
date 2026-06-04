# =====================================================================
# 01_build_seurat.R — 把 UCSC 文件组装成 Seurat 对象(沿用作者的聚类/坐标)
# 依赖: Seurat (>=4), data.table, Matrix
#   install.packages(c("Seurat","data.table","Matrix"))
# 输出: data/zfdev.rds
# =====================================================================
## 运行前请将工作目录设为本仓库根目录(含 data/ figs/ tables/ 的文件夹)。
## RStudio: Session > Set Working Directory > To Source File Location;
## 或取消下一行注释、改成你的实际路径:
# setwd("path/to/Rscript")
library(data.table)
library(Matrix)
library(Seurat)

## 1) 读表达矩阵 (第一列是基因名，其余列是细胞) -----------------------
# fread 解压 .gz 并快速读入；矩阵约 95M，需要 ~8G 内存
expr <- fread("data/exprMatrix.tsv.gz", sep = "\t", header = TRUE,
              data.table = FALSE)
genes <- expr[[1]]
expr  <- as.matrix(expr[, -1])
rownames(expr) <- genes
# 本数据基因名为 "ENSDARG...|symbol" 形式，取 | 后面的基因符号(如 arap3)
rownames(expr) <- sub("^.*\\|", "", rownames(expr))
rownames(expr) <- make.unique(rownames(expr))   # 防止符号重复
expr <- as(expr, "dgCMatrix")     # 稀疏化省内存

## 2) 读元数据 --------------------------------------------------------
meta <- read.table("data/meta.tsv", sep = "\t", header = TRUE,
                   row.names = 1, comment.char = "", quote = "",
                   check.names = FALSE)
# 对齐细胞顺序
common <- intersect(colnames(expr), rownames(meta))
expr <- expr[, common]
meta <- meta[common, , drop = FALSE]

## 3) 建 Seurat 对象。UCSC 矩阵通常已是 log 归一化值，放进 data 槽 ----
obj <- CreateSeuratObject(counts = expr, meta.data = meta)
obj <- SetAssayData(obj, slot = "data", new.data = expr)  # 直接用作 data

## 4) 设细胞类型为 Idents ---------------------------------------------
# 注意: 本数据 ClusterNames 只是 0~220 的编号；真正的细胞类型缩写在 Cluster 列
# (如 Endothe1, Vessel2, Peri1a, RBC1a, Macro5, FBProgAlla ...)
stopifnot("Cluster" %in% colnames(obj@meta.data))
Idents(obj) <- "Cluster"

## 5) 载入作者的 UMAP / tSNE 坐标 ------------------------------------
load_coords <- function(path, key) {
  # UCSC 坐标文件无表头: 列为 cell, x, y
  co <- read.table(path, sep = "\t", header = FALSE, row.names = 1,
                   check.names = FALSE)
  co <- as.matrix(co[colnames(obj), 1:2])
  colnames(co) <- paste0(key, "_", 1:2)
  CreateDimReducObject(embeddings = co, key = paste0(key, "_"),
                       assay = DefaultAssay(obj))
}
obj[["umap"]] <- load_coords("data/UMAP.coords.tsv.gz", "UMAP")
if (file.exists("data/Seurat_tsne.coords.tsv.gz"))
  obj[["tsne"]] <- load_coords("data/Seurat_tsne.coords.tsv.gz", "tSNE")

saveRDS(obj, "data/zfdev.rds")
message("Seurat object saved: data/zfdev.rds")
print(obj)
message("arap3 in matrix? ", "arap3" %in% rownames(obj))
