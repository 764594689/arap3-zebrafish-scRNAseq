# =====================================================================
# 00_download.R  —  从 UCSC Cell Browser 下载 zebrafish-dev 数据
# 数据来源: Farnsworth et al. 2019/2020, A Single-Cell Transcriptome
#           Atlas for Zebrafish Development. PMID 31782996
# 在 RStudio 中: 先 setwd 到本脚本所在目录，或用 Session > Set Working Directory
# =====================================================================

dir.create("data", showWarnings = FALSE)

base <- "https://cells.ucsc.edu/zebrafish-dev/"
files <- c("exprMatrix.tsv.gz",          # 基因 x 细胞 表达矩阵(已 log 归一化)
           "meta.tsv",                    # 细胞元数据(含 ClusterNames / Cluster)
           "UMAP.coords.tsv.gz",          # UMAP 坐标
           "Seurat_tsne.coords.tsv.gz")   # tSNE 坐标(可选)

options(timeout = 3600)  # 95M 矩阵下载需放宽超时
for (f in files) {
  dest <- file.path("data", f)
  if (!file.exists(dest)) {
    message("downloading ", f, " ...")
    download.file(paste0(base, f), dest, mode = "wb")
  } else {
    message("exists, skip: ", f)
  }
}
message("done. files in data/:")
print(list.files("data"))
