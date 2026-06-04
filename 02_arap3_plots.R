# =====================================================================
# 02_arap3_plots.R — 标准组合: Fig1 (UMAP/FeaturePlot, panel a) +
#                    Fig4 (DotPlot, panel d)
# 依赖: Seurat, ggplot2, patchwork, dplyr
# 输出: figs/Fig1_UMAP_arap3.{pdf,tiff}, figs/Fig4_DotPlot.{pdf,tiff}
# =====================================================================
## 运行前请将工作目录设为本仓库根目录(含 data/ figs/ tables/ 的文件夹)。
## RStudio: Session > Set Working Directory > To Source File Location;
## 或取消下一行注释、改成你的实际路径:
# setwd("path/to/Rscript")
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
dir.create("figs", showWarnings = FALSE)

obj  <- readRDS("data/zfdev.rds")
gene <- "arap3"
stopifnot(gene %in% rownames(obj))

## ---- 候选共局在 marker(斑马鱼内皮系为主，按需增删) --------------
markers <- c("kdrl","cdh5","fli1a","flt1","tek","tie1","pecam1","cldn5b")
markers <- markers[markers %in% rownames(obj)]   # 只留数据中存在的

## ====================================================================
## 219 个簇过多：按 Cluster 缩写前缀折叠成"大类(lineage)"用于②③。
## 分类依据 Farnsworth 命名约定。case_when 自上而下，先匹配先生效。
## ====================================================================
Idents(obj) <- "Cluster"
ct <- as.character(obj$Cluster)
lineage <- dplyr::case_when(
  grepl("^Endothe|^Vessel", ct)                         ~ "Endothelial",
  grepl("^Peri|MusVascSmo", ct)                         ~ "Pericyte/Mural",
  grepl("^RBC", ct)                                     ~ "Erythroid",
  grepl("^Macro|^Neutrophil|^Thymus|^Spleen", ct)       ~ "Immune",
  grepl("^Mus|^Myo|Myotome|^Heart|NCcardiac", ct)       ~ "Muscle/Heart",
  grepl("^Ret|^RPE|^Lens", ct)                          ~ "Retina/Eye",
  grepl("^NC", ct)                                      ~ "Neural crest",
  grepl("^FB|^MB|^HB|^MHB|^SC|^antSC|RadialGlia|Oligo|Floorplate|Roofplate|Epiphysis|Olfactory|CranGang|HairCell|^Otic|LatLine|TasteBud", ct) ~ "Neural",
  grepl("^Basal|PharEpi|IntestineEpi|^Gill|^Iono", ct)  ~ "Epithelial",
  grepl("^Hepato|^Panc|Intestine|^Kidney|^HG|Hypophysis|Hypochord", ct) ~ "Endoderm/Organ",
  grepl("cart|^Osteo|^Sclero|Notochord|Parachord", ct)  ~ "Cartilage/Skeletal",
  grepl("^FinBud", ct)                                  ~ "Fin bud",
  TRUE ~ "Other")
obj$lineage <- factor(lineage)
cat("lineage counts:\n"); print(sort(table(obj$lineage), decreasing = TRUE))

## =================== 图① UMAP 三联(细胞型 + 内皮高亮 + arap3) =======
## 内皮细胞仅约 800/44020，所以(1)给 FeaturePlot 加大点径 order=TRUE 让阳性
## 细胞置顶可见; (2)增加一张"内皮系高亮"UMAP，方便读者对位。
endo_cells <- colnames(obj)[obj$lineage == "Endothelial"]
p_ct  <- DimPlot(obj, reduction = "umap", group.by = "lineage",
                 label = TRUE, repel = TRUE, raster = TRUE, label.size = 3) +
         ggtitle("Cell lineages") + theme(legend.position = "none")
p_hl  <- DimPlot(obj, reduction = "umap", cells.highlight = endo_cells,
                 cols.highlight = "#D7191C", sizes.highlight = 1.2,
                 raster = TRUE) + ggtitle("Endothelial cells") +
         theme(legend.position = "none")
p_g   <- FeaturePlot(obj, features = gene, reduction = "umap",
                     order = TRUE, pt.size = 1.4, raster = TRUE,
                     raster.dpi = c(600, 600),
                     min.cutoff = 0, max.cutoff = "q99") +  # 截断个别极高值
         # 灰->黄->橙->红->暗红 多段饱和色阶，阳性细胞立刻脱离灰色背景
         scale_color_gradientn(
           colours = c("grey88", "#FED976", "#FD8D3C",
                       "#E31A1C", "#800026")) +
         ggtitle("arap3") +
         theme(plot.title = element_text(face = "italic"))  # 基因名斜体
ggsave("figs/Fig1_UMAP_arap3.pdf", p_ct + p_hl + p_g, width = 18, height = 6)
ggsave("figs/Fig1_UMAP_arap3.tiff", p_ct + p_hl + p_g, width = 18, height = 6,
       dpi = 300, compression = "lzw")

## =================== 图④ DotPlot: arap3 + markers =================
## (标准组合 panel d；Fig2 热图、Fig3 小提琴由 05_standard_panels.R 生成)
p_dot <- DotPlot(obj, features = c(gene, markers), group.by = "lineage") +
  scale_color_gradient(low = "lightgrey", high = "red") +
  ggtitle("Marker expression by lineage") +
  # x 轴基因名斜体(同时保留旋转角度)
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
ggsave("figs/Fig4_DotPlot.pdf", p_dot, width = 8, height = 5)
ggsave("figs/Fig4_DotPlot.tiff", p_dot, width = 8, height = 5,
       dpi = 300, compression = "lzw")

message("figures written to figs/  (Fig1_UMAP_arap3, Fig4_DotPlot)")
