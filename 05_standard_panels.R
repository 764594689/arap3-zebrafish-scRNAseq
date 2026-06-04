# =====================================================================
# 05_standard_panels.R — 标准组合的两块:
#   Fig2 平均表达热图  scaled heatmap   (参考图1-b 风格)
#   Fig3 堆叠小提琴图  StackedVlnPlot   (参考图1-c 风格)
# 配合 02 的 Fig1(FeaturePlot, a) + Fig4(DotPlot, d) 即为完整标准组合(a~d)。
# 依赖: Seurat, ggplot2, dplyr, tidyr, patchwork
# =====================================================================
## 运行前请将工作目录设为本仓库根目录(含 data/ figs/ tables/ 的文件夹)。
## RStudio: Session > Set Working Directory > To Source File Location;
## 或取消下一行注释、改成你的实际路径:
# setwd("path/to/Rscript")
suppressMessages({library(Seurat); library(ggplot2); library(dplyr)
                  library(tidyr); library(patchwork)})
dir.create("figs", showWarnings = FALSE)

obj <- readRDS("data/zfdev.rds"); Idents(obj) <- "Cluster"
ct  <- as.character(obj$Cluster)
obj$lineage <- factor(dplyr::case_when(
  grepl("^Endothe|^Vessel", ct) ~ "Endothelial",
  grepl("^Peri|MusVascSmo", ct) ~ "Pericyte/Mural",
  grepl("^RBC", ct) ~ "Erythroid",
  grepl("^Macro|^Neutrophil|^Thymus|^Spleen", ct) ~ "Immune",
  grepl("^Mus|^Myo|Myotome|^Heart|NCcardiac", ct) ~ "Muscle/Heart",
  grepl("^Ret|^RPE|^Lens", ct) ~ "Retina/Eye",
  grepl("^NC", ct) ~ "Neural crest",
  grepl("^FB|^MB|^HB|^MHB|^SC|^antSC|RadialGlia|Oligo|Floorplate|Roofplate|Epiphysis|Olfactory|CranGang|HairCell|^Otic|LatLine|TasteBud", ct) ~ "Neural",
  grepl("^Basal|PharEpi|IntestineEpi|^Gill|^Iono", ct) ~ "Epithelial",
  grepl("^Hepato|^Panc|Intestine|^Kidney|^HG|Hypophysis|Hypochord", ct) ~ "Endoderm/Organ",
  grepl("cart|^Osteo|^Sclero|Notochord|Parachord", ct) ~ "Cartilage/Skeletal",
  grepl("^FinBud", ct) ~ "Fin bud", TRUE ~ "Other"))

gene    <- "arap3"
markers <- c("kdrl","cdh5","fli1a","flt1","tek","tie1","pecam1","cldn5b")
markers <- markers[markers %in% rownames(obj)]
feats   <- c(gene, markers)

## 让内皮排在最前，便于阅读
lin_order <- c("Endothelial","Pericyte/Mural","Muscle/Heart","Immune",
               "Erythroid","Neural","Neural crest","Retina/Eye","Epithelial",
               "Endoderm/Organ","Cartilage/Skeletal","Fin bud","Other")
obj$lineage <- factor(obj$lineage, levels = lin_order)

## ============ Fig3 堆叠小提琴图 (参考图1-c) =======================
p_vln <- VlnPlot(obj, features = feats, group.by = "lineage",
                 stack = TRUE, flip = TRUE, fill.by = "feature") +
         NoLegend() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               strip.text.y = element_text(face = "italic")) +  # 右侧基因名斜体
         ggtitle("Expression by lineage")
ggsave("figs/Fig3_StackedViolin.pdf", p_vln, width = 8, height = 7)
ggsave("figs/Fig3_StackedViolin.tiff", p_vln, width = 8, height = 7,
       dpi = 300, compression = "lzw")

## ============ Fig2 平均表达热图 (参考图1-b) =======================
# 每个 lineage 的平均 log 表达, 再按基因做 z-score, 画 tile 热图
dat <- GetAssayData(obj, layer = "data")[feats, ]
avg <- sapply(levels(obj$lineage), function(g)
         Matrix::rowMeans(dat[, obj$lineage == g, drop = FALSE]))
z   <- t(scale(t(avg)))                       # 按基因 z-score
hm  <- as.data.frame(z) |>
  tibble::rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "lineage", values_to = "z")
hm$gene    <- factor(hm$gene, levels = rev(feats))
hm$lineage <- factor(hm$lineage, levels = lin_order)

p_hm <- ggplot(hm, aes(lineage, gene, fill = z)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "z-score") +
  labs(x = NULL, y = NULL,
       title = "Scaled mean expression across lineages") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),   # 基因名斜体
        panel.grid = element_blank())
ggsave("figs/Fig2_Heatmap.pdf", p_hm, width = 8, height = 4)
ggsave("figs/Fig2_Heatmap.tiff", p_hm, width = 8, height = 4,
       dpi = 300, compression = "lzw")

cat("done: Fig2_Heatmap, Fig3_StackedViolin\n")
