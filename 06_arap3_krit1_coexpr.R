# =====================================================================
# 06_arap3_krit1_coexpr.R — arap3 与 krit1 在内皮细胞的共表达展示
#   Fig5  FeaturePlot 并排 (arap3 | krit1)
#   Fig6  内皮细胞共表达状态 UMAP (4 类) + 堆叠柱状图
#   Fig7  内皮内 arap3 vs krit1 散点 (FeatureScatter 风格)
# 诚实呈现: arap3 内皮富集; krit1 广泛表达; 内皮内存在少量双阳细胞。
# 注意: krit1 低表达 + dropout, 双阳比例为下限估计, 不做"共富集"主张。
# 依赖: Seurat, ggplot2, patchwork, dplyr
# =====================================================================
## 运行前请将工作目录设为本仓库根目录(含 data/ figs/ tables/ 的文件夹)。
## RStudio: Session > Set Working Directory > To Source File Location;
## 或取消下一行注释、改成你的实际路径:
# setwd("path/to/Rscript")
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
dir.create("figs", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

obj  <- readRDS("data/zfdev.rds"); Idents(obj) <- "Cluster"
g1 <- "arap3"; g2 <- "krit1"
stopifnot(all(c(g1, g2) %in% rownames(obj)))

## ---- 系统归类(同前) ------------------------------------------------
ct <- as.character(obj$Cluster)
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

dat <- GetAssayData(obj, layer = "data")

## ====================================================================
## Fig5  FeaturePlot 并排 (arap3 | krit1)
## ====================================================================
ft <- function(g) FeaturePlot(obj, features = g, reduction = "umap",
                              order = TRUE, pt.size = 1.2, raster = TRUE,
                              raster.dpi = c(600, 600),
                              min.cutoff = 0, max.cutoff = "q99") +
  scale_color_gradientn(colours = c("grey88","#FED976","#FD8D3C",
                                    "#E31A1C","#800026")) +
  ggtitle(g) + theme(plot.title = element_text(face = "italic"))
p5 <- ft(g1) | ft(g2)
ggsave("figs/Fig5_FeaturePlot_arap3_krit1.pdf", p5, width = 12, height = 6)
ggsave("figs/Fig5_FeaturePlot_arap3_krit1.tiff", p5, width = 12, height = 6,
       dpi = 300, compression = "lzw")

## ====================================================================
## Fig6  内皮细胞共表达状态 (4 类) UMAP + 堆叠柱状图
## ====================================================================
pos1 <- dat[g1, ] > 0
pos2 <- dat[g2, ] > 0
status <- dplyr::case_when(
   pos1 &  pos2 ~ "arap3+krit1+",
   pos1 & !pos2 ~ "arap3+ only",
  !pos1 &  pos2 ~ "krit1+ only",
  TRUE          ~ "double-negative")
obj$coexp <- factor(status,
  levels = c("arap3+krit1+","arap3+ only","krit1+ only","double-negative"))

ec <- obj$lineage == "Endothelial"
ec_obj <- obj[, ec]
cols4 <- c("arap3+krit1+"="#D7191C","arap3+ only"="#FDAE61",
           "krit1+ only"="#2C7BB6","double-negative"="grey85")

# (左) 仅内皮细胞的 UMAP, 按共表达状态着色; 双阳点置顶
ord <- order(factor(ec_obj$coexp,
        levels = rev(c("arap3+krit1+","arap3+ only","krit1+ only","double-negative"))))
p6a <- DimPlot(ec_obj, reduction = "umap", group.by = "coexp",
               cols = cols4, order = TRUE, pt.size = 1.6, raster = TRUE) +
  ggtitle("Endothelial cells: arap3 / krit1 status") +
  theme(legend.position = "right")

# (右) 内皮细胞各状态比例 堆叠柱状图
df6 <- as.data.frame(prop.table(table(ec_obj$coexp)) * 100)
colnames(df6) <- c("status","pct")
df6$status <- factor(df6$status, levels = rev(levels(obj$coexp)))
p6b <- ggplot(df6, aes(x = "Endothelial", y = pct, fill = status)) +
  geom_col(width = .6, color = "white") +
  scale_fill_manual(values = cols4) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = .5), size = 3) +
  labs(x = NULL, y = "Percentage of endothelial cells", fill = NULL) +
  theme_classic()
p6 <- p6a + p6b + patchwork::plot_layout(widths = c(2, 1))
ggsave("figs/Fig6_endo_coexp_status.pdf", p6, width = 12, height = 5.5)
ggsave("figs/Fig6_endo_coexp_status.tiff", p6, width = 12, height = 5.5,
       dpi = 300, compression = "lzw")

## ====================================================================
## Fig7  内皮内 arap3 vs krit1 散点 (每个内皮细胞一点)
## ====================================================================
d7 <- data.frame(arap3 = dat[g1, ec], krit1 = dat[g2, ec],
                 status = ec_obj$coexp)
rho <- suppressWarnings(cor(d7$arap3, d7$krit1, method = "spearman"))
p7 <- ggplot(d7, aes(arap3, krit1, color = status)) +
  geom_jitter(width = .12, height = .12, size = 1.1, alpha = .8) +
  scale_color_manual(values = cols4) +
  labs(x = "arap3 (log-norm)", y = "krit1 (log-norm)",
       color = NULL,
       title = sprintf("Endothelial cells (n=%d), Spearman rho=%.2f",
                       sum(ec), rho)) +
  theme_classic() +
  theme(axis.title = element_text(face = "italic"))
ggsave("figs/Fig7_endo_scatter_arap3_krit1.pdf", p7, width = 7, height = 5.5)
ggsave("figs/Fig7_endo_scatter_arap3_krit1.tiff", p7, width = 7, height = 5.5,
       dpi = 300, compression = "lzw")

## ====================================================================
## 统计表 + Fisher 检验(供正文/图注引用; 双阳是否高于随机重叠)
## ====================================================================
a <- dat[g1, ec] > 0; k <- dat[g2, ec] > 0
tab <- table(arap3 = ifelse(a,"+","-"), krit1 = ifelse(k,"+","-"))
ftest <- fisher.test(tab)
sink("tables/endo_coexp_stats.txt")
cat("Endothelial cells n =", sum(ec), "\n")
cat(sprintf("arap3+ : %.1f%%\n", mean(a)*100))
cat(sprintf("krit1+ : %.1f%%\n", mean(k)*100))
cat(sprintf("arap3+krit1+ (observed): %.1f%%\n", mean(a&k)*100))
cat(sprintf("expected if independent: %.1f%%\n", mean(a)*mean(k)*100))
cat(sprintf("Fisher exact OR=%.2f, p=%.3g (95%%CI %.2f-%.2f)\n",
            ftest$estimate, ftest$p.value, ftest$conf.int[1], ftest$conf.int[2]))
cat("\n2x2 table:\n"); print(tab)
cat("\n注: krit1 低表达且存在 dropout, 双阳比例为下限估计; 共表达未显著高于随机重叠。\n")
sink()

message("done: Fig5/Fig6/Fig7 + tables/endo_coexp_stats.txt")
