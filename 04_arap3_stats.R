# =====================================================================
# 04_arap3_stats.R — arap3 定量统计(供图注/正文与确定共定位用)
# 输出: tables/arap3_by_lineage.csv, tables/arap3_top_clusters.csv
# =====================================================================
## 运行前请将工作目录设为本仓库根目录(含 data/ figs/ tables/ 的文件夹)。
## RStudio: Session > Set Working Directory > To Source File Location;
## 或取消下一行注释、改成你的实际路径:
# setwd("path/to/Rscript")
suppressMessages({library(Seurat); library(dplyr)})
dir.create("tables", showWarnings = FALSE)

obj <- readRDS("data/zfdev.rds")
Idents(obj) <- "Cluster"
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
  grepl("^FinBud", ct) ~ "Fin bud",
  TRUE ~ "Other"))

v <- GetAssayData(obj, layer = "data")["arap3", ]
cat("arap3 overall: max log-expr =", round(max(v),2),
    "; % cells expressing =", round(mean(v>0)*100,2), "%\n")

by_lin <- data.frame(lineage = obj$lineage, expr = v) |>
  group_by(lineage) |>
  summarise(n_cells = n(), pct_expressing = round(mean(expr>0)*100,2),
            mean_logexpr = round(mean(expr),3)) |>
  arrange(desc(mean_logexpr))
write.csv(by_lin, "tables/arap3_by_lineage.csv", row.names = FALSE)
print(as.data.frame(by_lin))

by_clu <- data.frame(cluster = obj$Cluster, expr = v) |>
  group_by(cluster) |>
  summarise(n_cells = n(), pct_expressing = round(mean(expr>0)*100,2),
            mean_logexpr = round(mean(expr),3)) |>
  arrange(desc(mean_logexpr)) |> head(20)
write.csv(by_clu, "tables/arap3_top_clusters.csv", row.names = FALSE)
cat("\nsaved: tables/arap3_by_lineage.csv, tables/arap3_top_clusters.csv\n")
