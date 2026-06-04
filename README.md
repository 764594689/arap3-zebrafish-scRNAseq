# arap3 co-expression analysis in the zebrafish single-cell atlas

Reproducible R code for the supplementary single-cell analysis showing that
***arap3*** is selectively co-expressed with endothelial markers, based on a
re-analysis of a public zebrafish whole-embryo scRNA-seq atlas.

斑马鱼单细胞图谱中 **arap3** 与内皮 marker 共表达分析的可复现 R 代码。

---

## Data source / 数据来源

- **Dataset:** A Single-Cell Transcriptome Atlas for Zebrafish Development
- **Reference:** Farnsworth DR, Saunders LM, Miller AC. *Dev. Biol.* 2020;459(2):100-108.
  **PMID: 31782996**
- **Accession:** NCBI BioProject **PRJNA564810**
- **Browser / download:** UCSC Cell Browser, dataset `zebrafish-dev`
  (https://cells.ucsc.edu/?ds=zebrafish-dev)
- 44,020 cells; 1, 2 and 5 dpf; 219 annotated clusters (`Cluster` field).
- Original UMAP/tSNE coordinates and cluster annotations from the authors are reused;
  no re-clustering or re-embedding is performed.

The expression matrix distributed by the UCSC Cell Browser is log-normalized;
gene identifiers are `ENSDARGxxxxxxxxxxx|symbol` and the symbol after `|` is used.

---

## Requirements / 运行环境

- **R 4.4.3** (Windows x86_64; also runs on macOS/Linux)
- CRAN packages (versions used here; see `sessionInfo.txt` for the full list):

  | package | version |
  |---|---|
  | Seurat | 5.2.1 |
  | SeuratObject | 5.0.2 |
  | data.table | 1.17.0 |
  | Matrix | 1.7-2 |
  | ggplot2 | 4.0.1 |
  | patchwork | 1.3.2 |
  | dplyr | 1.1.4 |
  | tidyr | 1.3.1 |

Install:

```r
install.packages(c("Seurat","data.table","Matrix","ggplot2",
                   "patchwork","dplyr","tidyr"))
```

> **Before running, set the working directory to this folder** (the repo root that
> contains `data/`, `figs/`, `tables/`). In RStudio: *Session > Set Working Directory >
> To Source File Location*, or uncomment and edit the `setwd(...)` line at the top of
> each script.

---

## Repository structure / 文件结构

```
Rscript/
├── 00_download.R                 # 从 UCSC Cell Browser 下载数据到 data/
├── 01_build_seurat.R             # 组装 Seurat 对象 -> data/zfdev.rds
├── 02_arap3_plots.R              # Fig1 (UMAP/FeaturePlot) + Fig4 (DotPlot)
├── 03_coexpression_correlation.R # 共表达相关性表 (支撑数据, 不入正图)
├── 04_arap3_stats.R              # arap3 定量统计表 (供图注/正文)
├── 05_standard_panels.R          # Fig2 (Heatmap) + Fig3 (StackedViolin)
├── README.md                     # 本文件
├── sessionInfo.txt               # 完整 R 会话/包版本
├── data/                         # 下载的原始文件 + zfdev.rds (体积大, 通常不随稿件上传)
├── figs/                         # 输出图 (PDF + TIFF)
└── tables/                       # 输出统计表 (CSV)
```

---

## How to run / 运行顺序

```r
source("00_download.R")                  # 下载数据 (~100 MB, 较慢)
source("01_build_seurat.R")              # 构建 Seurat 对象 (需 ~8 GB 内存)
source("02_arap3_plots.R")               # Fig1, Fig4
source("05_standard_panels.R")           # Fig2, Fig3
source("03_coexpression_correlation.R")  # 相关性表 (可选)
source("04_arap3_stats.R")               # 统计表 (可选)
```

---

## Outputs / 输出

**Figures (standard combination, panels a–d), each as PDF + 300-dpi RGB TIFF:**

| file | panel | content |
|---|---|---|
| `figs/Fig1_UMAP_arap3` | a | UMAP: lineages / endothelial highlight / *arap3* expression |
| `figs/Fig2_Heatmap`    | b | z-scored mean expression of *arap3* + 8 endothelial markers |
| `figs/Fig3_StackedViolin` | c | stacked violins of *arap3* + markers across lineages |
| `figs/Fig4_DotPlot`    | d | dot plot of *arap3* + markers across lineages |

**Tables:**

| file | content |
|---|---|
| `tables/arap3_by_lineage.csv` | per-lineage % expressing & mean expression |
| `tables/arap3_top_clusters.csv` | top clusters by *arap3* expression |
| `tables/arap3_coexpression_correlation.csv` | genome-wide Spearman correlation with *arap3* |

Endothelial markers used: *kdrl, cdh5, fli1a, flt1, tek, tie1, pecam1, cldn5b*.
The 219 clusters are grouped into 13 broad lineages (see the `case_when` block in the scripts).

---

## Suggested availability statement / 可用性声明模板

> **Data and resource availability.** This study re-analyses the publicly available
> zebrafish single-cell atlas of Farnsworth et al. (2020; NCBI BioProject PRJNA564810),
> obtained from the UCSC Cell Browser (dataset `zebrafish-dev`). All custom analysis
> code is available at [repository URL / DOI]. All other relevant data are within the
> article and its supplementary information.

The underlying dataset is subject to its original terms (Farnsworth et al., 2020).
