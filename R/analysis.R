# Preliminary analysis of the data

# Global options ---------------------------------------------------------------
ds <- "data/combined.Rds"
dest <- paste0("analysis/", "combined/", Sys.Date(), "/")
dir.create(dest, recursive = TRUE)

# Load data set ----------------------------------------------------------------
ds <- readRDS(ds)

# Quality control --------------------------------------------------------------

ds$percent.mito
mt_genes <- rownames(ds)[grep("^MT-", rownames(ds))]
ds$percent.mito <- round(
  sparseMatrixStats::colSums2(ds[["RNA"]]@counts[mt_genes, ]) /
    sparseMatrixStats::colSums2(ds[["RNA"]]@counts), 3
) * 100

# Plot
ggplot2::ggplot(
  ds@meta.data,
  ggplot2::aes(nCount_RNA, percent.mito, col = condition)
  ) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~condition) +
  ggplot2::theme_classic(20) +
  ggplot2::theme(
    legend.position = "",
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    )

fn <- paste0(dest, "qc", "_", "libsize", "-", "pMt", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5)

# Plot
ggplot2::ggplot(
  ds@meta.data,
  ggplot2::aes(nCount_RNA, nFeature_RNA, col = condition)
) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~condition) +
  ggplot2::theme_classic(20) +
  ggplot2::theme(
    legend.position = "",
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ggplot2::scale_x_continuous(trans = "log10") +
  ggplot2::scale_y_continuous(trans = "log10")

fn <- paste0(dest, "qc", "_", "libsize", "-", "features", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5)

# Normalization ----------------------------------------------------------------
ds <- Seurat::NormalizeData(ds)

# Feature Selection ------------------------------------------------------------
ds <- Seurat::FindVariableFeatures(ds, nfeatures = 5000)

# Linear dimensional reduction -------------------------------------------------

# Continue without scaling
hvgs <- ds[["RNA"]]@var.features
ds[["RNA"]]@scale.data <- as.matrix(ds[["RNA"]]@data[hvgs, ])

# Compute PCA
npcs <- c(10, 20, 30, 40, 60, 80)
for (i in npcs) {
  key <- paste0("N", length(hvgs), "_", "pca", i)
  ds <- Seurat::RunPCA(
    ds, npcs = i, reduction.name = key, reduction.key = paste0(key, "_")
    )
}

# fastMNN batch correction -----------------------------------------------------

for (i in npcs) {
  key <- paste0("N", length(hvgs), "_", "mnn", i)
  pca <- batchelor::fastMNN(
    ds[["RNA"]]@data[hvgs, ], batch = ds$condition, d = i
  )
  pca <- pca@int_colData$reducedDims$corrected
  colnames(pca) <- paste0(key, "_", 1:i)
  ds[[key]] <- Seurat::CreateDimReducObject(
    embeddings = pca, assay = "RNA"
  )
}

# Cell embedding ---------------------------------------------------------------

# Select reductions to use for cell embeddings
index <- names(ds@reductions)

for (i in index) {
  set.seed(42)
  umap <- uwot::umap(
    ds[[i]]@cell.embeddings,
    metric = "cosine", n_neighbors = 30L, min_dist = 0.3
    )
  colnames(umap) <- paste0(i, "_UMAP_", 1:ncol(umap))
  key <- paste0(i, "_umap")
  ds[[key]] <- Seurat::CreateDimReducObject(
    embeddings = umap, assay = "RNA", misc = list(from = i)
    )
}

# Plot
Seurat::DimPlot(ds, group.by = "condition", reduction = "N5000_pca40_umap") +
  Seurat::DimPlot(ds, group.by = "condition", reduction = "N5000_mnn40_umap")

fn <- paste0(dest, "umap", "_", "condition", ".", "png")
ggplot2::ggsave(fn, width = 10, height = 5)

ce <- "mnn40_umap"
gene <- "CD3D"
Seurat::FeaturePlot(ds, gene, reduction = ce) +
  viridis::scale_color_viridis(option = "A", direction = -1) +
  ggplot2::coord_fixed()

# Clustering -------------------------------------------------------------------
ds <- Seurat::FindNeighbors(ds, reduction = ce, dims = 1:ncols)
ds <- Seurat::FindClusters(ds)

# Exploration ------------------------------------------------------------------

# Select embedding
ce <- "mnn50_umap"

# Plot
Seurat::FeaturePlot(ds, "percent.mito", pt.size = 0.1, reduction = ce) +
  viridis::scale_color_viridis(option = "A", direction = -1) +
  ggplot2::coord_fixed() +
  ggplot2::facet_wrap(~ds@meta.data$condition) +
  ggplot2::theme_classic(20) +
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = "grey94")) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(barheight = 15, ticks = FALSE)
  )

fn <- paste0(dest, ce, "_", "pMt", "-", "condition", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 5.5)

# Plot
gene <- "CD163"
Seurat::FeaturePlot(ds, gene, pt.size = 0.1, reduction = ce) +
  viridis::scale_color_viridis(option = "A", direction = -1) +
  ggplot2::coord_fixed() +
  ggplot2::facet_wrap(~ds@meta.data$condition) +
  ggplot2::theme_classic(20) +
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = "grey93")) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(barheight = 15, ticks = FALSE)
  )

fn <- paste0(dest, "genes/", "umap", "_", gene, "-", "condition", ".", "png")
dir.create(dirname(fn), showWarnings = FALSE)
ggplot2::ggsave(fn, width = 12, height = 5.5)

# Differential expression ------------------------------------------------------

# Fine-tune cluster resolution
ds <- Seurat::FindClusters(ds, algorithm = 3, resolution = 0.9)
Seurat::DimPlot(ds, label = TRUE, label.box = TRUE, repel = TRUE) +
  ggplot2::coord_fixed() + ggplot2::theme_void(20) +
  ggplot2::labs(title = "Cluster")

fn <- paste0(dest, "umap", "_", "cluster", ".", "png")
ggplot2::ggsave(fn, width = 5.65, height = 6)

# Calculate differential expression
markers <- Seurat::FindAllMarkers(ds, only.pos = TRUE)

# Plot
genes <- unlist(lapply(split(markers$gene, markers$cluster), head))
Seurat::DoHeatmap(ds, features = genes)

fn <- paste0(dest, "heatmap", "_", "cluster", ".", "png")
ggplot2::ggsave(fn, width = 7, height = 9)
