# Analysis of the complete immune cell data

# Global options ---------------------------------------------------------------
file <- "data/complete.Rds"
dest <- paste0("analysis/complete/analysis/")
dir.create(dest, recursive = TRUE)

# Read data --------------------------------------------------------------------
ds <- readRDS(file)

# Normalize --------------------------------------------------------------------

ds <- Seurat::NormalizeData(ds)
ds <- Seurat::FindVariableFeatures(ds, nfeatures = 5000)

# Dimensional reduction  -------------------------------------------------------

pca <- irlba::irlba(ds@assays$RNA@data[ds@assays$RNA@var.features, ], 50)
pca$u

ds[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = matrix(
    data = pca$v, ncol = ncol(pca$v),
    dimnames = list(colnames(ds), paste0("PC_", 1:ncol(pca$u)))
  ),
  loadings = matrix(
    data = pca$u, ncol = ncol(pca$u),
    dimnames = list(ds@assays$RNA@var.features, paste0("PC_", 1:ncol(pca$u)))
  )
)

ds <- Seurat::RunUMAP(ds, dims = 1:40)

Seurat::DimPlot(ds, group.by = "timepoint") + ggplot2::coord_fixed()

# Clustering -------------------------------------------------------------------

# Seurat clustering
ds <- Seurat::FindNeighbors(ds, )
ds <- Seurat::FindClusters(ds, resolution = 1)

Seurat::DimPlot(ds, group.by = "seurat_clusters") + ggplot2::coord_fixed()

# Save file --------------------------------------------------------------------
saveRDS(ds, file)
