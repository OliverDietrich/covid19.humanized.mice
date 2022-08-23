# Exploration of the complete immune cell data

# Global options ---------------------------------------------------------------
file <- "data/complete.Rds"
dest <- paste0("analysis/complete/", Sys.Date(), "/")
dir.create(dest, recursive = TRUE)

# Read data --------------------------------------------------------------------
ds <- readRDS(file)

# Calculate differential expression --------------------------------------------

cells <- names(which(ds$Mo_Mac))

Seurat::Idents(ds) <- ds$seurat_clusters
markers <- Seurat::FindAllMarkers(
  subset(ds, cells = cells), logfc.threshold = 0.5, only.pos = TRUE
  )

markers$logP <- -log10(markers$p_val_adj)
markers$logP[markers$logP == Inf] <- 300

markers$mean <- sparseMatrixStats::rowMeans2(ds@assays$RNA@data[markers$gene, ])

p <- ggplot2::ggplot(
  markers, ggplot2::aes(avg_log2FC, logP, label = gene, col = cluster)
  ) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~cluster) +
  ggplot2::theme_classic(10)

plotly::ggplotly(p)

fn <- paste0(dest, "cluster-markers", ".html")
htmlwidgets::saveWidget(plotly::ggplotly(p), fn)

# Explore gene expression ------------------------------------------------------

# Timepoint & Macrophages
plot_embedding(ds, "timepoint") + plot_embedding(ds, "Mo_Mac") +
  plot_embedding(ds, "seurat_clusters", label = "label")
fn <- paste0(dest, "umap_timepoint-mono-cluster", ".png")
ggplot2::ggsave(fn, width = 12, height = 4.8, bg = "white")

# Timepoint & genes
genes <- c("TREM2", "IL22")
plot_embedding(ds, "timepoint") + plot_embedding(ds, genes[1]) +
  plot_embedding(ds, genes[2])
fn <- paste0(dest, "umap_timepoint-", stringr::str_flatten(genes, "-"), ".png")
ggplot2::ggsave(fn, width = 12, height = 5, bg = "white")

# Enrichment of COVID-19 gene sets ---------------------------------------------

# Retrieve gene set dictionary
dict <- tempfile()
download.file(
  "https://nubes.helmholtz-berlin.de/s/25ZGCXAHdBffZCP/download", dict
)
dict <- readxl::read_excel(dict)
names(dict) <- c("ref", "term", "disease", "gene")

unique(dict$term)
ds$geneset <- Seurat::AddModuleScore(
  ds, list(dict$gene[dict$term == "IPFeMp"])
)$Cluster1

plot_embedding(ds, "geneset")
