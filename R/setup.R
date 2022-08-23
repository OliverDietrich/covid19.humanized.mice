# Setup of the Seurat object

# Global options ---------------------------------------------------------------
dir <- "data/raw/GSE200562/"
dest <- "data/combined.Rds"

# Read files -------------------------------------------------------------------
cont <- list.files(dir)

# Separate bulk from single-cell
bulk <- cont[stringr::str_detect(cont, "txt")]
sc <- cont[!stringr::str_detect(cont, "txt")]

# Read matrix and supply row (genes) and colnames (cells)
sc.sets <- unique(stringr::str_split_fixed(sc, "_", 3)[, 1])
ds <- list()
for (i in sc.sets) {
  # Matrix
  file <- sc[stringr::str_detect(sc, i) & stringr::str_detect(sc, "mtx")]
  file <- paste0(dir, file)
  ds[[i]] <- Matrix::readMM(file)
  # Genes
  file <- sc[stringr::str_detect(sc, i) & stringr::str_detect(sc, "genes")]
  file <- paste0(dir, file)
  rownames(ds[[i]]) <- read.table(file)[, 1]
  # Cell barcodes
  file <- sc[stringr::str_detect(sc, i) & stringr::str_detect(sc, "barcodes")]
  file <- file[!stringr::str_detect(file, "mono")]
  file <- paste0(dir, file)
  colnames(ds[[i]]) <- read.table(file)[, 1]
}

# Create meta/coldata, add mono/macro annotation
coldata <- list()
for (i in names(ds)) {
  # Mono/macro annotation
  file <- sc[stringr::str_detect(sc, i) & stringr::str_detect(sc, "barcodes")]
  file <- file[stringr::str_detect(file, "mono")]
  file <- paste0(dir, file)

  coldata[[i]] <- data.frame(
    row.names = colnames(ds[[i]]),
    barcode   = colnames(ds[[i]]),
    dataset   = i,
    condition = stringr::str_split(basename(file), "_", 3, simplify = T)[, 2],
    Mo_Mac    = colnames(ds[[i]]) %in% read.table(file)[, 1]
  )
}

# Combine data sets ------------------------------------------------------------

# Identify set of shared features
fset <- character()
for (i in lapply(ds, rownames)) {
  fset <- c(fset, i[!i %in% fset])
}

# Add missing features (as zeros)
for (i in names(ds)) {
  mf <- fset[which(!fset %in% rownames(ds[[i]]))]
  nc <- ncol(ds[[i]])
  mm <- as(matrix(
    data = 0, nrow = length(mf), ncol = nc,
    dimnames = list(mf, colnames(ds[[i]]))
    ), "dgCMatrix")
  ds[[i]] <- rbind(ds[[i]], mm)
  ds[[i]] <- ds[[i]][fset, ]
}

# Check dimensions
lapply(ds, dim)

# Combine matrices
ds <- do.call(cbind, ds)

# Combine metadata
coldata <- do.call(rbind, coldata)
coldata <- coldata[match(coldata$barcode, colnames(ds)), ]
row.names(coldata) <- colnames(ds)

# Create Seurat object ---------------------------------------------------------

ds <- Seurat::CreateSeuratObject(
  counts = ds, meta.data = coldata, project = "mmCovLung", assay = "RNA"
)

# Save object ------------------------------------------------------------------
saveRDS(ds, dest)
