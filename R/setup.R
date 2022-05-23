# Analysis of COVID-19 in a humanized mouse model

# Global options ---------------------------------------------------------------
dir <- "data/raw/GSE200562/"

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
    Mo_Mac    = colnames(ds[[i]]) %in% read.table(file)[, 1]
  )
}

# Combine datasets
lapply(ds, rownames)

# ds <- do.call(cbind, ds)

coldata <- do.call(rbind, coldata)
