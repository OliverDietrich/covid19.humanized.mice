# Download of the raw data from Gene Expresssion Omnibus (GEO)

# Global options ---------------------------------------------------------------
url <- "https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/download/?acc=GSE200562&format=file"
file <- "data/raw/GSE200562.tar"
dir.create(dirname(file), recursive = TRUE)

# Download files ---------------------------------------------------------------
download.file(url, file)
untar(file, exdir = stringr::str_remove(file, ".tar"))
file.remove(file)
