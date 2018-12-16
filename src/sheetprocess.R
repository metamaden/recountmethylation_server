library(minfi)

# sheetprocess.R
# Script to process data sheets for Recount Methylation.
# Notes and TODO:
# fix basename files, one row per sample, one basename per row (no grn/red idat info)
# fix read-in so colnames recognized properly
# fix colnames, use 'Basename' instead of idat_fn

sheet_fn = list.files()
sheetidat <- read.table(sheet_fn[grep(".*rsheet.idats$", sheet_fn)], 
                        sep = " ", stringsAsFactors = F)
colnames(sheetidat) <- as.character(as.matrix(sheetidat)[1,])
sheetidat <- sheetidat[2:nrow(sheetidat),]
colnames(sheetidat)[6] <- "Basename"
sheetidat <- sheetidat[!duplicated(sheetidat$Basename),]
sheetidat <- sheetidat[!duplicated(sheetidat$gsmid),]


mexp = read.metharray.exp(base=".",targets=sheetidat[c(1:2),], force=TRUE)