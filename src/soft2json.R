#!user/bin/env R

# soft2json.R
# Script to convert GSM soft file to JSON format
# Arguments
# 1. gsmsoft_fn (str) : GSM SOFT filename. 
# 2. gsm_softpath (str) : Filepath to search for gsmsoft_fn().
# 3. gsm_json_destdir (str) : Destination directory to store new JSON file.

require(jsonlite)

soft2json <- function(gsmsoft_fn, gsm_softpath, gsm_json_destdir){
  linefile <- read.table(file.path(gsm_softpath, gsmsoft_fn),sep="\n")
  dffile <- as.data.frame(matrix(nrow=1,ncol=nrow(linefile)))
  colnames(dffile) <- gsub(" =.*","",linefile[,1])
  dffile[1,] <- gsub("!.*= ","",linefile[,1])
  jsoni <- toJSON(dffile, pretty=T)
  write(jsoni, file = file.path(gsm_json_destdir,paste0(gsmsoft_fn,".json")))
}

suppressWarnings(soft2json(gsmsoft_fn = commandArgs(T)[1],
          gsm_softpath = commandArgs(T)[2],
          gsm_json_destdir = commandArgs(T)[3]))


