#!user/bin/env R

# Author: Sean Maden
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
  # scan for nested metadata and expand into new appended variables
  for(c in 1:nrow(linefile)){
    linec = as.character(linefile[c,])
    if(grepl(".*:.*",linec)){
      linecsplitlist = unlist(strsplit(linec,":"))
      if(length(linecsplitlist)==2){
        linecsplitlist[1] = gsub("^.*= ","",linecsplitlist[1])
        linecsplitlist[1] = gsub("^ | $","",linecsplitlist[1])
        linecsplitlist[1] = gsub(" ","_",linecsplitlist[1])
        linecsplitlist[2] = gsub("^ | $","",linecsplitlist[2])
        linecsplitlist[2] = gsub(" ","_",linecsplitlist[2])
        if(!(grepl("^http$|^https$|^ftp$",linecsplitlist[1])) & 
           nchar(as.character(linecsplitlist[2]))<20){
          newcval = as.vector(gsub(".*:","",linecsplitlist[2]))
          names(newcval) <- linecsplitlist[1]
          dffile <- cbind(dffile,t(as.matrix(newcval)))
        }
      }
    }
  }
  jsoni <- toJSON(dffile, pretty=T)
  write(jsoni, file = file.path(gsm_json_destdir,paste0(gsmsoft_fn,".json")))
}

suppressWarnings(soft2json(gsmsoft_fn = commandArgs(T)[1],
          gsm_softpath = commandArgs(T)[2],
          gsm_json_destdir = commandArgs(T)[3]))


