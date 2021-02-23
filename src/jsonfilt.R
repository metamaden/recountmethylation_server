#!/usr/bin/env R

# Author: Sean Maden
# 
# This script shows how sample/GSM metadata JSON files were filtered prior
# to metadata mapping or running in MetaSRA-pipeline. This filtering step
# is intended to remove study-specific data (e.g. common to all samples) 
# while retaining key sample-specific data that can aid sample 
# identification (e.g. tags specific in the keys.list argument).

library(readr); library(jsonlite)

jsonfilt <- function(jsonext = ".*\\.json$", jsonfiltext = ".*\\.json\\.filt$", 
                     json.dn = "gsm_json", jfilt.dn = "gsm_json_filt", 
                     verbose = TRUE, files.dname = "recount-methylation-files",
                     keys.list = c("!Sample_characteristics_ch1", 
                                   "!Sample_source_name_ch1","!Sample_title")){
  readpath <- file.path(files.dname, json.dn)
  destpath <- file.path(files.dname, jfilt.dn)
  if(!dir.exists(readpath)){stop("Couldn't find read path at: ", readpath)}
  if(!dir.exists(destpath)){message("Making new destpath: ", destpath)
    dir.create(destpath)}
  lf.json <- list.files(readpath); lf.json <- lf.json[grepl(jsonext, lf.json)]
  if(verbose){message("Apply JSON filter for ",length(lf.json)," files...")}
  for(i in seq(length(lf.json))){
    fni <- lf.json[i]; ts <- unlist(strsplit(fni,"\\."))[1]
    gsmi <- unlist(strsplit(fni,"\\."))[2]
    writefn <- paste0(paste(ts, gsmi, sep = "."), jsonfiltext)
    writepath <- file.path(destpath, writefn)
    rjsoni <- jsonlite::fromJSON(file.path(readpath, fni))
    if(verbose){message("Filtering keys for file: ",i)}
    rjsoni.keys <- colnames(rjsoni); rf <- list()
    for(k in keys.list){
      rekf <- rjsoni.keys[grepl(k, rjsoni.keys)]
      for(f in rekf){
        rf[[f]] <- as.character(unlist(rjsoni[f]))}}
    if(verbose){message("Writing filtered JSON data for file ",i, "...")}
    jsoni <- jsonlite::toJSON(rf, pretty=T, auto_unbox = T)
    write_lines("[", writepath); write_lines(jsoni, writepath, append=T)
    write_lines("]", writepath, append=T);message("finished file ",i)
  };return(NULL)
}

jsonfilt()