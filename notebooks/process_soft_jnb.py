#!/usr/bin/env python3

"""

Authors: Sean Maden, Abhi Nellore

To process the experiment soft files, Recount Methylation begins by downloading
the soft metadata file for each experiment. Each downloaded soft file thus has a 
corresponding GSE experiment ID. 

Experiment soft file downloads are managed along with GSM/sample idat downloads
from the celery job queue manager. GSE soft files are downloaded first, followed
by the valid GSM sample idats. Files are downloaded to a temp dir then validated 
against existing files in the 'gse_soft' subdir of 'recount-methylation-files'. 
After successful download, followed by validation showing a new soft file or 
that an existing soft file was updated, a call to add a new GSE soft file record 
to RMDB is made. This record includes the relevant file and download metadata.

The filenames of experiment soft files, like other Recount Methylation files,
are versioned using NTP timestamps in the filename, along with the experiment id
and the name of the original file that was downloaded from the GEO FTP call. 

"""

import subprocess
import glob
import sys
import os
import re
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import edirect_query
from edirect_query import gsm_query, gse_query, gsequery_filter  
from utilities import gettime_ntp, getlatest_filepath, querydict
from utilities import get_queryfilt_dict  
import settings
settings.init()

import process_soft

"""

To get from the GSE soft file to the mapped GSM metadata, several steps are 
taken: 1. expand GSE soft files; 2. extract valid GSM subsections from GSE soft
files; 3. convert extracted GSM soft data to JSON format; 4. pass the valid 
JSON-formatted GSM soft metadata to MetaSRA-pipeline , using the specified
fork of this application.

"""

"""

First, expand GSE soft files. This unpacks downloaded experiment soft files with
the option of removing compressed files after expansion:

"""

process_soft.expand_soft()

"""

Then extract GSM sample metadata. This will create new sample files in the 
gsm_soft files directory:

"""

process_soft.extract_gsm_soft(gse_softpath = settings.gsesoftpath, 
    gsm_softpath = settings.gsmsoftpath, gsmsoft_destpath = settings.gsmsoftpath)

extract_gsm_soft(gse_softpath = settings.gsesoftpath, 
    gsm_softpath = settings.gsmsoftpath, 
    gsmsoft_destpath = settings.gsmsoftpath)


# if restart necessary before extract_gsm_soft completes, do the following
# copy files from temp dir
eqfiltdict = get_queryfilt_dict()
gsmsoftlist = os.listdir("recount-methylation-files/gsm_soft") # get processed gsm ids
gsmidfilt = [fn.split(".")[1] for fn in gsmsoftlist]
gseid_processnew = []
vall = [v for v in eqfiltdict.values()]
keyl = [k for k in eqfiltdict.keys()]
for i, val in enumerate(vall):
    valf = [v for v in val if not v in gsmidfilt]
    if len(valf) > 0:
        gseid_processnew.append(keyl[i])
    print(str(i))

# get corresponding filenames for filtered study IDs
gsesoftl = os.listdir("recount-methylation-files/gse_soft")
gsesoftl = [fn for fn in gsesoftl if not "gz" in fn]
len(gsesoftl)
gsesoftlid = [gse.split(".")[0] for gse in gsesoftl]
gsefn_processnew = []
for i, gseid in enumerate(gsesoftlid):
    if gseid in gseid_processnew:
        gsefn_processnew.append(gsesoftl[i])

# rerun the function
extract_gsm_soft(gsesoft_flist = gsefn_processnew, gse_softpath = settings.gsesoftpath, 
    gsm_softpath = settings.gsmsoftpath, gsmsoft_destpath = settings.gsmsoftpath)

# from bash, move GSM SOFT files from temp dir
# mv -v recount-methylation-epic/recount-methylation-files/temp/tmp85udqr4o/* recount-methylation-epic/recount-methylation-files/gsm_soft/


"""

Next, convert GSM sample metadata from the soft file format to valid JSON format. 
This function passes sample soft files with an R function to convert to JSON.

"""
gsm_jsonpath = os.path.join("recount-methylation-files", "gsm_json")
gsm_softpath = os.path.join("recount-methylation-files", "gsm_soft")

process_soft.gsm_soft2json(gsm_jsonpath = gsm_jsonpath, gsm_softpath = gsm_softpath)

gsm_soft2json(gsm_jsonpath = gsm_jsonpath, gsm_softpath = gsm_softpath)



# run the filter script
#!/usr/bin/env R

# Author: Sean Maden
# This script shows how JSON-formatted sample/GSM metadata, 
# extracted from GSE SOFT files, was filtered prior to running in
# MetaSRA-pipeline. Most likely/highest-confidence sample type 
# predictions from this run were stored under the "sampletype"
# column in the main metadata table (Table S1).

library(readr)
library(jsonlite)

#------------------------
# assign dir's and path's
#------------------------
# input dir/path
json.dn <- "gsm_json"
readpath <- file.path("recount-methylation-files", json.dn)

# output dir/path
jfilt.dn <- paste0(json.dn, "_filt") # json.dn # paste0(json.dn, "_filt")
destpath <- file.path("recount-methylation-files", jfilt.dn)
if(!dir.exists(destpath)){dir.create(destpath)}

# new files stem
filestem <- "json.filt"

# keys of interest
keys.list <- c("!Sample_characteristics_ch1", "!Sample_source_name_ch1", 
               "!Sample_title")

#------------------
# filter json files
#------------------
# get unfiltered soft files
lf.json <- list.files(readpath)
lf.json <- lf.json[!grepl(".*filt.*", lf.json)]
# write filtered data as new json files
for(i in seq(length(lf.json))){
  fni <- lf.json[i]
  ts <- unlist(strsplit(fni,"\\."))[1] # timestamp
  gsmi <- unlist(strsplit(fni,"\\."))[2] # gsm id
  writefn <- paste(ts, gsmi, filestem, sep = ".")
  writepath <- file.path(destpath, writefn)
  rjsoni <- jsonlite::fromJSON(file.path(readpath,fni))
  # filter on valid sample-specific keys
  message("filtering keys for file ",i)
  rjsoni.keys <- colnames(rjsoni); rf <- list()
  for(k in keys.list){rekf <- rjsoni.keys[grepl(k, rjsoni.keys)]
    for(f in rekf){rf[[f]] <- as.character(unlist(rjsoni[f]))}
  }
  # write formatted json with top and bottom outside brackets
  message("writing filtered json data for file ",i)
  jsoni <- jsonlite::toJSON(rf, pretty=T, auto_unbox = T)
  write_lines("[", writepath); write_lines(jsoni, writepath, append=T)
  write_lines("]", writepath, append=T)
  message("finished file ",i)
}


"""

Finally map the JSON-formatted GSM metadata using MetaSRA-pipeline. 

Note, because the mapping process can take quite a while, depending on file size 
and format, we highly recommend deploying the pipeline across several background 
or detached screen sessions using the wrapper function. 

You can experiment with the number of screens and samples per screen that works 
best on your system. 

"""

gsm_jsonpath = os.path.join('recount-methylation-files','gsm_json_filt')
process_soft.msrap_launchproc(nprocsamp=100, nmaxproc=5, timelim=1000000, 
    fnpatt = '.*json.filt$', gsmjsonpath = gsm_jsonpath)



msrap_launchproc(nprocsamp=50, nmaxproc=20, timelim=2800, 
    fnpatt = '.*json\\.filt$', gsmjsonpath = jsonfilt_path)

