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

extract_gsm_soft(gse_softpath = settings.gsesoftpath, gsm_softpath = settings.gsmsoftpath, 
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
    if len(valf) > 10:
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


# OLD METHOD

gsm_jsonpath = os.path.join('recount-methylation-files','gsm_json_filt')
process_soft.msrap_launchproc(nprocsamp=100, nmaxproc=5, timelim=1000000, 
    fnpatt = '.*json.filt$', gsmjsonpath = gsm_jsonpath)

msrap_launchproc(nprocsamp=50, nmaxproc=20, timelim=2800, 
    fnpatt = '.*json\\.filt$', gsmjsonpath = jsonfilt_path)


# NEW METHOD

cd recount-methylation-hm450k/
python3


import os
import sys
import re
import gzip
import shutil
import subprocess
import filecmp
import tempfile
import pickle
from datetime import datetime
import time
from random import shuffle
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
from utilities import monitor_processes
import settings
settings.init()

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



def msrap_getsamples(json_flist=[], fnpatt=".*json.filt$", 
    gsmjsonpath=os.path.join("recount-methylation-files", "gsm_json_filt"), 
    nprocsamp=50, nmaxproc=20):
    """ msrap_getsamples
        
        Get the validated samples file list

        Arguments:
            * json_flist (list) : List of GSM JSON filenames to process. If not 
                provided, function automatically detects any new GSM JSON files
                without available MetaSRA-pipeline outfiles.
            * fnpatt (str): Filename pattern of valid json files to identify.
            * gsmjsonpath (path): Path to JSON formatted sample SOFT data.
            * nprocsamp (int) : Number of samples to process per screen deployed.
            * nmaxproc (int) : Maximum processes to launch
            * timelim (int) : time limit (minutes) for monitoring processes.
            * statint (int) : time (seconds) to sleep before next status update.
        Returns:
            (Null) Generates >=1 processes for file sublists

    """
    print("Checking dirs for msrapout and msrap logs...")
    os.makedirs(settings.gsmmsrapoutpath, exist_ok=True)
    os.makedirs(settings.msraplogspath, exist_ok=True)
    # detect gsm soft files
    psoftpath = settings.psoftscriptpath
    if os.path.exists(psoftpath):
        print("Process soft script found at: "+str(psoftpath))
    gsmsoftpath = settings.gsmsoftpath
    gsmmsrapoutpath = settings.gsmmsrapoutpath
    jsonfnpattern = fnpatt
    rjson = re.compile(jsonfnpattern)
    msrapoutfnpattern = settings.msrapoutfnpattern
    rmsrapout = re.compile(msrapoutfnpattern)
    # generate fl list of valid json files that haven't been processed yet
    fl = []
    if json_flist and len(json_flist)>0:
        jsonfnlist = list(filter(rjson.match, json_flist)) 
        jsongsmlist = [x.split('.')[1] for x in jsonfnlist]
    else:
        json_flist = os.listdir(gsmjsonpath)
        jsonfnlist = list(filter(rjson.match, json_flist)) 
        jsongsmlist = [x.split('.')[1] for x in jsonfnlist]
    msrapoutfnlist = os.listdir(gsmmsrapoutpath) 
    msrapoutfnlist = list(filter(rmsrapout.match, msrapoutfnlist))
    print("Found "+str(len(msrapoutfnlist))+" files with pattern "
        +msrapoutfnpattern+". Continuing...")
    msrapgsmlist = [x.split('.')[2] for x in msrapoutfnlist]
    gsmprocess = [g for g in jsongsmlist 
            if not g in msrapgsmlist and
            g[0:3]=='GSM'
        ]
    for index, gsmid in enumerate(gsmprocess):
        gjsonfpath = getlatest_filepath(filepath=gsmjsonpath,
                filestr=gsmid, embeddedpattern=True, tslocindex=0,
                returntype='returnlist'
            )
        if gjsonfpath and len(gjsonfpath)==1:
            gjsonfn = [os.path.basename(gjsonfpath[0])]
        else:
            gjsonfn = [os.path.basename(fn) for fn in gjsonfpath]
        gjsonfn = gjsonfn[0]
        fl.append(gjsonfn)
        numi = 100*(index/len(gsmprocess))
        perci = str(round(numi,2))
        print("Appended file "+gjsonfn+" to files list to process. "
            +"Progress: "+str(index)+"/"+str(len(gsmprocess))+"="
            +perci+"%. Continuing...")
    # form list of fn lists based on nscreensi and indices/slices
    if fl:
        print("Forming list of fn lists for screen deployment...")
        ll = []
        rangelist = [i for i in range(0, len(fl), nprocsamp)]
        for enum, i in enumerate(rangelist[:-1]):
            ll.append(fl[i:rangelist[enum+1]])
        if len(fl[rangelist[-1]::]) > 0:
            ll.append(fl[rangelist[-1]::])
    else:
        print("Error, no files list object to process. Returning...")
        return None
    print('screens ll list, len = ' + str(len(ll)))
    print('nmax screens = '+str(nmaxproc))
    return ll







ll = msrap_getsamples(nprocsamp=2100, nmaxproc=22)





lall = ll
for ll in lall:
    json_flist=ll
    jsonpatt=".*json.filt$"
    gsm_jsonpath = os.path.join("recount-methylation-files", "gsm_json_filt")
    timestamp=gettime_ntp()
    print("Validating files on GSM equery results...")
    eqfiltdict=get_queryfilt_dict()
    validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ]
    msrap_runpath = settings.msraprunscriptpath
    # if filenames provided, form list, else list json dir contents
    print("Validating JSON files on filename patterns...")
    if json_flist and len(json_flist)>0:
        rjson = re.compile(jsonpatt)
        gsm_json_fn_list = list(filter(rjson.match, json_flist)) 
    else:
        gsm_json_fn_list = os.listdir(gsm_jsonpath)
        rjson = re.compile(jsonpatt)
        gsm_json_fn_list = list(filter(rjson.match, gsm_json_fn_list))
    # make write path
    msrap_destpath = settings.gsmmsrapoutpath
    os.makedirs(msrap_destpath, exist_ok=True)
    print("Getting read and write paths...")
    # define file read paths
    gsm_read_pathl = [os.path.join(gsm_jsonpath, fn) for fn in gsm_json_fn_list]
    # define write paths, same file order as read paths
    msrap_fn = settings.msrapfnstem
    ts_gsm_fnl = [".".join(fn.split(".")[0:2]) for fn in gsm_json_fn_list]
    gsm_write_fnl = [".".join([timestamp,outfn,msrap_fn]) for outfn in ts_gsm_fnl]
    gsm_write_pathl = [os.path.join(msrap_destpath,fn) for fn in gsm_write_fnl]
    # note -- pass 1 long args list once
    args0 = ";".join(gsm_read_pathl)
    args1 = ";".join(gsm_write_pathl)
    cmdlist = " ".join(["python2", msrap_runpath, "--fnvread", '"'+args0+'"', "--fnvwrite", '"'+args1+'"'])
    print("Calling the subprocess...")
    proc = subprocess.call(cmdlist, shell=True)
    print("Finished calling subprocess. Returning...") 



#-----
# old
#-----



import os
import sys
import re
import gzip
import shutil
import subprocess
import filecmp
import tempfile
import pickle
from datetime import datetime
import time
from random import shuffle
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
from utilities import monitor_processes
import settings
settings.init()

json_flist=[]
fnpatt='.*json.filt$'
gsmjsonpath = os.path.join('recount-methylation-files','gsm_json_filt')
timestamp=gettime_ntp()
nprocsamp=500
nmaxproc=15
timelim=2800
statint=5

if not os.path.exists(settings.msraprunscriptpath):
        print("Error: MetaSRA-pipeline script not found. Please check your "
            +"local download of MetaSRA-pipeline.")
        return None
print("Checking dirs for msrapout and msrap logs...")
os.makedirs(settings.gsmmsrapoutpath, exist_ok=True)
os.makedirs(settings.msraplogspath, exist_ok=True)
# detect gsm soft files
psoftpath = settings.psoftscriptpath
if os.path.exists(psoftpath):
    print("Process soft script found at: "+str(psoftpath))
gsmsoftpath = settings.gsmsoftpath
gsmmsrapoutpath = settings.gsmmsrapoutpath
jsonfnpattern = fnpatt
rjson = re.compile(jsonfnpattern)
msrapoutfnpattern = settings.msrapoutfnpattern
rmsrapout = re.compile(msrapoutfnpattern)
# generate fl list of valid json files that haven't been processed yet
fl = []
if json_flist and len(json_flist)>0:
    jsonfnlist = list(filter(rjson.match, json_flist)) 
    jsongsmlist = [x.split('.')[1] for x in jsonfnlist]
else:
    json_flist = os.listdir(gsmjsonpath)
    jsonfnlist = list(filter(rjson.match, json_flist)) 
    jsongsmlist = [x.split('.')[1] for x in jsonfnlist]
msrapoutfnlist = os.listdir(gsmmsrapoutpath) 
msrapoutfnlist = list(filter(rmsrapout.match, msrapoutfnlist))
print("Found "+str(len(msrapoutfnlist))+" files with pattern "
    +msrapoutfnpattern+". Continuing...")
msrapgsmlist = [x.split('.')[2] for x in msrapoutfnlist]
gsmprocess = [g for g in jsongsmlist 
        if not g in msrapgsmlist and
        g[0:3]=='GSM'
    ]
for index, gsmid in enumerate(gsmprocess):
    gjsonfpath = getlatest_filepath(filepath=gsmjsonpath,
            filestr=gsmid, embeddedpattern=True, tslocindex=0,
            returntype='returnlist'
        )
    if gjsonfpath and len(gjsonfpath)==1:
        gjsonfn = [os.path.basename(gjsonfpath[0])]
    else:
        gjsonfn = [os.path.basename(fn) for fn in gjsonfpath]
    gjsonfn = gjsonfn[0]
    fl.append(gjsonfn)
    numi = 100*(index/len(gsmprocess))
    perci = str(round(numi,2))
    print("Appended file "+gjsonfn+" to files list to process. "
        +"Progress: "+str(index)+"/"+str(len(gsmprocess))+"="
        +perci+"%. Continuing...")
# form list of fn lists based on nscreensi and indices/slices
if fl:
    print("Forming list of fn lists for screen deployment...")
    ll = []
    rangelist = [i for i in range(0, len(fl), nprocsamp)]
    for enum, i in enumerate(rangelist[:-1]):
        ll.append(fl[i:rangelist[enum+1]])
    if len(fl[rangelist[-1]::]) > 0:
        ll.append(fl[rangelist[-1]::])
else:
    print("Error, no files list object to process. Returning...")
    return None
print('screens ll list, len = ' + str(len(ll)))
print('nmax screens = '+str(nmaxproc))

ll = ll[0:nmaxproc] # slice ll based on screen count max
ts = timestamp # single timestamp call shared across screens

# parse a single session
lall = ll
for ll in lall:
    json_flist=ll
    jsonpatt=".*json.filt$"
    gsm_jsonpath = os.path.join("recount-methylation-files", "gsm_json_filt")
    timestamp=gettime_ntp()
    print("Validating files on GSM equery results...")
    eqfiltdict=get_queryfilt_dict()
    validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ]
    msrap_runpath = settings.msraprunscriptpath
    # if filenames provided, form list, else list json dir contents
    print("Validating JSON files on filename patterns...")
    if json_flist and len(json_flist)>0:
        rjson = re.compile(jsonpatt)
        gsm_json_fn_list = list(filter(rjson.match, json_flist)) 
    else:
        gsm_json_fn_list = os.listdir(gsm_jsonpath)
        rjson = re.compile(jsonpatt)
        gsm_json_fn_list = list(filter(rjson.match, gsm_json_fn_list))
    # make write path
    msrap_destpath = settings.gsmmsrapoutpath
    os.makedirs(msrap_destpath, exist_ok=True)
    print("Getting read and write paths...")
    # define file read paths
    gsm_read_pathl = [os.path.join(gsm_jsonpath, fn) for fn in gsm_json_fn_list]
    # define write paths, same file order as read paths
    msrap_fn = settings.msrapfnstem
    ts_gsm_fnl = [".".join(fn.split(".")[0:2]) for fn in gsm_json_fn_list]
    gsm_write_fnl = [".".join([timestamp,outfn,msrap_fn]) for outfn in ts_gsm_fnl]
    gsm_write_pathl = [os.path.join(msrap_destpath,fn) for fn in gsm_write_fnl]
    # note -- pass 1 long args list once
    args0 = ";".join(gsm_read_pathl)
    args1 = ";".join(gsm_write_pathl)
    cmdlist = " ".join(["python2", msrap_runpath, "--fnvread", '"'+args0+'"', "--fnvwrite", '"'+args1+'"'])
    print("Calling the subprocess...")
    proc = subprocess.call(cmdlist, shell=True)
    print("Finished calling subprocess. Returning...")    



