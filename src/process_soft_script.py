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
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
from utilities import monitor_processes
import settings
settings.init()

json_flist=[]
fnpatt='.*json.filt$'
gsmjsonpath = os.path.join('recount-methylation-files','gsm_json_filt')
timestamp=gettime_ntp()
nprocsamp=500
nmaxproc=20
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





if len(ll)>1:
    for loc, sublist in enumerate(ll, 1):
        # each loc in ll represents a new screen index, check vs. screen max
        if loc <= nmaxproc:
            cmdlist0 = ['screen','-S',"msrapsession"+str(loc), '-dm', 
                'python3', psoftpath, '--msraplist', 
                ' '.join(str(item) for item in ll[0]),'--ntptime', ts,
                '--gsm_jsonpath', gsmjsonpath
            ]
            proc = subprocess.Popen(cmdlist0, stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE)
            process_list.append(proc)
else:
    cmdlist0 = ['screen','-S',"msrapsession"+str(loc), '-dm', 
                'python3', psoftpath, '--msraplist', 
                ' '.join(str(item) for item in ll[0]),'--ntptime', ts,
                '--gsm_jsonpath', gsmjsonpath
            ]
    proc = subprocess.Popen(cmdlist0, stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE)
    process_list.append(proc)