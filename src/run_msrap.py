#!/usr/bin/env python3

""" process_soft.py

    Authors: Sean Maden, Abhi Nellore
    
    Run the MetaSRA-pipeline.
    
    Notes:

    Functions:
"""

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

def run_metasrapipeline(json_flist=[], jsonpatt=".*json.filt$", 
    gsm_jsonpath=settings.gsmjsonfiltpath, timestamp=gettime_ntp()):
    """ run_metasrapipeline
        
        Run MetaSRA-pipeline on GSM JSON files. It is highly recommended to 
        implement this with parallelization using msrap_screens, instead of 
        instantiating directly!
        
        Arguments:
            * json_flist (list, optional) : List of JSON filename(s) to process. 
                If not provided, automatically targets all JSON files at 
                gsm_jsondir.
            * timestamp (str) : NTP timestamp version for expanded files.       
        
        Returns:
            * msrap_statlist (list), Int (1) or error: Whether MetaSRA-pipeline 
                successfully ran, generating a new MetaSRA file as side effect. 

    """
eqfiltdict=get_queryfilt_dict()
validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
    for gsmid in gselist
]
msrap_runpath = settings.msraprunscriptpath
    # if filenames provided, form list, else list json dir contents
if json_flist and len(json_flist)>0:
    rjson = re.compile(jsonpatt)
    gsm_json_fn_list = list(filter(rjson.match, json_flist)) 
else:
    gsm_json_fn_list = os.listdir(gsm_jsonpath)
    rjson = re.compile(jsonpatt)
    gsm_json_fn_list = list(filter(rjson.match, gsm_json_fn_list))

msrap_destpath = settings.gsmmsrapoutpath
os.makedirs(msrap_destpath, exist_ok=True)
msrap_statlist = []
msrap_fn = settings.msrapfnstem
# iterate over valid filenames at gsm json dir
process_list = []
args_list = []

for gsm_json_fn in gsm_json_fn_list:
    
    gsmid = gsm_json_fn.split('.')[1]
    
    if gsmid in validgsmlist:
        
outfn = os.path.splitext(gsm_json_fn)[0] # fn without extension
gsmjson_readpath = os.path.join(gsm_jsonpath, gsm_json_fn)
gsm_msrapout_writepath = os.path.join(msrap_destpath,
    ".".join([timestamp,outfn,msrap_fn]))
cmdlist = ['python2',
    msrap_runpath,
    gsmjson_readpath,
    gsm_msrapout_writepath
    ]
proc = subprocess.call(cmdlist, shell=False)
process_list.append(proc)
args_list.append([gsmjson_readpath, gsm_msrapout_writepath])

    else:
        msrap_statlist.append(None)
        print("GSM id : "+gsmid+" is not a valid HM450k sample. "
            +"Continuing...")
    return msrap_statlist

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
    return ll
    

def msrap_launchproc(json_flist=[], fnpatt=settings.jsonfnpattern, 
    gsmjsonpath=settings.gsmjsonpath, timestamp=gettime_ntp(), nprocsamp=50, 
    nmaxproc=20, timelim=2800, statint=5):
    """ msrap_launchproc
        Preprocess subsets of GSM JSON files in MetaSRA-pipeline in background, 
        with process monitoring
        Notes:
            *If no GSM JSON files list supplied to 'json_flist', then a new list 
                of GSMs is generated for valid GSM JSON files that don't already 
                have msrapout files available.
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
    process_list = []
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
    print("Finished launching background processes. Check screens for ongoing "
        +"status reports.")
    # monitor processes
    print("Beginning process monitoring...")
    monitor_processes(process_list=process_list, logpath=settings.msraplogspath)
    print("Finished with process monitoring.")
    print("Returning...")
    return None

def run_metasrapipeline2(json_flist=[], jsonpatt=".*json.filt$", 
    gsm_jsonpath = os.path.join("recount-methylation-files", "gsm_json_filt"), 
    timestamp=gettime_ntp()):
    """ run_metasrapipeline2
        
        Designed to run with modified script "run_pipeline.py"

        Run MetaSRA-pipeline on GSM JSON files. It is highly recommended to 
        implement this with parallelization using msrap_screens, instead of 
        instantiating directly!
        Arguments:
            * json_flist (list, optional) : List of JSON filename(s) to process. 
                If not provided, automatically targets all JSON files at 
                gsm_jsondir.
            * timestamp (str) : NTP timestamp version for expanded files.       
        Returns:
            * msrap_statlist (list), Int (1) or error: Whether MetaSRA-pipeline 
                uccessfully ran, generating a new MetaSRA file as side effect. 
    """
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
    args0 = gsm_read_pathl.join(";")
    args1 = gsm_write_pathl.join(";")
    print("Calling the subprocess...")
    cmdlist = " ".join(["python2", msrap_runpath, "--fnvread", '"'+args0+'"', 
        "--fnvwrite", '"'+args1+'"'])
    proc = subprocess.call(cmdlist, shell=True)
    print("Finished calling subprocess. Returning...")
    return True

if __name__ == "__main__":
    # the following is called by msrap_screens()
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--msraplist", type=str, required=True,
        default=None, help='Files list to process with MetaSRA-pipeline.')
    parser.add_argument("--ntptime", type=str, required=True,
        default=gettime_ntp(), help='NTP timestamp, as a string.')
    parser.add_argument("--gsm_jsonpath", type=str, required=False,
        default=gettime_ntp(), help='GSM JSON file path')
    args = parser.parse_args()
    # parse filename strings into list 
    flmsrap = [file for file in args.msraplist.split(' ')]
    run_metasrapipeline2(json_flist=flmsrap, timestamp=args.ntptime)