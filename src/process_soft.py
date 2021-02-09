#!/usr/bin/env python3

""" process_soft.py

    Authors: Sean Maden, Abhi Nellore
    
    Preprocess GSE/GSM Soft files for recount methylation samples. Includes GSM 
    Soft metadata extraction, methods to convert to JSON, and methods to process 
    samples in MetaSRA-pipeline.
    
    Notes:
    * To avoid propagating invalid (e.g. non-HM450k) experiments or samples
        through the pipeline, an equeryfiltdict dictionary object is called, 
        with valid GSE ids as keys and valid/filtered GSM ids listed as values.
    
    Functions:
        * expand_soft : Expand/extract a compressed GSE Soft file from GEO.
        * extract_gsm_soft : Extract GSM-level sample metadata from GSE soft 
            file.
        * gsm_soft2json : Convert GSM Soft metadata file to JSON for passage to
            MetaSRA-pipeline.
        * msrap_prepare_json : Concatenate multiple JSON sample files for 
            passage to MetaSRA-pipeline
        * run_metasrapipeline : run MetaSRA-pipeline in Python2 with screen.
        * msrap_screens : deploy multiple screen sessions running MetaSRA-
            pipeline at varying starting indices in the files list.
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

def expand_soft(rmcompressed=False):
    """ expand_soft
        Expand compressed GSE soft files, including filter on valid GSE ids.
        Arguments:
            * rmcompressed (T/F,Bool.) : Whether to remove compressed soft files
                once they have been successfully expanded.
        Returns:
            * rsoftd (list) : List of filenames and statuses, produces expanded 
                soft files as side effect.
    """
    eqfiltdict=get_queryfilt_dict()
    validgselist = list(eqfiltdict.keys())
    gsesoft_fpath = settings.gsesoftpath
    gsesoft_fnlist = os.listdir(gsesoft_fpath)  
    gseidlist =  [fn.split('.')[0] for fn in gsesoft_fnlist]
    gsepatt = settings.gsepatt
    rgse = re.compile(gsepatt)
    gseidlist = list(filter(rgse.match, gseidlist)) # valid GSE IDs
    print(str(gseidlist))
    softpatt = settings.softallpatt
    comppatt = settings.compsoftpatt
    expsoftpatt = settings.expsoftpatt
    rsoft = re.compile(softpatt)
    rcompsoft = re.compile(comppatt)
    rexpsoft = re.compile(expsoftpatt)
    gsesoft_compflist = [] # filtered fn list of compressed soft files
    for gseid in gseidlist:
        if gseid in validgselist:
            gse_filtfn = []
            flatest = getlatest_filepath(
                    filepath=settings.gsesoftpath, 
                    filestr=str(gseid), tslocindex=1, returntype='returnlist'
                )
            if flatest:
                if len(flatest)==1:
                    gsesoft_latestfn = [os.path.basename(flatest[0])]
                else:
                    gsesoft_latestfn = [os.path.basename(fpath) for fpath in 
                        flatest
                    ]
            else:
                print("Error getting latest files list for gseid : "+gseid)
                break
            # valid expanded soft
            gse_expsoftfn = []
            gse_expsoftfn = list(filter(rexpsoft.match, gsesoft_latestfn))
            # if no latest expanded soft file, search for valid compsoft file
            if len(gse_expsoftfn)==0:
                # valid soft
                gse_filtfn = list(filter(rsoft.match, gsesoft_latestfn)) 
                # valid cmpsoft
                gse_filtfn = list(filter(rcompsoft.match, gse_filtfn))
                if gse_filtfn and len(gse_filtfn)==1:
                    gsesoft_compflist.append(gse_filtfn[0])
        else:
            print("GSE id :"+gseid+" is not a valid HM450k experiment. "
                +"Continuing...")    
    rsoftd = {} # soft file statues dictionary to be returned
    softl_compressed = gsesoft_compflist
    if len(softl_compressed)>0:
        for softcompfile in softl_compressed:
            rsoftd[softcompfile] = [] # instantiate new dict val as list 
            with gzip.open(os.path.join(gsesoft_fpath, softcompfile), 'rb') as f_in:
                with open(os.path.join(gsesoft_fpath, softcompfile[:-3]), 'wb') as f_out:
                    try:
                        shutil.copyfileobj(f_in, f_out)
                        rsoftd[softcompfile].append(True) # if success
                    except shutil.Error as se:
                        rsoftd[softcompfile].append(se) # if failure
        statusindices = [i for i, x in enumerate(rsoftd[softcompfile]) 
                if x == True
            ]
        rmsuccess = [softl_compressed[i] for i in statusindices]
        if rmcompressed and len(rmsuccess) > 0:
            for compfilename in rmsuccess:
                os.remove(os.path.join(gsesoft_fpath,compfilename))
                rsoftd[compfilename].append(True) # if comp. file removed
    else: 
        print("No valid compressed soft files found at specified gse_softpath. "
            +"Are all valid soft files already expanded?")
        return None
    return rsoftd

def extract_gsm_soft(gsesoft_flist=[], softopenindex='.*!Sample_title.*', 
    softcloseindex='.*!Sample_data_row_count.*', timestamp=gettime_ntp(), 
    gse_softpath = settings.gsesoftpath, gsm_softpath = settings.gsmsoftpath, 
    gsmsoft_destpath = settings.gsmsoftpath, rmtempdir = True, validate=True):
    """ extract_gsm_soft
        Extract GSM soft file sections from GSE soft files.
        Arguments: 
            * gsesoft_flist (list, optional) : List of gse soft files to process
            * softopenindex (str) : Index of label/tag to open entry, defaults
                to sample title section.
            * softcloseindex (str) : Index of label/tag to close entry, defaults 
                to close just before possible by-CpG methylation table. To 
                include possible methylation data table, change to 
                '!sample_table_end'.
            * timestamp (str) : NTP timestamp version for expanded files.
            * rmtempdir (Bool.) : Whether to remove temp directory.
            * validate (Bool.) : Validate extracted GSM files against files in 
                gsm_soft directory?
        Returns:
            * newfilesd (dictionary), or error (null), generates GSM soft files 
                as a side effect.
    """
    eqfiltdict=get_queryfilt_dict()
    validgsmlist = list(set([gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ]))
    print("length validgsmlist : "+str(len(validgsmlist)))
    rvalidsoft = re.compile(".*soft$") # identify expanded GSE soft files
    gsmsoft_temppath = settings.temppath
    os.makedirs(gsm_softpath, exist_ok=True)
    os.makedirs(gsmsoft_temppath, exist_ok=True)
    temp_dir_make = tempfile.mkdtemp(dir=gsmsoft_temppath)
    if not gsesoft_flist or len(gsesoft_flist)==0:
        gse_soft_dirlist = os.listdir(gse_softpath)
    else:
        gse_soft_dirlist = gsesoft_flist
        gse_soft_dirlist = [gsefile for gsefile in gse_soft_dirlist
            if os.path.exists(os.path.join(gse_softpath, gsefile))
        ]
    gse_softlist = list(filter(rvalidsoft.match, gse_soft_dirlist))
    shuffle(gse_softlist)
    newfilesd = {} # new files, status dictionary to return
    print("new tempdir for writing soft files : "+str(temp_dir_make))
    print("length gse_softlist: "+str(len(gse_softlist)))
    rxopen = re.compile(softopenindex)
    rxclose = re.compile(softcloseindex)
    rxgsm = re.compile('GSM[0-9]*')
    rxgsmfile = re.compile('.*GSM.*')
    for gse_softfile in gse_softlist:
        print("Beginning gse softfile : "+gse_softfile)
        newfilesd[gse_softfile] = []
        openindex = []
        closeindex = []
        lsoft = []
        gse_softfile_path = os.path.join(gse_softpath, gse_softfile)
        with open(gse_softfile_path) as file:
            for num, line in enumerate(file, 0):
                if rxclose.search(line):
                    closeindex.append(num)
                if rxopen.search(line):
                    openindex.append(num)
                lsoft.append(line)
        print("for gse, found n = "+str(len(openindex))+" valid open "
                +"indices..")
        print("for gse, found n = "+str(len(lsoft))+" soft file lines. "
            +"Continuing...")
        for num, openi in enumerate(openindex,0):
            print("num : "+str(num))
            print("openi : "+str(openi))
            print("closei : "+str(closeindex[num]))
            try:
                gsm_softlines = lsoft[openi:closeindex[num]] # read gsm lines
            except:
                break
            gsmid_lines = [line for line in gsm_softlines
                if '!Sample_geo_accession' in line
            ]
            if len(gsmid_lines)==1:
                gsmid = str(rxgsm.findall(gsmid_lines[0])[0])
                print("Found GSM id : "+gsmid)
                gsm_softfn = ".".join([timestamp, gsmid, 'soft'])
                print("GSM id found : "+gsmid)
                print(gsmid+" in valid list... "+str(gsmid in validgsmlist))
                if gsmid in validgsmlist:
                    newfilesd[gse_softfile].append(gsm_softfn)
                    gsm_newfile_path = os.path.join(temp_dir_make, gsm_softfn)
                    write = [gsmfile.write(line) for line in open(gsm_newfile_path,"w+")]
                    open(gsm_newfile_path,"w+").write("\n".join(gsm_softlines))
                else: 
                    print("GSM id :"+gsmid+" is not a valid HM450k sample. "
                        +"Continuing...")
            else:
                print("GSM soft lines malformed! Continuing...")
    print("newfilesd : "+str(newfilesd))
    if validate:
        print("Beginning validation for files: ", list(newfilesd.keys()))
        for gse_softfn in list(newfilesd.keys()):
            gsmfilelist = list(filter(rxgsmfile.match, newfilesd[gse_softfn]))
            if gsmfilelist and len(gsmfilelist)>0:
                print(str(gsmfilelist))
                for gsmfile in gsmfilelist:
                    gsm_oldfile_path = ""
                    gsm_newfile_path = ""
                    gsm_softfn = gsmfile
                    gsmstr = gsm_softfn.split(".")[1]
                    print("gsmfile: "+str(gsmfile))
                    print("gsmstr : "+gsmstr)
                    gsm_newfile_path = os.path.join(temp_dir_make, gsm_softfn)
                    gsm_oldfile_path = getlatest_filepath(
                            filepath=gsmsoft_destpath, filestr=gsmstr, 
                            embeddedpattern=True, tslocindex=0
                        )
                    print("gsm_oldfile_path : "+str(gsm_oldfile_path))
                    print("gsm_newfile_path : "+str(gsm_newfile_path))
                    if os.path.exists(gsm_newfile_path):
                        if gsm_oldfile_path:
                            if filecmp.cmp(gsm_oldfile_path, gsm_newfile_path):
                                print("Identical GSM soft file detected, removing...")
                                os.remove(gsm_newfile_path)
                                newfilesd[gsmfile] = False
                            else:
                                print("New GSM soft file detected, moving from temp...")
                                shutil.move(gsm_newfile_path, os.path.join(
                                        gsmsoft_destpath, 
                                        os.path.basename(gsm_newfile_path))
                                    )
                                newfilesd[gsmfile] = True
                        else: 
                            print("New GSM soft file detected, moving from temp...")
                            shutil.move(gsm_newfile_path, os.path.join(
                                        gsmsoft_destpath, 
                                        os.path.basename(gsm_newfile_path))
                                    )
                            newfilesd[gsmfile] = True
                    else:
                        print("GSM soft file unavailable. Continuing...")
                        newfilesd[gsmfile] = False
    else:
        for fn in gsmfile:
            print("Moving file ", str(fn), "...")
            shutil.move(os.path.join(temp_dir_make, fn), 
                os.path.join(gsmsoft_destpath, fn))
    if rmtempdir:
        print("Removing tempdir...")
        shutil.rmtree(temp_dir_make) # remove tempdir
    return newfilesd 

def gsm_soft2json(gsm_softlist=[], scriptpath=settings.s2jscriptpath,
    gsm_jsonpath = settings.gsmjsonpath, gsm_softpath = settings.gsmsoftpath):
    """ gsm_soft2json
        Convert GSM soft file to JSON format Calls R script to coerce GSM soft 
        files (XML-like format) to valid JSON format.
        Arguments:
            * gsm_softlist (list, optional) : List of GSM soft filenames to 
                process.
            * scriptpath (str) : Path to R script for JSON conversion. If not 
                provided, automatically checks current working directory.
        Returns:
            * rlist object (list) of converted files and statuses, or error, 
                generates GSM JSON files as a side effect.
    """
    eqfiltdict = get_queryfilt_dict()
    validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ]
    if not os.path.exists(scriptpath):
        print("Error: Soft-to-JSON conversion script not detected at path. "
                +"Returning...")
        return None
    os.makedirs(gsm_jsonpath, exist_ok = True)
    rjsonlist_current = os.listdir(gsm_jsonpath) # current json dir contents
    # form list of gsm soft filenames, as gsm
    if gsm_softlist and len(gsm_softlist)>0:
        gsm_softfn_list = gsm_softlist
    else:
        gsm_softfn_list = os.listdir(gsm_softpath)
    # status dict, to return
    statd = {}
    if gsm_softfn_list and len(gsm_softfn_list)>0:
        # check existant json files before attempting conversion
        for gsmi, gsm_softfn in enumerate(gsm_softfn_list, 1):
            gsmid = gsm_softfn.split('.')[1]
            statd[gsm_softfn] = []
            if gsmid in validgsmlist:
                softts = gsm_softfn.split('.')[0] # soft file timestamp
                rgsm = re.compile('.*GSM.*')
                gsmid = gsm_softfn.split('.')[1]
                gsmid = str(rgsm.findall(gsmid)[0])
                gsmjson_latestfpath = getlatest_filepath(
                        filepath=gsm_jsonpath, filestr=gsmid, tslocindex=0,  
                        returntype='returnlist', embeddedpattern=True,
                    )
                if gsmjson_latestfpath and len(gsmjson_latestfpath)==1:
                    gsmjson_latestfpath = os.path.basename(
                            gsmjson_latestfpath[0]
                        )
                    jsonlatestts = gsmjson_latestfpath.split('.')[0] # json timestamp
                    if int(softts)>int(jsonlatestts):
                        try:
                            cmdlist = ['Rscript', scriptpath, gsm_softfn, 
                                gsm_softpath, gsm_jsonpath
                            ]
                            subprocess.call(cmdlist,shell=False)
                            statd[gsm_softfn].append(True)
                            print("R session launched for sample "+str(gsmi), 
                                end="\r")
                        except subprocess.CalledProcessError as e:
                            statd[gsm_softfn].append(None)
                            statd[gsm_softfn].append(e)
                else:
                    try:
                        cmdlist = ['Rscript', scriptpath, gsm_softfn, 
                            gsm_softpath, gsm_jsonpath
                        ]
                        subprocess.call(cmdlist,shell=False)
                        statd[gsm_softfn].append(True)
                        print("R session launched for sample "+str(gsmi), 
                            end="\r")
                    except subprocess.CalledProcessError as e:
                        statd[gsm_softfn].append(None)
                        statd[gsm_softfn].append(e)
            else:
                statd[gsm_softfn].append(None)
                print("Sample num "+str(gsmi)+" is not a valid HM450k sample. ",
                    end="\r")
    else:
        print("Error: No valid GSM Soft files to process from list. Returning.")
        return None
    # tally new json files generated
    rjsonlist_new = os.listdir(gsm_jsonpath)
    rjsonlist_return = [jfile for jfile in rjsonlist_new
        if not jfile in rjsonlist_current
    ]
    print("Generated "+str(len(rjsonlist_return))+" new json files in dir.")
    rlist = [statd, rjsonlist_return]
    return rlist

def run_metasrapipeline(json_flist=[], jsonpatt=".*json.filt$", 
    gsm_jsonpath = settings.gsmjsonpath, timestamp=gettime_ntp()):
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
                uccessfully ran, generating a new MetaSRA file as side effect. 
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
