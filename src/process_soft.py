#!/usr/bin/env python3

""" process_soft.py

    Authors: Sean Maden, Abhi Nellore
    
    Preprocess GSE/GSM SOFT files for recount methylation samples. Includes GSM 
    SOFT metadata extraction, methods to convert to JSON, and methods to process 
    samples in MetaSRA-pipeline.
    
    Notes:
    * To avoid propagating invalid (e.g. non-HM450k) experiments or samples
        through the pipeline, an equeryfiltdict dictionary object is called, 
        with valid GSE ids as keys and valid/filtered GSM ids listed as values.
    
    Functions:
        * expand_soft : Expand/extract a compressed GSE SOFT file from GEO.
        * extract_gsm_soft : Extract GSM-level sample metadata from GSE SOFT 
            file.
        * gsm_soft2json : Convert GSM SOFT metadata file to JSON for passage to
            MetaSRA-pipeline.
        * msrap_prepare_json : Concatenate multiple JSON sample files for 
            passage to MetaSRA-pipeline
        * run_metasrapipeline : run MetaSRA-pipeline in Python2 with screen.
        * msrap_screens : deploy multiple screen sessions running MetaSRA-
            pipeline at varying starting indices in the files list.
"""

import os, sys, re, gzip, shutil, subprocess, filecmp, tempfile, pickle, time
from datetime import datetime; from random import shuffle
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
from utilities import monitor_processes; import settings; settings.init()

def expand_soft(rmcompressed=False):
    """ expand_soft
        
        Expand compressed GSE SOFT files, including filter on valid GSE ids.
        
        Arguments:
        * rmcompressed (T/F,Bool.) : Whether to remove compressed SOFT files
            once they have been successfully expanded.
        
        Returns:
        * rsoftd (list) : List of filenames and statuses, produces expanded 
            SOFT files as side effect.
                
    """
    eqfiltdict=get_queryfilt_dict(); validgselist = list(eqfiltdict.keys())
    gsesoft_fpath = settings.gsesoftpath
    gsesoft_fnlist = os.listdir(gsesoft_fpath)  
    gseidlist =  [fn.split('.')[0] for fn in gsesoft_fnlist]
    gsepatt = settings.gsepatt; rgse = re.compile(gsepatt)
    gseidlist = list(filter(rgse.match, gseidlist)) # valid GSE IDs
    print('Processing SOFT files for GSE list: '+str(gseidlist))
    softpatt = settings.softallpatt
    comppatt = settings.compsoftpatt; expsoftpatt = settings.expsoftpatt
    rsoft = re.compile(softpatt); rcompsoft = re.compile(comppatt)
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
        
        Extract GSM sample metadata from GSE SOFT files.
        
        Arguments: 
        * gsesoft_flist (list, optional) : List of gse SOFT files to process
        * softopenindex (str) : Index of label/tag to open entry, defaults
            to sample title section.
        * softcloseindex (str) : Index of label/tag to close entry, defaults 
            to close just before possible by-CpG methylation table. To include 
            possible methylation data table, change to '!sample_table_end'.
        * timestamp (str) : NTP timestamp version for expanded files.
        * rmtempdir (Bool.) : Whether to remove temp directory.
        * validate (Bool.) : Validate extracted GSM files against files in 
            gsm_soft directory?
        
        Returns:
        * newfilesd (dictionary), or error (null), generates GSM SOFT files 
            as a side effect.
    """
    eqfiltdict=get_queryfilt_dict()
    validgsmlist = list(set([gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ])); print("length validgsmlist : "+str(len(validgsmlist)))
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
    shuffle(gse_softlist); newfilesd = {} # new files, status dict to return
    print("new tempdir for writing soft files : "+str(temp_dir_make))
    print("length gse_softlist: "+str(len(gse_softlist)))
    rxopen = re.compile(softopenindex); rxclose = re.compile(softcloseindex)
    rxgsm = re.compile('GSM[0-9]*'); rxgsmfile = re.compile('.*GSM.*')
    for gse_softfile in gse_softlist:
        print("Beginning gse softfile : "+gse_softfile)
        newfilesd[gse_softfile] = []; openindex = []; closeindex = [];lsoft = []
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
            print("num : "+str(num)); print("openi : "+str(openi))
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
                    gsm_oldfile_path = ""; gsm_newfile_path = ""
                    gsm_softfn = gsmfile; gsmstr = gsm_softfn.split(".")[1]
                    print("gsmfile: "+str(gsmfile)); print("gsmstr : "+gsmstr)
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
        print("Removing tempdir..."); shutil.rmtree(temp_dir_make)
    return newfilesd 

def gsm_soft2json(gsm_softlist = [], scriptpath = settings.s2jscriptpath,
    gsm_jsonpath = settings.gsmjsonpath, gsm_softpath = settings.gsmsoftpath):
    """ gsm_soft2json
        
        Convert GSM soft file to JSON format Calls R script to coerce GSM soft 
        files (XML-like format) to valid JSON format.
        
        Arguments:
        * gsm_softlist (list, optional) : List of GSM soft filenames to process.
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
                rgsm = re.compile('.*GSM.*'); gsmid = gsm_softfn.split('.')[1]
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

if __name__ == "__main__":
    """ process_soft.py

    Process downloaded SOFT files
    
    """
    expand_soft(); extract_gsm_soft(); gsm_soft2json()