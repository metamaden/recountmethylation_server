#!/usr/bin/env python3

import os
import sys
import re
import gzip
import shutil
import subprocess
import filecmp
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from serverutilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict

def expand_soft(gse_softdir = 'gse_soft', rmcompressed = False, 
    files_dir = 'recount-methylation-files'):
    """ Automatically expand detected GSE soft files.
        Arguments
            * gse_softdir (str) : Directory to search for GSE soft files.
            * rmcompressed (T/F,Bool.) : Remove compressed files after expanded?
            * files_dir (str) : Base files directory containing GSE gse_softdir.
        Returns
            * rsoftd (list) : Compressed filenames and statuslist.
    """
    gse_softpath = os.path.join(files_dir, gse_softdir)
    r1 = re.compile(".*soft.*")
    dirlist = os.listdir(gse_softpath)
    soft_list = list(filter(r1.match, dirlist)) 
    r2 = re.compile(".*\.gz$")
    rsoftd = {}
    softl_compressed = list(filter(r2.match, soft_list))
    if len(softl_compressed)>0:
        statuslist = []
        for softcompfile in softl_compressed:
            rsoftd[softcompfile] = [] # instantiate new dict val as list 
            with gzip.open(os.path.join(gse_softpath, softcompfile), 'rb') as f_in:
                with open(os.path.join(gse_softpath, softcompfile[:-3]), 'wb') as f_out:
                    try:
                        shutil.copyfileobj(f_in, f_out)
                        rsoftd[softcompfile].append(1)
                    except shutil.Error as se:
                        rsoftd[softcompfile].append(se)
        statusindices = [i for i, x in enumerate(rsoftd[softcompfile]) 
            if x == 1
            ]
        rmsuccess = [softl_compressed[i] for i in statusindices]
        if rmcompressed and len(rmsuccess) > 0:
            for compfilename in rmsuccess:
                os.remove(os.path.join(gse_softpath,compfilename))
                rsoftd[compfilename].append('removedcompressedfile')
    else: 
        print("Error: no compressed soft files found at specified gse_softpath.")
        return 
    # return dictionary of compressed GSE soft files and statuses
    return rsoftd

def extract_gsm_soft(gse_softdir = 'gse_soft',gsm_softdir = 'gsm_soft',
    filesdir = 'recount-methylation-files', temp_dir = 'temp', 
    timestamp = gettime_ntp(), softcloseindex = '!sample_table_begin',
    validate = True, tempdir = 'temp'):
    """ Extract GSM soft file sections from GSE soft files.
        Arguments 
            * gse_softdir (str) : Directory to search for GSE soft files.
            * gsm_softdir (str) : Directory to write new GSM soft files.
            * filesdir (str) : Name of root files directory.
            * temp_dir (str) : Name of temp directory for new file downloads.
            * timestamp (str) : NTP timestamp version for expanded files.
            * softcloseindex (str) : Index of label/tag to close entry, defaults 
                to just before possible by-CpG methylation table. To include 
                this possible methylation data from the soft file, change to 
                '!sample_table_end'.
        Returns
            * newfilesd (dictionary), or error (null), generates GSM soft files 
                as a side effect.
    """
    r1 = re.compile(".*soft$") # identify expanded GSE soft files
    # gse_soft files path
    gse_softpath = os.path.join(filesdir, gse_softdir)
    gse_soft_dirlist = os.listdir(gse_softpath)
    gse_softlist = list(filter(r1.match, gse_soft_dirlist))
    gsm_softpath = os.path.join(filesdir,gsm_softdir)
    os.makedirs(gsm_softpath, exist_ok = True)
    # temp dir
    gsmsoft_tempdir = os.path.join(filesdir,temp_dir)
    os.makedirs(gsmsoft_tempdir, exist_ok=True)
    # gsm soft dest path
    gsmsoft_destpath = os.path.join(filesdir, gsm_softdir)
    newfilesd = {}
    for gse_softfile in gse_softlist:
        newfilesd[gse_softfile] = []
        openindex = []
        closeindex = []
        rxopen = re.compile('!Sample_title')
        rxclose = re.compile(softcloseindex)
        lsoft = []
        gse_softfile_path = os.path.join(gse_softpath, gse_softfile)
        with open(gse_softfile_path) as file:
            for num, line in enumerate(file, 0):
                if rxclose.search(line):
                    closeindex.append(num)
                if rxopen.search(line):
                    openindex.append(num)
                lsoft.append(line)
        rxgsm = re.compile('GSM[0-9]*')
        # may add filter here, using edirect query file
        for num, openi in enumerate(openindex,0):
            # read gsm lines
            gsm_softlines = lsoft[openi:closeindex[num]]
            gsm_softfn = ".".join([timestamp, 
                str(rxgsm.findall(gsm_softlines[1])[0]),
                'soft'])
            newfilesd[gse_softfile].append(gsm_softfn)
            # first write new files to temp dir
            gsm_newfile_path = os.path.join(gsmsoft_tempdir, gsm_softfn)
            with open(gsm_newfile_path,"w+") as gsmfile:
                for line in gsm_softlines:
                    gsmfile.write(line)
                gsmfile.close()
    if validate:
        print(list(newfilesd.keys()))
        for gsmfile in list(newfilesd.keys()):
            gsm_oldfile_path = ""
            gsm_newfile_path = ""
            gsm_softfn = newfilesd[gsmfile][0]
            gsmstr = gsm_softfn.split(".")[1]
            print("gsmstr : "+gsmstr)
            gsm_newfile_path = os.path.join(gsmsoft_tempdir, 
                gsm_softfn
                )
            gsm_oldfile_path = getlatest_filepath(filepath = gsmsoft_destpath,
                filestr = gsmstr, embeddedpattern = True, tslocindex = 0)
            print("gsm_oldfile_path : "+str(gsm_oldfile_path))
            print("gsm_newfile_path : "+str(gsm_newfile_path))
            if gsm_oldfile_path and not gsm_oldfile_path == 0:
                if filecmp.cmp(gsm_oldfile_path, gsm_newfile_path):
                    print("Redundant GSM soft file detected, removing...")
                    os.remove(gsm_newfile_path)
                    newfilesd[gsmfile].append(False)
                else:
                    print("New GSM soft file detected, moving from temp...")
                    shutil.move(gsm_newfile_path, os.path.join(
                            gsmsoft_destpath, 
                            os.path.basename(gsm_newfile_path))
                        )
                    newfilesd[gsmfile].append(True)
            else: 
                print("New GSM soft file detected, moving from temp...")
                shutil.move(gsm_newfile_path, os.path.join(
                            gsmsoft_destpath, 
                            os.path.basename(gsm_newfile_path))
                        )
                newfilesd[gsmfile].append(True)
        shutil.rmtree(gsmsoft_tempdir) # remove tempdir if validate true
    # return dict, keys GSE soft files, vals are new GSM soft files
    return newfilesd 

def gsm_soft2json(gsm_softlist = [], gsm_softdir = 'gsm_soft', 
    gsm_jsondir = 'gsm_json', filesdir = 'recount-methylation-files',
    scriptpath = os.path.join('recount-methylation-server','src','soft2json.R')
    ):
    """ Convert GSM soft file to JSON format
        Call R script to process GSM soft (XML) to JSON format
        Arguments
            * gsm_softlist (list) : Optional list of GSM soft filenames.
            * gsm_softdir (str): GSM soft files directory to search.
            * gsm_jsondir (str): Destination directory to store new JSON files.
            * filesdir (str): Base files directory to search.
            * scriptpath (str) : Path to R script for JSON conversion. If not 
                provided, automatically checks current working directory.
            * gsm_json_destdir (str) : Name of destination directory for new 
                JSON files.
        Returns
            * rlist object (list) of converted files and statuses, or error, 
                generates GSM JSON files as a side effect.
    """
    gsm_jsonpath = os.path.join(filesdir, gsm_jsondir)
    os.makedirs(gsm_jsonpath, exist_ok = True)
    rjsonlist_current = os.listdir(gsm_jsonpath) # get current json dir contents
    gsm_softpath = os.path.join(filesdir, gsm_softdir)
    # form list of gsm soft filenames
    if len(gsm_softlist)>0:
        gsm_softfn_list = [os.path.splitext(gsmname)[0] for gsmname 
            in gsm_softlist
        ]
    else:
        gsm_fn_list = os.listdir(gsm_softpath)
        gsm_softfn_list = [os.path.splitext(gsmname)[0] for gsmname 
            in gsm_fn_list
        ]
    statd = {} # return status dictionary from attempts at json conversion
    for gsm_softfn in gsm_softfn_list:
        # run R script from provided path or present directory
        # R function takes 3 args: 1. gsm soft filename; 2. gsm soft path,
        # and 3. gsm json path
        statd[gsm_softfn] = []
        try:
            if scriptpath:
                cmdlist = ['Rscript', scriptpath, gsm_softfn, gsm_softpath, 
                gsm_jsonpath
                ]
            else:
                cmdlist = ['Rscript', 'soft2json.R', gsm_softfn, gsm_softpath, 
                gsm_jsonpath
                ]
            subprocess.call(cmdlist,shell=False)
            statd[gsm_softfn].append(1)
            #subprocess.check_call(cmdlist,executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            statd[gsm_softfn].append(0)
            statd[gsm_softfn].append(e)
    # double check gsm json dir for new files
    rjsonlist_new = os.listdir(gsm_jsonpath)
    rjsonlist_return = [jfile for jfile in rjsonlist_new
        if not jfile in rjsonlist_current
    ]
    print("Detected N = "+str(len(rjsonlist_return))+" new json files in dir.")
    rlist = [statd, rjsonlist_return]
    return rlist

def msrap_prepare_json(gsm_json_filelist, gsm_jsondir = 'gsm_json',
    msrapin_dest_dir = 'msrap_infiles', filesdir = 'recount-methylation-files',
    dest_fnstem = 'new_prepared_file'):
    """ Prepare GSM JSON metadata files for input to MetaSRA-pipeline.
        Aggregates multiple specified GSM JSON files into a single composite 
        file that is valid for input to MetaSRA-pipeline.
        Arguments
            * gsm_json_filelist (list): List of GSM JSON files to prepare.
            * gsm_jsondir (str): Directory of GSM JSON files to read.
            * msrapin_dest_dir (str) : Destination directory for new prepared 
                JSON files.
            * dest_dir (str): Location to store new prepared file.
            * filesdir (str): Base files directory.
            * dest_filename (str): Filename stem for newly written files.
        Returns
            * null, or error, generates new aggregated JSON files as a side 
                effect.
    """
    msrap_infile_path = os.path.join(filesdir, msrapin_dest_dir)
    os.makedirs(filesdir, msrapin_dest_dir, exist_ok=True)
    # manually add brackets '[' or ']', exclude from JSON concatenation
    with open(os.path.join(msrap_infile_path,dest_fnstem), "w+") as preparedfile:
        preparedfile.write("[\n")
        for num, jsonfn in enumerate(gsm_json_filelist, 0):
            with open(os.path.join(gsm_jsondir, jsonfn)) as jsonfile:
                for line in jsonfile:
                    if not ( line[0] == "]" or line[0] == "[" ):
                        if line == "  }\n" and num < (len(gsm_json_filelist)-1):
                            preparedfile.write("},\n")
                        else:
                            preparedfile.write(line)
        preparedfile.write("]")

def run_metasrapipeline(json_flist = [], gsm_jsondir = 'gsm_json',
    filesdir = 'recount-methylation-files', msrap_fn = 'msrapout', 
    msrap_destdir = 'gsm_msrap_outfiles', msrap_dir = '.'):
    """ Run MetaSRA-pipeline on available GSM JSON files.
        Arguments
            * json_flist (list) : List of JSON filename(s) to process. If not 
                provided, automatically targets all JSON files at gsm_jsondir.
            * gsm_jsondir (str): Directory of JSON file(s) to map.
            * msrap_fn (str): Stem name of new files written.
            * msrap_destdir (str): Destination directory of final mapped output.
            * msrap_path (str): Path to MetaSRA-pipeline app, defaults to pwd
        Returns
            * msrap_statlist (list), Int (1) or error: Whether MetaSRA-pipeline 
                uccessfully ran, generating a new MetaSRA file as side effect. 
    """
    # path to run_pipeline file
    msrap_path = os.path.join(msrap_dir, 'MetaSRA-pipeline', 'run_pipeline.py')
    # path to read valid gsm json files from
    gsm_jsonpath = os.path.join(filesdir, gsm_jsondir)
    # if filenames provided, form list, else list json dir contents
    if len(json_flist)>0:
        rjson = re.compile(".*json$")
        gsm_json_fn_list = list(filter(rjson.match, json_flist)) 
    else:
        gsm_json_fn_list = os.listdir(gsm_jsonpath)
        rjson = re.compile(".*json$")
        gsm_json_fn_list = list(filter(rjson.match, gsm_json_fn_list)) 
    # path to write output files to
    msrap_destpath = os.path.join(filesdir, msrap_destdir)
    os.makedirs(msrap_destpath, exist_ok=True)
    msrap_statlist = []
    # iterate over valid filenames at gsm json dir
    for gsm_json_fn in gsm_json_fn_list:
        outfn = os.path.splitext(gsm_json_fn)[0] # fn without extension
        gsmjson_readpath = os.path.join(gsm_jsonpath, gsm_json_fn)
        gsm_msrapout_writepath = os.path.join(msrap_destpath,
            ".".join([outfn,msrapout]))
        try:
            cmdlist = ['python2',
                msrap_path,
                gsmjson_readpath,
                gsm_msrapout_writepath
                ]
            subprocess.call(cmdlist,shell=False)
            msrap_statlist.append(1)
        except subprocess.CalledProcessError as e:
            msrap_statlist.append(e)
    return msrap_statlist

""" Notes and Tutorial

import os
import sys
import re
import gzip
import shutil
import subprocess
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import process_soft

process_soft.gsm_soft2json()
process_soft.run_metasrapipeline()

"""
