#!/usr/bin/env python3

import os
import sys
import re
import gzip
import shutil
import subprocess
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from server import get_queryfilt_dict
from dl import gettime_ntp

def expand_soft(gse_softdir='gse_soft', files_dir = 'recount-methylation-files',
    rmcompressed = True):
    """ Automatically expand detected soft files
        Arguments
            * gse_softdir (str): directory to search for GSE soft files
            * files_dir (str): base files dir. to search for GSE gse_softdir
            * rmcompressed (T/F,Bool.): remove compressed files after expanded?
        Returns
            * list : compressed filenames and statuslist
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
    timestamp = gettime_ntp(), softcloseindex = '!sample_table_begin'):
    """ Extract GSM soft file sections from GSE soft file
        Arguments 
            * gse_softdir : directory to search for gse soft files
            * gsm_softdir : dir to write new gsm soft files
            * filesdir : base files dir to search/make soft dir(s)
            * timestamp : timestamp version for expanded files
            * softcloseindex (str): index of row to close entry, defaults to 
                just before possible by-CpG methylation table. To include this
                possible methylation data from the soft file, change to 
                '!sample_table_end'.
        Returns
            * null, generates gsm soft files as side effect
    """
    r1 = re.compile(".*soft$") # identify expanded GSE soft files
    gse_softpath = os.path.join(filesdir, gse_softdir)
    gse_soft_dirlist = os.listdir(gse_softpath)
    gse_softlist = list(filter(r1.match, gse_soft_dirlist))
    gsm_softpath = os.path.join(filesdir,gsm_softdir)
    os.makedirs(gsm_softpath, exist_ok = True)
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
        # may add filter based on edirect query file
        for num, openi in enumerate(openindex,0):
            # read gsm lines
            gsm_softlines = lsoft[openi:closeindex[num]]
            gsm_softfn = ".".join([timestamp, 
                str(rxgsm.findall(gsm_softlines[1])[0]),
                'soft'])
            newfilesd[gse_softfile].append(gsm_softfn)
            gsm_newfile_path = os.path.join(gsm_softpath, gsm_softfn)
            with open(gsm_newfile_path,"w+") as gsmfile:
                for line in gsm_softlines:
                    gsmfile.write(line)
                gsmfile.close()
    # return dict, keys GSE soft files, vals are new GSM soft files
    return newfilesd 

def gsm_soft2json(gsm_softdir = 'gsm_soft', gsm_jsondir = 'gsm_json',
    filesdir = 'recount-methylation-files',
    scriptpath = os.path.join('recount-methylation-server','src','soft2json.R')
    ):
    """ Convert GSM soft file to JSON format
        Call R script to process GSM soft (XML) to JSON format
        Arguments
            * gsm_softdir (str): GSM soft files dir to search
            * gsm_jsondir (str): Destination dir to store GSM JSON files
            * filesdir (str): base files dir to search
            * scriptpath : path to R script for JSON conversion
            * gsm_json_destdir : name of dest dir for new JSON files
        Returns
            * null, or error
    """
    gsm_jsonpath = os.path.join(filesdir, gsm_jsondir)
    os.makedirs(gsm_jsonpath, exist_ok = True)
    gsm_softpath = os.path.join(filesdir, gsm_softdir)
    # form list of gsm soft filenames
    gsm_fn_list = os.listdir(gsm_softpath)
    gsm_fn_list = [os.path.splitex(gsmname)[0] for gsmname in gsm_fn_list]
    for gsm_softfn in gsm_softfn_list:
        # run R script from provided path or present directory
        # R function takes 3 args: 1. gsm soft filename; 2. gsm soft path,
        # and 3. gsm json path
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
            #subprocess.check_call(cmdlist,executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            raise e

def msrap_prepare_json(gsm_json_filelist, gsm_jsondir = 'gsm_json',
    dest_dir = 'msrap_infiles', filesdir = 'recount-methylation-files',
    dest_filename = 'new_prepared_file'):
    """ Prepare GSM JSON metadata files for input to MetaSRA-pipeline
        Arguments
            * gsm_json_filelist (list): list of files to prepare
            * gsm_jsondir (str): directory of GSM JSON files to read
            * dest_dir (str): location to store new prepared file
            * filesdir (str): base files directory
            * dest_filename (str): filename stem for new written files
        Returns
            * null or error
    """
    os.makedirs(dest_dir, exist_ok=True)
    # manually add brackets '[' or ']', exclude from JSON concatenation
    with open(os.path.join(dest_dir,dest_filename), "w+") as preparedfile:
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

def run_metasrapipeline(gsm_jsondir = 'gsm_json',
    filesdir = 'recount-methylation-files', msrap_fn = 'msrapout', 
    msrap_destdir = 'gsm_msrap_outfiles', msrap_dir = '.'):
    """ Run MetaSRA-pipeline on available GSM JSON files
        Arguments
            * gsm_json_file (str): individual or group GSM JSON file
            * gsm_jsondir (str): dir of JSON file to map
            * msrap_fn (str): name of file to write new mapped MSRA-p. output
            * msrap_destdir (str): dest. dir. of final mapped output
            * msrap_path (str): path to MetaSRA-pipeline app, defaults to pwd
        Returns
            * Int (1) or error: Whether MetaSRA-pipeline successfully ran, 
                generating a new MetaSRA file as side effect 
    """
    # path to run_pipeline file
    msrap_path = os.path.join(msrap_dir, 'MetaSRA-pipeline', 'run_pipeline.py')
    # path to read valid gsm json files from
    gsm_jsonpath = os.path.join(filesdir, gsm_jsondir)
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
