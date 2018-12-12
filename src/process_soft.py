import os
import sys
import re
import gzip
import shutil
import subprocess
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from server import get_queryfilt_dict
from dl import gettime_ntp

def expand_soft(soft_dir='gse_soft',rmcompressed=False):
    """ Automatically expand detected soft files
        Arguments
            * soft_dir : directory to search for soft files
            * rmcompressed : remove compressed files successfully expanded
        Returns
            * list : compressed filenames and statuslist
    """
    r1 = re.compile(".*soft.*")
    dirlist = os.listdir(soft_dir)
    soft_list = list(filter(r1.match, dirlist)) 
    r2 = re.compile(".*\.gz$")
    softl_compressed = list(filter(r2.match,soft_list))
    if len(softl_compressed)>0:
        statuslist = []
        for softcompfile in softl_compressed:
            with gzip.open(os.path.join(soft_dir,softcompfile), 'rb') as f_in:
                with open(os.path.join(soft_dir,softcompfile[:-3]), 'wb') as f_out:
                    try:
                        shutil.copyfileobj(f_in, f_out)
                        statuslist.append(True)
                    except shutil.Error as se:
                        statuslist.append(se)
        statusindices = [i for i, x in enumerate(statuslist) if x == True]
        rmsuccess = [softl_compressed[i] for i in statusindices]
        if rmcompressed & len(rmsuccess)>0:
            for compfilename in rmsuccess:
                os.remove(os.path.join(soft_dir,compfilename))
    else: 
        print("Error: no compressed soft files found at specified soft_dir.")
        return 
    return [softl_compressed,statuslist]

def extract_gsm_soft(gse_soft_dir = 'gse_soft',gsm_soft_dir = 'gsm_soft',
    temp_dir = 'temp', timestamp = gettime_ntp()):
    """ Extract GSM soft file sections from GSE soft file
        Arguments 
            * gse_soft_dir : directory to search for gse soft files
            * gsm_soft_dir : dir to write new gsm soft files
        Returns
            * null, generates gsm soft files as side effect
    """
    r1 = re.compile(".*soft$")
    dirlist = os.listdir(gse_soft_dir)
    soft_list = list(filter(r1.match, dirlist))
    for soft in soft_list:
        openindex = []
        closeindex = []
        rxopen = re.compile('!Sample_title')
        rxclose = re.compile('!sample_table_end')
        lsoft = []
        with open(softfile) as file:
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
            gsmsoftlines = lsoft[openi:closeindex[num]]
            gsm_newfile_path = os.path.join(
                gsm_soft_dir,
                timestamp+str(rxgsm.findall(gsmsoft[1])[0])+'.soft'
                )
            with open(gsm_newfile_path,"w+") as gsmfile:
                for line in gsmsoftlines:
                    gsmfile.write(line)
                gsmfile.close()
    return

def gsm_soft2json(gsm_soft_fn, 
    scriptpath = os.path.join('recount-methylation-server','src','soft2json.R'), 
    gsmsoft_dir = 'gsm_soft',
    filesdir = 'recount-methylation-files',
    gsm_json_destdir='gsm_json'):
    """ Convert GSM soft file to JSON format
        Call R script to process GSM soft (XML) to JSON format
        Arguments
            * gsm_soft_fn : name of GSM soft file to process
            * scriptpath : path to R script for JSON conversion
            * gsm_json_destdir : name of dest dir for new JSON files
        Returns
            * null, or error
    """
    gsmsoftpath = os.path.join(filesdir,gsmsoft_dir)
    gsmi_softpath = os.path.join(gsmsoftpath,gsm_soft_fn)
    jsonpath = os.path.join(filesdir,gsm_json_destdir)
    os.makedirs(jsonpath, exist_ok = True)
    try:
        if scriptpath:
            cmdlist = ['Rscript',scriptpath,gsm_soft_fn,jsonpath]
        else:
            cmdlist = ['Rscript','soft2json.R',gsm_soft_fn,jsonpath]
        subprocess.call(cmdlist,shell=False)
        #subprocess.check_call(cmdlist,executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        raise e

def msrap_prepare_json(gsm_json_filelist,gsm_json_dir='gsm_json',
    dest_dir='msrap_infiles',dest_filename='new_prepared_file'):
    """ Prepare GSM JSON metadata files for input to MetaSRA-pipeline
        Arguments
            * gsm_json_filelist : list of files to prepare
            * gsm_json_dir : directory of GSM JSON files to read
            * dest_dir : location to store new prepared file
            * dest_filename : prepared file name, to write
        Returns
            * null or error
    """
    os.makedirs(dest_dir, exist_ok=True)
    # manually add brackets '[' or ']', exclude from JSON concatenation
    with open(os.path.join(dest_dir,dest_filename),"w+") as preparedfile:
        preparedfile.write("[\n")
        for num, jsonfn in enumerate(gsm_json_filelist,0):
            with open(os.path.join(gsm_json_dir,jsonfn)) as jsonfile:
                for line in jsonfile:
                    if not ( line[0] == "]" or line[0] == "[" ):
                        if line == "  }\n" and num < (len(gsm_json_filelist)-1):
                            preparedfile.write("},\n")
                        else:
                            preparedfile.write(line)
        preparedfile.write("]")

def run_metasrapipeline(gsm_json_file,json_file_dir='msrap_torun',
    msrap_fn='new_msrap_outfile',
    msrap_dest_dir='msrap_output'):
    """ Run MetaSRA-pipeline on a GSM json file or list of JSON files
        Arguments
            * gsm_json_file : individual or group GSM JSON file
            * json_file_dir : dir of JSON file to map
            * msrap_fn : name of file to write new mapped MSRA-p. output
            * msrap_dest_dir : dest. dir. of final mapped output
        Returns
            * null or error, generating a new MetaSRA file as side effect 
    """
    os.makedirs(msrap_dest_dir, exist_ok=True)
    try:
        cmdlist = ['python',os.path.join('MetaSRA-pipeline','run_pipeline.py'),
        os.path.join(json_file_dir,gsm_json_file),'>',
        os.path.join(msrap_dest_dir,msrap_fn)]
        subprocess.call(cmdlist,shell=False)
    except subprocess.CalledProcessError as e:
        raise e
