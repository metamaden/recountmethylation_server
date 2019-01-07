#!/usr/bin/env python3

import os
import sys
import re
import gzip
import shutil
import subprocess
import filecmp
import tempfile
from random import shuffle
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict

""" process_soft.py
    Preprocess GSE/GSM Soft files for recount methylation samples. Includes GSM 
    Soft metadata extraction, methods to convert to JSON, and methods to process 
    samples in MetaSRA-pipeline.
    
    Notes:
    * To avoid propagating invalid (non-HM450k) experiments or samples, an 
        equeryfiltdict dictionary object is called, with valid GSE ids as keys
        and valid/filtered GSM ids listed as values.
    
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


def expand_soft(gse_softdir='gse_soft', softpatt='.*\.soft.*', comppatt='.*\.gz$', 
    expsoftpatt='.*\.soft$', rmcompressed=False, gsepatt='^GSE.*',
    files_dir='recount-methylation-files', eqfiltdict=get_queryfilt_dict()):
    """ Expand GSE soft files that aren't expanded yet.
        Arguments
            * gse_softdir (str) : Directory to search for GSE soft files.
            * rmcompressed (T/F,Bool.) : Remove compressed files after expanded?
            * files_dir (str) : Base files directory containing GSE gse_softdir.
        Returns
            * rsoftd (list) : Compressed filenames and statuslist.
    """
    validgselist = list(eqfiltdict.keys())
    # get all gse ids in gse soft directory
    gsesoft_fpath = os.path.join(files_dir, gse_softdir)
    gsesoft_fnlist = os.listdir(gsesoft_fpath)  
    gseidlist =  [fn.split('.')[0] for fn in gsesoft_fnlist]
    rgse = re.compile(gsepatt)
    gseidlist = list(filter(rgse.match, gseidlist)) # valid GSE IDs
    # print(str(gseidlist))
    # filter terms
    rsoft = re.compile(softpatt)
    rcompsoft = re.compile(comppatt)
    rexpsoft = re.compile(expsoftpatt)
    gsesoft_compflist = [] # filtered fn list of compressed soft files
    for gseid in gseidlist:
        if gseid in validgselist:
            gse_filtfn = []
            flatest = getlatest_filepath(
                    filepath=os.path.join(files_dir, gse_softdir), 
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
    # form the soft files dictionary
    rsoftd = {}
    softl_compressed = gsesoft_compflist
    if len(softl_compressed)>0:
        for softcompfile in softl_compressed:
            rsoftd[softcompfile] = [] # instantiate new dict val as list 
            with gzip.open(os.path.join(gsesoft_fpath, softcompfile), 'rb') as f_in:
                with open(os.path.join(gsesoft_fpath, softcompfile[:-3]), 'wb') as f_out:
                    try:
                        shutil.copyfileobj(f_in, f_out)
                        # if file successfully 
                        rsoftd[softcompfile].append(True)
                    except shutil.Error as se:
                        rsoftd[softcompfile].append(se)
        statusindices = [i for i, x in enumerate(rsoftd[softcompfile]) 
                if x == True
            ]
        rmsuccess = [softl_compressed[i] for i in statusindices]
        if rmcompressed and len(rmsuccess) > 0:
            for compfilename in rmsuccess:
                os.remove(os.path.join(gse_softpath,compfilename))
                # if compressed file removed, append True
                rsoftd[compfilename].append(True) 
    else: 
        print("No valid compressed soft files found at specified gse_softpath. "
            +"Are all valid soft files already expanded?")
        return None
    # return dictionary of compressed GSE soft files and statuses
    return rsoftd

def extract_gsm_soft(gsesoft_flist=[], gse_softdir='gse_soft',
    gsm_softdir='gsm_soft', filesdir='recount-methylation-files', 
    temp_dir='temp', timestamp=gettime_ntp(), softopenindex='.*!Sample_title.*', 
    softcloseindex='.*!Sample_data_row_count.*', validate=True, 
    tempdir='temp', eqfiltdict=get_queryfilt_dict(), qcprint=False):
    """ Extract GSM soft file sections from GSE soft files.
        Arguments 
            * gse_softdir (str) : Directory to search for GSE soft files.
            * gsm_softdir (str) : Directory to write new GSM soft files.
            * filesdir (str) : Name of root files directory.
            * temp_dir (str) : Name of temp directory for new file downloads.
            * timestamp (str) : NTP timestamp version for expanded files.
            * softopenindex (str) : Index of label/tag to open entry, defaults
                to sample title section.
            * softcloseindex (str) : Index of label/tag to close entry, defaults 
                to just before possible by-CpG methylation table. To include 
                this possible methylation data from the soft file, change to 
                '!sample_table_end'.
            * eqfiltdict (dict) : Dictionary from equery filter, where keys are
                valid HM450k experiment ids (GSE ids) and values are lists of 
                valid/filtered HM450k samples (GSM ids). 
            * qcprint (Bool.) : whether to print status texts for QC.
        Returns
            * newfilesd (dictionary), or error (null), generates GSM soft files 
                as a side effect.
    """
    validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ]
    if qcprint:
        print("length validgsmlist : "+str(len(validgsmlist)))
    rvalidsoft = re.compile(".*soft$") # identify expanded GSE soft files
    gse_softpath = os.path.join(filesdir, gse_softdir) # GSE softpath
    gsm_softpath = os.path.join(filesdir, gsm_softdir) # GSM softpath
    gsmsoft_tempdir = os.path.join(filesdir,temp_dir)
    os.makedirs(gsm_softpath, exist_ok=True)
    os.makedirs(gsmsoft_tempdir, exist_ok=True)
    temp_dir_make = tempfile.mkdtemp(dir=gsmsoft_tempdir)
    if not gsesoft_flist or len(gsesoft_flist)==0:
        gse_soft_dirlist = os.listdir(gse_softpath)
    else:
        gse_soft_dirlist = gsesoft_flist
        gse_soft_dirlist = [gsefile for gsefile in gse_soft_dirlist
            if os.path.exists(os.path.join(gse_softpath), gsefile)
        ]
    gse_softlist = list(filter(rvalidsoft.match, gse_soft_dirlist))
    shuffle(gse_softlist)
    gsmsoft_destpath = os.path.join(filesdir, gsm_softdir) # file dest path
    newfilesd = {} # new files, status dictionary to return
    if qcprint:
        print("new tempdir for writing soft files : "+str(temp_dir_make))
        print("length gse_softlist: "+str(len(gse_softlist)))
    rxopen = re.compile(softopenindex)
    rxclose = re.compile(softcloseindex)
    rxgsm = re.compile('GSM[0-9]*')
    rxgsmfile = re.compile('.*GSM.*')
    for gse_softfile in gse_softlist:
        if qcprint:
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
        if qcprint:
            print("for gse, found n = "+str(len(openindex))+" valid open "
                +"indices..")
            print("for gse, found n = "+str(len(lsoft))+" soft file lines. "
                +"Continuing...")
        for num, openi in enumerate(openindex,0):
            if qcprint:
                print("num : "+str(num))
                print("openi : "+str(openi))
                print("closei : "+str(closeindex[num]))
            gsm_softlines = lsoft[openi:closeindex[num]] # read gsm lines
            gsmid_lines = [line for line in gsm_softlines
                if '!Sample_geo_accession' in line
            ]
            if len(gsmid_lines)==1:
                gsmid = str(rxgsm.findall(gsmid_lines[0])[0])
                print("Found GSM id : "+gsmid)
                gsm_softfn = ".".join([timestamp, gsmid, 'soft'])
                if qcprint:
                    print("GSM id found : "+gsmid)
                    print(gsmid+" in valid list... "+str(gsmid in validgsmlist))
                if gsmid in validgsmlist:
                    newfilesd[gse_softfile].append(gsm_softfn)
                    gsm_newfile_path = os.path.join(temp_dir_make, gsm_softfn)
                    write = [gsmfile.write(line) for line in open(gsm_newfile_path,"w+")]
                    open(gsm_newfile_path,"w+").write("\n".join(gsm_softlines))
                    #with open(gsm_newfile_path,"w+") as gsmfile:
                    #    for line in gsm_softlines:
                    #        gsmfile.write(line)
                    #    writeobj =  [gsmfile.write(line) for line in gsmfile]
                else: 
                    print("GSM id :"+gsmid+" is not a valid HM450k sample. "
                        +"Continuing...")
            else:
                print("GSM soft lines malformed! Continuing...")
    if qcprint:
        print("newfilesd : "+str(newfilesd))
    if validate:
        if qcprint:
            print(list(newfilesd.keys()))
        for gse_softfn in list(newfilesd.keys()):
            gsmfilelist = list(filter(rxgsmfile.match, newfilesd[gse_softfn]))
            if gsmfilelist and len(gsmfilelist)>0:
                if qcprint:
                    print(str(gsmfilelist))
                for gsmfile in gsmfilelist:
                    gsm_oldfile_path = ""
                    gsm_newfile_path = ""
                    gsm_softfn = gsmfile
                    gsmstr = gsm_softfn.split(".")[1]
                    if qcprint:
                        print("gsmfile: "+str(gsmfile))
                        print("gsmstr : "+gsmstr)
                    gsm_newfile_path = os.path.join(temp_dir_make, gsm_softfn)
                    gsm_oldfile_path = getlatest_filepath(
                            filepath=gsmsoft_destpath, filestr=gsmstr, 
                            embeddedpattern=True, tslocindex=0
                        )
                    if qcprint:
                        print("gsm_oldfile_path : "+str(gsm_oldfile_path))
                        print("gsm_newfile_path : "+str(gsm_newfile_path))
                    if gsm_oldfile_path and not gsm_oldfile_path == 0:
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
        shutil.rmtree(temp_dir_make) # remove tempdir if validate true
    # return dict, keys GSE soft files, vals are new GSM soft files
    return newfilesd 

def gsm_soft2json(gsm_softlist=[], gsm_softdir='gsm_soft', 
    gsm_jsondir='gsm_json', filesdir='recount-methylation-files',
    scriptpath=os.path.join('recount-methylation-server','src','soft2json.R'),
    eqfiltdict = get_queryfilt_dict()):
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
            * eqfiltdict (dict) : Dictionary from equery filter, where keys are
                valid HM450k experiment ids (GSE ids) and values are lists of 
                valid/filtered HM450k samples (GSM ids). 
        Returns
            * rlist object (list) of converted files and statuses, or error, 
                generates GSM JSON files as a side effect.
    """
    validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ]
    if not os.path.exists(scriptpath):
        print("Error: JSON conversion script not detected at path. Returning.")
        return None
    gsm_jsonpath = os.path.join(filesdir, gsm_jsondir)
    gsm_softpath = os.path.join(filesdir, gsm_softdir)
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
        for gsm_softfn in gsm_softfn_list:
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
                        except subprocess.CalledProcessError as e:
                            statd[gsm_softfn].append(None)
                            statd[gsm_softfn].append(e)
                else:
                    print("No latest JSON files found for GSM id: "+gsmid
                        +". Continuing...")
                    try:
                        cmdlist = ['Rscript', scriptpath, gsm_softfn, 
                            gsm_softpath, gsm_jsonpath
                        ]
                        subprocess.call(cmdlist,shell=False)
                        statd[gsm_softfn].append(True)
                    except subprocess.CalledProcessError as e:
                        statd[gsm_softfn].append(None)
                        statd[gsm_softfn].append(e)
            else:
                statd[gsm_softfn].append(None)
                print("GSM id : "+gsmid+" is not a valid HM450k sample. "
                    +"Continuing...")
    else:
        print("Error: No valid GSM Soft files to process from list. Returning.")
        return None
    # tally new json files generated
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

def run_metasrapipeline(json_flist=[], gsm_jsondir='gsm_json',
    filesdir='recount-methylation-files', msrap_fn='msrapout', 
    msrap_destdir='gsm_msrap_outfiles', msrap_dir='.',
    timestamp=gettime_ntp(), eqfiltdict=get_queryfilt_dict()):
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
    validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist
    ]
    # path to run_pipeline file
    msrap_path = os.path.join(msrap_dir, 'MetaSRA-pipeline', 'run_pipeline.py')
    # path to read valid gsm json files from
    gsm_jsonpath = os.path.join(filesdir, gsm_jsondir)
    # if filenames provided, form list, else list json dir contents
    if json_flist and len(json_flist)>0:
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
        gsmid = gsm_json_fn.split('.')[1]
        if gsmid in validgsmlist:
            outfn = os.path.splitext(gsm_json_fn)[0] # fn without extension
            gsmjson_readpath = os.path.join(gsm_jsonpath, gsm_json_fn)
            gsm_msrapout_writepath = os.path.join(msrap_destpath,
                ".".join([timestamp,outfn,msrap_fn]))
            try:
                cmdlist = ['python2',
                    msrap_path,
                    gsmjson_readpath,
                    gsm_msrapout_writepath
                    ]
                subprocess.call(cmdlist,shell=False)
                msrap_statlist.append(True)
            except subprocess.CalledProcessError as e:
                msrap_statlist.append(e)
        else:
            msrap_statlist.append(None)
            print("GSM id : "+gsmid+" is not a valid HM450k sample. "
                +"Continuing...")
    return msrap_statlist

def msrap_screens(json_flist=[], nscreensi=50, nmaxscreens=20, srcdir='src', 
    filesdir='recount-methylation-files', gsmjsondir='gsm_json',
    gsmsoftdir='gsm_soft', serverfilesdir='recount-methylation-server', 
    psoftfn='process_soft.py', gsmmsrapoutdir='gsm_msrap_outfiles',
    jsonfnpattern='.*json$', msrapoutfnpattern='.*\.msrapout$', qcprint=False):
    """ Deploy multiple screens processing GSM JSON files with MetaSRA-pipeline.
        Only process GSMs for which no MetaSRA-pipeline outfiles already exist.
        Notes:
            * Deployed screens should automatically quit once process is finished
            * Arg 'nscreensi' is number of samples to process, per screen. The
                corresponding number of screens launched depends on the files 
                list processed.
            * Arg 'nmaxscreens' should be informed by the memory resources of 
                the individual system, to avoid resource overuse and diminishing
                compute spped (e.g. for slower systems a lower total nmaxscreens 
                should be used).
        Arguments:
            * json_flist (list) : List of GSM JSON filenames to process. If not 
                provided, function automatically detects any new GSM JSON files
                without available MetaSRA-pipeline outfiles.
            * nscreensi (int) : Number of samples to process per screen deployed.
            * nmaxscreens (int) : Limit to total screen session deployed. 
            * srcdir (str) : Source files (e.g. scripts) directory name.
            * filesdir (str) : Recount methylation files based directory name.
            * gsmjsondir (str) : Files directory name containing GSM JSON files.
            * gsmsoftdir (str) : Files directory name containing GSM Soft files.
            * serverfilesdir (str) : Recount methylation server directory name.
            * psoftfn (str) : Name of script for preprocessing soft files.
            * gsmmsrapoutdir (str) : MetaSRA-pipeline outfiles directory name.
            * qcprint (Bool) : Whether to print QC texts.
        Returns:
            (Null) Generates >=1 screens for preprocessing files list(s), as a 
                side effect.
    """
    # define filepaths
    psoftpath = os.path.join(serverfilesdir, srcdir, psoftfn)
    gsmsoftpath = os.path.join(filesdir, gsmsoftdir)
    gsmjsonpath = os.path.join(filesdir, gsmjsondir)
    gsmmsrapoutpath = os.path.join(filesdir, gsmmsrapoutdir)
    # define the file patterns for filtering
    rjson = re.compile(jsonfnpattern)
    rmsrapout = re.compile(msrapoutfnpattern)
    # generate fl list of valid json files (haven't been processed yet)
    fl = []
    if json_flist and len(json_flist)>0:
        jsonfnlist = list(filter(rjson.match, json_flist)) 
        jsongsmlist = [x.split('.')[1] for x in jsonfnlist]
    else:
        json_flist = os.listdir(gsmjsonpath)
        jsonfnlist = list(filter(rjson.match, json_flist)) 
        jsongsmlist = [x.split('.')[1] for x in jsonfnlist]
    # list all msrapout files, filtering on valid filenames
    msrapoutfnlist = os.listdir(gsmmsrapoutpath) 
    msrapoutfnlist = list(filter(rmsrapout.match, msrapoutfnlist))
    if qcprint:
        print("Found "+str(len(msrapoutfnlist))+" files with pattern "
            +msrapoutfnpattern+". Continuing...")
    # list msrap outfiles gsm ids, and filter json gsm ids
    msrapgsmlist = [x.split('.')[2] for x in msrapoutfnlist]
    gsmprocess = [g for g in jsongsmlist 
            if not g in msrapgsmlist and
            g[0:3]=='GSM'
        ]
    # grab json filenames from json dir for gsm ids in process list
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
        if qcprint:
            print("Appended file "+gjsonfn+" to files list to process. "
                +"Progress: "+str(index)+"/"+str(len(gsmprocess))+"="+
                +str(100*round(index/len(gsmprocess),2))+"%. Continuing...")
    # form list of fn lists based on nscreensi and indices/slices
    if fl:
        if qcprint:
            print("Forming list of fn lists for screen deployment...")
        ll = []
        rangelist = [i for i in range(0, len(fl), nscreensi)]
        for enum, i in enumerate(rangelist[:-1]):
            ll.append(fl[i:rangelist[enum+1]])
        if len(fl[rangelist[-1]::]) > 0:
            ll.append(fl[rangelist[-1]::])
    else:
        print("Error, no files list object to process. Returning...")
        return None
    if qcprint:
        print('screens ll list, len = ' + str(len(ll)))
    timestampi = gettime_ntp() # make single timestamp call, for all indices
    # deploy screen(s) running MetaSRA-pipeline
    if len(ll)>1:
        for loc, sublist in enumerate(ll):
            # each loc in ll represents a new screen index, check vs. screen max
            if loc <= nmaxscreens:
                strid = str(loc)
                cmdlist0 = ['screen', '-S', 'MetaSRApipeline'+strid, '-dm', 'python3',
                        psoftpath, '--msraplist', 
                        ' '.join(str(item) for item in sublist), '--ntptime', 
                        timestampi
                    ]
            subprocess.call(cmdlist0,shell=False)
    else:
        cmdlist0 = ['screen', '-S', 'MetaSRApipeline'+str(0), '-dm', 'python3',
                psoftpath, '--msraplist', ' '.join(str(item) for item in ll[0]),
                '--ntptime', timestampi
            ]
        subprocess.call(cmdlist0,shell=False)

def main(filesdir = 'recount-methylation-files'):
    """
    """
    return

if __name__ == "__main__":
    # code called by msrap_screens() function to run MetaSRA-pipeline on fn list
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--msraplist", type=str, required=True,
        default=None, help='Files list to process with MetaSRA-pipeline.')
    parser.add_argument("--ntptime", type=str, required=True,
        default=gettime_ntp(), help='NTP timestamp, as a string.')
    args = parser.parse_args()
    # process filenames into files list object
    flmsrap = [file for file in args.msraplist.split(' ')]
    # run msrap on files list
    run_metasrapipeline(json_flist=flmsrap, timestamp=args.ntptime)









