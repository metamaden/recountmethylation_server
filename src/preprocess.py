#!/usr/bin/env python3

import pymongo
import sys
import os
import datetime
import inspect
import re
import json
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from process_soft import expand_soft, extract_gsm_soft, gsm_soft2json
from process_soft import msrap_prepare_json, run_metasrapipeline
from process_idats import expand_idats
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
from utilities import querydict
import settings
settings.init()

""" preprocess.py
    Wrapper functions for preprocessing idats and soft files, and functions to
    prepare data and valid sheet files for R analysis.
    Prepares SOFT data for by-GSM input to MetaSRA-pipeline, and runs the 
    pipeline.
    Functions:
        * get_mongofiles_for_preprocessing: Access latest file info for idats 
            and soft files, from RMDB. Includes validation of filepaths.
        * process_gsesoft: Wrapper function serving as preprocessing pipeline
            from gse soft files to gsm MetaSRA-pipeline outfiles (mapped sample
            metadata).
        * compile_rsheet: Compile a valid sheet/table of samples with info for 
            metadata and idat filenames, and collapsed form of the MetaSRA-
            pipeline output. Can be read into R/minfi.
"""

# def get_mongofiles_for_preprocessing(idatsdir, softdir, filtresults = True):
def get_mongofiles_for_preprocessing(filtresults = True):
    """ Get GSE and GSM IDs from MongoDB 
        Get most recent records for relevant GSE and GSM MongoDB entries
        Arguments
            *idatsdir (str) : Path to dir containing idat files
            * softdir (str) : Path to dir containing soft files
            * filtresults (T/F, Bool.) : Whether to pre-filter returned records
                on valid file status (e.g. if path exists).
        Returns
            * doclist object (list): List of relevant docs
    """
    # connect to mongodb
    # client = pymongo.MongoClient('localhost', 27017)
    client = pymongo.MongoClient(settings.rmdbhost, settings.rmdbport)
    dbcon = client.recount_methylation
    idatscon = dbcon.gsm.idats
    softcon = dbcon.gse.soft
    # grab unique gsm ids
    idatslist = list(idatscon.find())
    idatslist = [record for record in idatslist if 'gsmid' in record.keys()]
    gsmindex = list(set([record['gsmid'] for record in idatslist]))
    # grab unique gse ids
    softlist = list(softcon.find())
    softlist = [record for record in softlist if 'gseid' in record.keys()]
    gseindex = set([record['gseid'] for record in softlist])
        # filter all records for gsm on most recent update datetime
    idatrecordsfilt = {}
    # grab and filter idats file list
    for gsm in gsmindex:
        idatrecordsfilt[gsm] = []
        idatspreprocess = []
        recordsgsm = [record for record in idatslist if record['gsmid']==gsm]
        # handle each color chan type, red and grn
        idatsgsmgrn = [record for record in recordsgsm 
            if isinstance(record['date'],datetime.datetime)
            and re.search('Grn',os.path.basename(record['filepath']))
            and re.search('idat',os.path.basename(record['filepath']))
        ]
        idatsgsmred = [record for record in recordsgsm 
            if isinstance(record['date'],datetime.datetime)
            and re.search('Red',os.path.basename(record['filepath']))
            and re.search('idat',os.path.basename(record['filepath']))
        ]
        if idatsgsmgrn and idatsgsmred:
            ifiltgrn = sorted(idatsgsmgrn, key=lambda k: k['date'])[-1]
            ifiltred = sorted(idatsgsmred, key=lambda k: k['date'])[-1]
            # check that filepaths exist
            osg = os.path.exists(os.path.join(settings.idatspath,
                os.path.basename(ifiltgrn['filepath']
                )
            )
            )
            osr = os.path.exists(os.path.join(settings.idatspath,
                os.path.basename(ifiltred['filepath']
                )
            )
            )
            if osg:
                idatrecordsfilt[gsm].append(ifiltgrn)
            else:
                idatrecordsfilt[gsm].append('invalidpathgrn')
            if osr:
                idatrecordsfilt[gsm].append(ifiltred)
            else:
                idatrecordsfilt[gsm].append('invalidpathred')
        else:
            idatrecordsfilt[gsm].append('missingidat')
    # grab and filter soft file list
    softrecordsfilt = {}
    for gse in gseindex:
        softrecordsfilt[gse] = []
        recordgse = [record for record in softlist if record['gseid']==gse]
        gsesoftlist = [record for record in recordgse 
            if isinstance(record['date'],datetime.datetime)
            and re.search('soft',os.path.basename(record['filepath']))
        ]
        if gsesoftlist:
            gsesoftfilt = sorted(gsesoftlist, key=lambda k: k['date'])[-1]
        else:
            softrecordsfilt[gse].append(True)
        if gsesoftfilt:
            ossoft = os.path.exists(os.path.join(settings.gsesoftpath,
                os.path.basename(gsesoftfilt['filepath'])
                )
            )
            ossoftfn = os.path.exists(gsesoftfilt['filepath'])
            if ossoft or ossoftfn:
                softrecordsfilt[gse].append(gsesoftfilt)
            else:
                softrecordsfilt[gse].append(False)
        else:
            softrecordsfilt[gse].append(None)
    # return filtered file lists as dictionary
    drfiles = {}
    if not filtresults:
        #idatrec_resultfilt = []
        idatrec_resultfilt = [record for record in idatrecordsfilt 
            if not 'invalidpathgrn' in record
            and not 'invalidpathred' in record
            and not 'missingidat' in record
        ]
        # softrec_resultfilt = []
        softrec_resultfilt = [record for record in softrecordsfilt
            if not 'novalidrecordsoft' in record
            and not 'invalidpathsoft' in record
        ]
        drfiles['idats'] = idatrec_resultfilt
        drfiles['soft'] = softrec_resultfilt
    else:
        drfiles['idats'] = idatrecordsfilt
        drfiles['soft'] = softrecordsfilt
    return drfiles

#def process_gsesoft(filesdir = 'recount-methylation-files',
#    rmcompressed_gsesoft = False,
#    gse_softdir = 'gse_soft', gsm_softdir = 'gsm_soft', 
#    gsm_jsondir = 'gsm_json', expand_soft_opt = True, extract_gsm_opt = True, 
#    conv_json_opt = True, msrap_prepare_opt = True):
def process_gsesoft(rmcompressed_gsesoft=False, expand_soft_opt=True, 
    extract_gsm_opt=True, conv_json_opt=True, msrap_prepare_opt=True):
    """ Wrapper to preprocess GSE soft files and run MetaSRA-pipeline on GSM 
        data.
        Arguments
            * filesdir (str): Root files directory name.
            * rmcompressed_gsesoft (True/False, bool.) : Whether to remove the 
                compressed GSE soft files after they have been expanded.
            * gse_softdir (str): Directory to look for GSE soft files.
            * gsm_softdir (str): Destination directory for GSM soft files.
            * gsm_jsondir (str): Destination directory for GSM JSON files.
            * expand_soft_opt (True/False, Bool.): Option, whether to scan 
                for/expand compressed soft files.
            * extract_gsm_opt (True/False, Bool. ): Option, whether to extract 
                GSM-level soft data from GSE soft file(s).            
            * conv_json_opt (True/False, Bool.): Option, whether to convert GSM 
                soft files to JSON files.
            * msrap_prepare_opt (True/False, Bool.): Option, whether to prepare 
                GSM JSON files for MetaSRA-pipeline by aggregating sample JSON
                files.
        Returns
            * null
    """
    # filter softlist on valid files
    #gse_softpath = os.path.join(filesdir, gse_softdir)
    #gsm_softpath = os.path.join(filesdir, gsm_softdir)
    gse_softpath = settings.gsesoftpath
    gsm_softpath = settings.gsmsoftpath
    # expand all gse soft files at target dir
    if expand_soft_opt:
        expand_soft(rmcompressed = rmcompressed_gsesoft)
    # extract all gsm soft file metadata from expanded soft files
    if extract_gsm_opt:
        os.makedirs(gsm_softpath, exist_ok = True) # mkdir noclobber
        # get snapshot of current gsm soft destdir
        gsmsoft_currentdirlist = os.listdir(gsm_softpath) 
        extract_gsm_soft()
        # compare new snapshot to old, retain only changed files
        gsmsoft_newdirlist = os.listdir(gsm_softpath) 
        gsmsoft_difdirlist = [gsmfn for gsmfn in gsmsoft_newdirlist 
            if not gsmfn in gsmsoft_currentdirlist
        ]
        print(gsmsoft_difdirlist)
        # if new files extracted, continue with MetaSRA-pipeline prep/analysis
        if gsmsoft_difdirlist and len(gsmsoft_difdirlist)>0:
            # convert gsm soft to json
            gsm_softlist = gsmsoft_difdirlist
            if conv_json_opt:
                if gsm_softlist:
                    jsonconv_list = gsm_soft2json(gsm_softlist=gsm_softlist)
                else:
                    print("No GSM SOFT files detected, skipping JSON "
                        +"conversion...")
            jd = jsonconv_list[0]
            gsm_json_fnlist = [keyfn for keyfn in jd 
                if jd[keyfn] == 1
            ]
            # either prepare aggregated json files, or run each one processively
            # grab the list of new, successfully created json files
            if msrap_prepare_opt:
                if gsm_json_fnlist:
                    msrap_prepare_json(gsm_json_filelist = gsm_json_fnlist)
                else:
                    print("No GSM JSON files listed. Skipping MSRA-pipeline "
                        +"preparation...")
            else:
                if run_msrap_opt:
                    run_metasrapipeline(json_flist = gsm_json_fnlist)
    return 

def compile_rsheet(eqfiltd=get_queryfilt_dict(), sheetfn_ext='rsheet', 
    msrapfn_ext='msrapout', msrapfn='msrapout', idatsfn_ext='idat',
    timestamp=gettime_ntp()):
    """ Knits poised file data together into a sheet to be read into R using 
        minfi.
        Steps taken include: 
            1. Grab msrap file list
            2. Grab idats file list
            3. Intersect files lists
            4. Subset eqfilt dict on gse
            5. Form and write new sheet files, one per gse
        Arguments
            * eqfiltd (function or dictionary) : Equery filter dictionary object.
            * sheetsdir (str) : Directory to write new sheet files.
            * sheetfn_ext (str) : Filename extension for new sheet files.
            * msrapdir (str) : Directory containing MetaSRA-pipeline datafiles.
            * msrapfn_ext (str) : Filename extension of valid MetaSRA-pipeline
                datafiles.
            * idatsfn_ext (str) : Filename extension of valid idat files.
            * idatsdir (str) : Name of directory containing GSM idat files.
            * filesdir (str) : Root name of directory containing database files.
            * timestamp (str) : NTP timestamp for file versioning.
            * msrapfn (str) : File name stem for MetaSRA-pipeline files
        Returns:
            * null, produces sheet files as a side effect.
    """
    # form the sheet path and make dir as needed
    # sheetspath = os.path.join(filesdir, sheetsdir)
    sheetspath = settings.sheetspath
    os.makedirs(sheetspath, exist_ok = True)
    sheets_fpath = os.path.join(sheetspath,
        ".".join([timestamp, sheetfn_ext])
        )
    # form msrap and idat paths and get filenames
    # msrap_path = os.path.join(filesdir, msrapdir)
    msrap_path = settings.gsmmsrapoutpath
    rxmsrap = re.compile(".*"+msrapfn_ext+"$")
    msrap_fnlist = list(filter(rxmsrap.match, os.listdir(msrap_path)))
    print("msrap_fnlist : "+str(msrap_fnlist))
    # idats fn
    # idats_path = os.path.join(filesdir, idatsdir)
    idats_path = settings.idatspath
    rxidat = re.compile(".*"+idatsfn_ext+"$")
    idats_fnlist = list(filter(rxidat.match, os.listdir(idats_path)))
    # extract gsm ids
    rxgsm = re.compile(".*GSM[0-9]")
    idats_splitlist = [idatfn.split(".")[0]
        for idatfn in idats_fnlist
        if len(idatfn.split("."))>1
    ]
    idats_gsmlist_filt = list(set(filter(rxgsm.match,idats_splitlist))) # unique gsm ids
    msrap_splitlist = [msrapfn.split(".")[1]
        for msrapfn in msrap_fnlist
        if len(msrapfn.split("."))>1
    ]
    msrap_gsmlist_filt = list(set(filter(rxgsm.match,msrap_splitlist))) # unique gsm ids
    print("idats_gsmlist_filt : "+str(idats_gsmlist_filt))
    print("msrap_gsmlist_filt : "+str(msrap_gsmlist_filt))
    gsmvalid = [gsmid for gsmid in msrap_gsmlist_filt if gsmid in idats_gsmlist_filt]
    if len(gsmvalid)>0:
        rxgrn = re.compile(".*Grn.idat$")
        rxred = re.compile(".*Red.idat$")
        lsheet = [] # list object to write rsheet, one row per gsmid
        # append colnames
        lsheet.append(" ".join(["gsmid",
            "gseid",
            "idats_fn",
            "msrapmd_fn",
            "msrapmd_flatjson",
            "SENTRIX_ID",
            "ARRAY_ID",
            "Basename"]))
        lsheet[0] = lsheet[0]+"\n"
        for gsmid in gsmvalid:
            # compile the file info for this gsm
            rxgsmi = re.compile(".*"+gsmid+".*")
            gsmi_idats = list(filter(rxgsmi.match,idats_fnlist))
            gsmi_red_idats = list(filter(rxred.match,gsmi_idats))
            gsmi_grn_idats = list(filter(rxgrn.match,gsmi_idats))
            # get the latest file versions
            gsmi_red_pattern = gsmi_red_idats[0].split(".")[2]
            gsmi_grn_pattern = gsmi_grn_idats[0].split(".")[2]
            gsmi_red_latest = getlatest_filepath(filepath=idats_path, 
                    filestr=gsmi_red_pattern, embeddedpattern=True
                )
            gsmi_grn_latest = getlatest_filepath(filepath=idats_path,
                    filestr=gsmi_grn_pattern,embeddedpattern=True
                )
            # get the latest msrap file
            gsmi_msrap_latest = getlatest_filepath(filepath=msrap_path,
                    filestr=gsmid,embeddedpattern=True
                )
            print(gsmi_msrap_latest)
            if (gsmi_red_latest and not gsmi_red_latest == 0 and gsmi_grn_latest 
                and not gsmi_grn_latest == 0 and gsmi_msrap_latest 
                and not gsmi_msrap_latest == 0):    
                # form the rsheets with valid gsm ids
                with open(gsmi_msrap_latest, 'r') as msrapmd:
                    gsmi_metadata_dict = json.load(msrapmd)
                gsmi_md = gsmi_metadata_dict[0] # weird dictionary
                grows = []
                for key in list(gsmi_md.keys()):
                    kval = gsmi_md[key]
                    if type(kval) is list:
                        grows.append(";".join(kval))
                    else:
                        grows.append(":".join([str(key),str(gsmi_md[key])]))
                gsmi_mdvar = "'"+";".join(grows)+"'"
                # grab the gse id for this gsm
                gseid = str([gsek for gsek in list(eqfiltd.keys()) 
                            if gsmid in eqfiltd[gsek]
                            ][0]
                        )
                # make the gsm arrays path Basename for minfi
                gsmi_bn = "_".join(gsmi_red_latest.split("_")[0:3])
                # one entry per gsm
                lgsmi = " ".join([gsmid, # gsm id
                    gseid, # gse id
                    ";".join([os.path.basename(gsmi_red_latest),
                        os.path.basename(gsmi_grn_latest)
                        ]
                    ), # idat filenames
                    os.path.basename(gsmi_msrap_latest), # metadata filename
                    gsmi_mdvar, # flattened json file
                    os.path.basename(gsmi_red_latest).split("_")[-2], # sentrix id
                    os.path.basename(gsmi_red_latest).split("_")[-3], # array id
                    gsmi_bn # minfi path Basename, for arrays
                ])
                lgsmi = lgsmi+"\n"
                lsheet.append(lgsmi)
    else:
        print("No valid GSM IDs detected. Check idats and MetaSRA-pipeline GSM "
            +"files directories.")
        return 0
    # write the final sheet files
    with open(sheets_fpath,'w') as fsheet:
        for item in lsheet:
            fsheet.write(item)
    
    return lsheet     

if __name__ == "__main__":
    #idir = os.path.join('recount-methylation-files','idats')
    #sdir = os.path.join('recount-methylation-files','gse_soft')
    idir = settings.idatspath
    sdir = settings.gsesoftdir
    pfd = get_mongofiles_for_preprocessing(idatsdir=idir, softdir = sdir)
    process_gsesoft() # processes all soft files
    gse_soft_records = pfd['soft'] # select valid soft files
