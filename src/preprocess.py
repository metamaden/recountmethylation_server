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

""" preprocess.py
    Wrapper functions for preprocessing idats and soft files, and functions to
    prepare data and valid sheet files for R analysis.
    Prepares SOFT data for by-GSM input to MetaSRA-pipeline, and runs the 
    pipeline.
"""

def get_mongofiles_for_preprocessing(idatsdir, softdir, filtresults = True):
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
    client = pymongo.MongoClient('localhost', 27017)
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
            osg = os.path.exists(os.path.join(idatsdir,
                os.path.basename(ifiltgrn['filepath']
                )
            )
            )
            osr = os.path.exists(os.path.join(idatsdir,
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
            softrecordsfilt[gse].append('novalidrecordsoft')
        if gsesoftfilt:
            ossoft = os.path.exists(os.path.join(softdir,
                os.path.basename(gsesoftfilt['filepath'])
                )
            )
            ossoftfn = os.path.exists(gsesoftfilt['filepath'])
            if ossoft or ossoftfn:
                softrecordsfilt[gse].append(gsesoftfilt)
            else:
                softrecordsfilt[gse].append('invalidpathsoft')
        else:
            softrecordsfilt[gse].append('novalidrecordsoft')
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

def process_gsesoft(filesdir = 'recount-methylation-files',
    rmcompressed_gsesoft = False,
    gse_softdir = 'gse_soft', gsm_softdir = 'gsm_soft', 
    gsm_jsondir = 'gsm_json', expand_soft_opt = True, extract_gsm_opt = True, 
    conv_json_opt = True, msrap_prepare_opt = True):
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
    gse_softpath = os.path.join(filesdir, gse_softdir)
    gsm_softpath = os.path.join(filesdir, gsm_softdir)
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
                    jsonconv_list = gsm_soft2json(gsm_softlist = gsm_softlist)
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
    
def compile_rsheet(eqfiltd = get_queryfilt_dict(), sheetsdir = 'sheetfiles',
    sheetfn_ext = 'rsheet',msrapdir = 'gsm_msrap_outfiles',
    msrapfn_ext = '.soft', idatsfn_ext = 'idat',
    idatsdir = 'idats', filesdir = 'recount-methylation-files',
    timestamp = gettime_ntp(), msrapfn = 'msrapout'):
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
    rxmsrap = re.compile(".*"+msrapfn+"$") # expanded msrap out file filter
    rxidat = re.compile(".*idat$") # expanded idat file filter
    msrap_path = os.path.join(filesdir, msrapdir)
    msrap_fnlist = os.listdir(msrap_path)
    idats_path = os.path.join(filesdir, idatsdir)
    idats_fnlist = os.listdir(idats_path)
    msrap_fnlist = list(filter(rxmsrap.match, msrap_fnlist))
    idats_fnlist = list(filter(rxidat.match, idats_fnlist))
    # extract gsm id from arbitrary filename position
    rxgsmall = re.compile(".*GSM.*")
    gsmid_idatfn = [str(list(filter(rxgsmall.match,fn.split('.')))) 
        for fn in idats_fnlist
        ] 
    gsmid_msrapfn = [str(list(filter(rxgsmall.match,fn.split('.')))) 
        for fn in msrap_fnlist
        ] 
    # get the final validated gsm ids
    gsmvalid = []
    rxgrn = re.compile(".*Grn.idat$")
    rxred = re.compile(".*Red.idat$")
    for gsmid in gsmid_idatfn:
        print(gsmid)
        # check that valid grn and red channel files exist
        # check that gsm has available mapped metadata
        rxgsm = re.compile("".join([".*",gsmid,".*"]))
        lgsmid_idat = list(filter(rxgsm.match, idats_fnlist))
        lgsmid_msrap = list(filter(rxgsm.match, msrap_fnlist))
        # add filter on latest file versions
        lgrn = list(filter(rxgrn.match, lgsmid_idat))
        lred = list(filter(rxred.match, lgsmid_idat))
        if len(lgrn) == 1 and len(lred) == 1 and len(lgsmid_msrap) == 1:
            gsmvalid.append(gsmid)
    if gsmvalid:
        # form the rsheets with valid gsm ids
        lsheet = [] # list object to write rsheet, one row per gsmid
        lsheet_idats = [] # list obj, one row per idat chan
        # append colnames
        lsheet.append(" ".join(["gsmid",
            "gseid",
            "idats_fn",
            "metadata_fn",
            "msrap_metadata_flatjson"]))
        lsheet_idats.append(" ".join([lsheet[0],
            "idat_fpaths",
            "Sentrix_Position",
            "Sentrix_ID"]
            )
        )
        lsheet[0] = lsheet[0]+"\n"
        lsheet_idats[0] = lsheet_idats[0]+"\n"
        for gsmid in gsmvalid:
            lgsm = ""
            rxgsm = re.compile("".join([".*",gsmid,".*"]))
            lgsmid_idat = list(filter(rxgsm.match, gsmid_idatfn))
            lgsmid_idat_fn = list(filter(rxgsm.match, idats_fnlist))
            lgsmid_msrap = list(filter(rxgsm.match, gsmid_msrapfn))[0]
            lgsmid_msrap_fn = list(filter(rxgsm.match, msrap_fnlist))[0]
            print("lgsmid_idat : ",lgsmid_idat)
            print("lgsmid_idat_fn : ",lgsmid_idat_fn)
            print("lgsmid_msrap : ",lgsmid_msrap)
            print("lgsmid_msrap_fn : ",lgsmid_msrap_fn)
            # extract the gse id from eqfiltd
            gseid = list(key for key in eqfiltd.keys() 
                if gsmid in eqfiltd[key]
                )[0]
            # grab the flattened metadata
            gsm_msrap_filepath = os.path.join(msrap_path,lgsmid_msrap_fn)
            with open(gsm_msrap_filepath) as f:
                gsm_metadata_dict = json.load(f)
            gsm_metadata_dict = gsm_metadata_dict[0] # weird dictionary
            dd = gsm_metadata_dict
            grows = []
            for key in list(dd.keys()):
                kval = dd[key]
                if type(kval) is list:
                    grows.append(";".join(kval))
                else:
                    grows.append(":".join([str(key),str(dd[key])]))
            gsm_mdat = "'"+";".join(grows)+"'"
            # one entry per gsm
            lgsm = " ".join([gsmid, # gsm id
                gseid, # gse id
                ";".join(lgsmid_idat_fn), # idat filenames
                lgsmid_msrap_fn, # metadata filename
                gsm_mdat # flattened json file
            ])
            print("lgsm : ",lgsm)
            # one entry per idat channel file
            print("grnmatch : ",str(list(filter(rxgrn.match, lgsmid_idat_fn))[0]))
            print("redmatch : ",str(list(filter(rxred.match, lgsmid_idat_fn))[0]))
            idatgrnfn = str(list(filter(rxgrn.match, lgsmid_idat_fn))[0])
            idatgrn_idatfn = idatgrnfn.split(".")[-2] # get the array full name
            idatgrn_pos = idatgrn_idatfn.split("_")[-2] # sentrix id
            idatgrn_id = idatgrn_idatfn.split("_")[-3] # array id
            idatredfn = str(list(filter(rxred.match, lgsmid_idat_fn))[0])
            idatred_idatfn = idatredfn.split(".")[-2] # get the array full name
            idatred_pos = idatred_idatfn.split("_")[-2] # sentrix id
            idatred_id = idatred_idatfn.split("_")[-3] # array id
            lgsm_grn = " ".join([lgsm, 
                os.path.join(idats_path, idatgrnfn),
                idatgrn_pos,
                idatgrn_id
                ]
            )
            lgsm_red = " ".join([lgsm, 
                os.path.join(idats_path, idatredfn),
                idatred_pos,
                idatred_id
                ]
            )
            lsheet.append(lgsm+"\n")
            lsheet_idats.append(lgsm_grn+"\n")
            lsheet_idats.append(lgsm_red+"\n")
        # write the final sheet files
        sheetspath = os.path.join(filesdir, sheetsdir)
        sheets_fpath = os.path.join(sheetspath,
            ".".join([timestamp, sheetfn_ext])
            )
        sheets_fidatpath = os.path.join(sheetspath,
            ".".join([timestamp, sheetfn_ext, "idats"])
            )
        os.makedirs(sheetspath, exist_ok = True)
        with open(sheets_fpath, 'w') as fsheet:
            for item in lsheet:
                fsheet.write(item)
        with open(sheets_fidatpath, 'w') as fsheet:
            for item in lsheet_idats:
                fsheet.write(item)
    else:
        print("No valid GSM IDs detected. Check idats and MetaSRA-pipeline GSM "
            +"files directories.")
        return

if __name__ == "__main__":
    idir = os.path.join('recount-methylation-files','idats')
    sdir = os.path.join('recount-methylation-files','gse_soft')
    pfd = preprocess.get_mongofiles_for_preprocessing(idatsdir=idir,
        softdir = sdir)
    process_gsesoft() # processes all soft files
    gse_soft_records = pfd['soft'] # select valid soft files

""" Notes and Tutorial

import pymongo

client = pymongo.MongoClient('localhost', 27017)
cidats = client.recount_methylation.gsm.idats
csoft = client.recount_methylation.gse.soft
softlist = list(client.recount_methylation.gse.soft.find())
softlist = [record for record in softlist if 'gseid' in record.keys()]
gseindex = set([record['gseid'] for record in softlist])

#

import sys
import os
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import preprocess

preprocess.process_gsesoft()

process_idats.expand_idats()
preprocess.compile_rsheet()


"""
