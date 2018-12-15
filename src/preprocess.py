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
from server import get_queryfilt_dict
from dl import gettime_ntp

""" preprocess.py
    Various wrappers to preprocessing utility functions.
    Prepares SOFT data for by-GSM input to MetaSRA-pipeline.
    Returns lists of existant SOFT and idat files to process.

"""

def get_mongofiles_for_preprocessing(idatsdir, softdir, filtresults = True):
    """ Get GSE and GSM IDs from MongoDB 
        Get most recent records for relevant GSE and GSM MongoDB entries
        Arguments
            *idatsdir (str): path to dir containing idat files
            * softdir (str): path to dir containing soft files
            * filtresults (T/F, Bool.): whether to pre-filter returned records
                on valid file status (e.g. path exists)
        Returns
            * doclist (list): list of relevant docs
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
    """ Runs preprocessing on soft files
        Arguments
            * filesdir (str): base files directory name
            * gse_softdir (str): dir to look for GSE soft files
            * gsm_softdir (str): dir for GSM soft files
            * gsm_jsondir (str): GSM JSON file dir
            * expand_soft_opt (Bool.): Option, whether to scan for/expand 
                compressed soft files
            * extract_gsm_opt (Bool. ): Option, whether to extract gsm-level 
                soft data
            * conv_json_opt (Bool.): Option, whether to convert gsm soft files 
                to json
            * msrap_prepare_opt (Bool.): Option, whether to prepare gsm json 
                files for MetaSRA-pipeline
        Returns
            * status (0 or 1): Whether soft list preprocessing completed
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
        extract_gsm_soft()
    # convert gsm soft to json
    try:
        gsm_softlist = os.listdir(gsm_softpath)
    except:
        print("Error finding GSM SOFT files. Returning...")
        return
    if conv_json_opt:
        if gsm_softlist:
            for gsmfile in gsm_softlist:
                gsm_soft2json()
        else:
            print('No GSM SOFT files detected, skipping JSON conversion...')
    # get the JSON filelist
    try:
        gsm_jsonpath = os.path.join(filesdir, gsm_jsondir)
        gsm_json_filelist = os.listdir(gsm_jsonpath)
    except:
        print("Error finding JSON files. Returning...")
        return
    if msrap_prepare_opt:
        if gsm_json_filelist:
            msrap_prepare_json(gsm_json_filelist = gsm_json_filelist)
        else:
            print('No GSM JSON files detected. Skipping MSRA-p preparation...')
    if run_msrap_opt:
        run_metasrapipeline()
    return

if __name__ == "__main__":
    idir = os.path.join('recount-methylation-files','idats')
    sdir = os.path.join('recount-methylation-files','gse_soft')
    pfd = preprocess.get_mongofiles_for_preprocessing(idatsdir=idir,
        softdir = sdir)
    process_gsesoft() # processes all soft files
    gse_soft_records = pfd['soft'] # select valid soft files 
    
def compile_rsheet(eqfiltd = get_queryfilt_dict(), sheetsdir = 'sheetfiles',
    sheetfn_ext = 'rsheet',msrapdir = 'gsm_msrap_outfiles',
    msrapfn_ext = '.soft', idatsfn_ext = 'idat',
    idatsdir = 'idats', filesdir = 'recount-methylation-files',
    timestamp = gettime_ntp()):
    """ Knits poised file data together into a minfi-ready sheet
        # grab msrap file list
        # grab idats file list
        # intersect files lists
        # subset eqfilt dict on gse
        # form and write new sheet, one per gse
    """
    rxmsrap = re.compile(".*soft$") # expanded soft file filter
    rxidat = re.compile(".*idat$") # expanded idat file filter
    msrap_path = os.path.join(filesdir, msrapdir)
    msrap_fnlist = os.listdir(msrap_path)
    idats_path = os.path.join(filesdir, idatsdir)
    idats_fnlist = os.listdir(idats_path)
    msrap_fnlist = list(filter(rxmsrap.match, msrap_fnlist))
    idats_fnlist = list(filter(rxidat.match, idats_fnlist))
    # extract gsm id
    gsmid_idatfn = [fn.split('.')[0] for fn in idats_fnlist] 
    gsmid_msrapfn = [fn.split('.')[1] for fn in msrap_fnlist] 
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
        lgrn = list(filter(rxgrn.match, lgsmid_idat))
        lred = list(filter(rxred.match, lgsmid_idat))
        if len(lgrn) == 1 and len(lred) == 1 and len(lgsmid_msrap) == 1:
            gsmvalid.append(gsmid)
    if gsmvalid:
        # form the rsheets with valid gsm ids
        lsheet = [] # list object to write rsheet, one row per gsmid
        lsheet_idats = [] # list obj, one row per idat chan
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
            gsm_mdat = ";".join(grows)
            # one entry per gsm
            lgsm = " ".join([gsmid, # gsm id
                gseid, # gse id
                ";".join(lgsmid_idat_fn), # idat names
                lgsmid_msrap,
                gsm_mdat # flattened json file
            ])
            print("lgsm : ",lgsm)
            # one entry per idat channel file
            print("grnmatch : ",str(list(filter(rxgrn.match, lgsmid_idat_fn))[0]))
            print("redmatch : ",str(list(filter(rxred.match, lgsmid_idat_fn))[0]))
            idatgrnfn = str(list(filter(rxgrn.match, lgsmid_idat_fn))[0])
            idatredfn = str(list(filter(rxred.match, lgsmid_idat_fn))[0])
            lgsm_grn = " ".join([lgsm, os.path.join(idats_path,idatgrnfn)])
            lgsm_red = " ".join([lgsm, os.path.join(idats_path, idatredfn)])
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
preprocess.compile_rsheet()


import process_idats
import re

process_idats.expand_idats()
preprocess.compile_rsheet()

filesdir = 'recount-methylation-files'
msrapdir = 'gsm_msrap_outfiles'
idatsdir = 'idats'

rxmsrap = re.compile(".*soft$")
rxidat = re.compile(".*idat$")
msrap_path = os.path.join(filesdir, msrapdir)
msrap_fnlist = os.listdir(msrap_path)
idats_path = os.path.join(filesdir, idatsdir)
idats_fnlist = os.listdir(idats_path)
msrap_fnlist = list(filter(rxmsrap.match, msrap_fnlist))
idats_fnlist = list(filter(rxidat.match, idats_fnlist))
# extract gsm id
gsmid_idatfn = [fn.split('.')[0] for fn in idats_fnlist] 
gsmid_msrapfn = [fn.split('.')[1] for fn in msrap_fnlist] 
# get the final validated gsm ids
gsmvalid = []
rxgrn = re.compile('.*Grn.idat$')
rxred = re.compile('.*Red.idat$')
gsmid = gsmid_idatfn[0]

rxgsm = re.compile(''.join(['.*',gsmid,'.*']))
lgsmid_idat = list(filter(rxgsm.match, idats_fnlist))
lgsmid_msrap = list(filter(rxgsm.match, gsmid_msrapfn))
lgrn = list(filter(rxgrn.match, lgsmid_idat))
lred = list(filter(rxred.match, lgsmid_idat))


preprocess.process_gsesoft()

idir = os.path.join('recount-methylation-files','idats')
sdir = os.path.join('recount-methylation-files','gse_soft')

pfd = preprocess.get_mongofiles_for_preprocessing(idatsdir=idir,
softdir = sdir)

preprocess.process_gsesoft()

pfd['idats'][list(pfd['idats'].keys())[0]]


# test eqfiltd filter
eqd = server.get_queryfilt_dict()
gseid = ""
gsmid = 'GSM999335'

list(key for key in eqd.keys() if gsmid in eqd[key])[0]

for key, values in eqd:
    if gsmid in values:
        gseid = key

for key, value in dd:

grows = []
for key in list(dd.keys()):
    kval = dd[key]
    if type(kval) is list:
        grows.append(";".join(kval))
    else:
        grows.append(":".join([k,dd[k]]))
gsm_mdat = ";".join(grows)



[":".join([key,dd[key]]) for key in dd]


"""
