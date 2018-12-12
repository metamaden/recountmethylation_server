#!/usr/bin/env python3
import pymongo
import sys
import os
import datetime
import inspect
import re
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from process_soft import expand_soft,extract_gsm_soft,gsm_soft2json
from process_soft import msrap_prepare_json,run_metasrapipeline
from server import get_queryfilt_dict

""" preprocess.py
    Various wrappers to preprocessing utility functions.
    Prepares SOFT data for by-GSM input to MetaSRA-pipeline.
    Returns lists of existant SOFT and idat files to process.

"""

def get_mongofiles_for_preprocessing(idatsdir,softdir):
    """ Get GSE and GSM IDs from MongoDB 
        Get most recent records for relevant GSE and GSM MongoDB entries
        Arguments
            *idatsdir: path to dir containing idat files
            * softdir: path to dir containing soft files
        Returns
            * doclist : list of relevant docs
    """
    # connect to mongodb
    client = pymongo.MongoClient('localhost', 27017)
    cidats = client.recount_methylation.gsm.idats
    csoft = client.recount_methylation.gse.soft
    # grab unique gsm ids
    idatslist = list(client.recount_methylation.gsm.idats.find())
    idatslist = [record for record in idatslist if 'gsmid' in record.keys()]
    gsmindex = list(set([record['gsmid'] for record in idatslist]))
    # grab unique gse ids
    softlist = list(client.recount_methylation.gse.soft.find())
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
    print(softrecordsfilt)
    for gse in gseindex:
        softrecordsfilt[gse] = []
        gsesoftlist = [record for record in recordsgsm 
            if isinstance(record['date'],datetime.datetime)
            and re.search('Grn',os.path.basename(record['filepath']))
            and re.search('idat',os.path.basename(record['filepath']))
        ]
        if gsesoftlist:
            gsesoftfilt = sorted(gsesoftlist, key=lambda k: k['date'])[-1]
        else:
            softrecordsfilt[gse].append('novalidrecordsoft')
        if gsesoftfilt:
            ossoft = os.path.exists(os.path.join(softdir,
                os.path.basename(gsesoftfilt['filepath']
                )
                )
            )
            if ossoft:
                softrecordsfilt[gse].append(gsesoftfilt)
            else:
                softrecordsfilt[gse].append('invalidpathsoft')
        else:
            softrecordsfilt[gse].append('novalidrecordsoft')
    # return filtered file lists as dictionary
    drfiles = {}
    drfiles['idats'] = idatrecordsfilt
    drfiles['soft'] = softrecordsfilt
    return drfiles

def preprocess_softlist(softlist,
    filesdir = 'recount-methylation-files',
    softdir = 'gse_soft',
    expand_soft = True, extract_gsm = True, conv_json = True, 
    msrap_preprare = True):
    """ Runs preprocessing on soft files
        Arguments
            * softdir: directory to look for soft files
            * softdir: dir to look for soft files
            * expand_soft: whether to scan for/expand compressed soft files
            * extract_gsm: whether to extract gsm-level soft data
            * conv_json: whether to convert gsm soft files to json
            * msrap_prepare: whether to prepare gsm json files for 
                MetaSRA-pipeline
        Returns
            * status (0 or 1)
    """
    # filter softlist on valid files

    softpath = os.path.join(filesdir,softdir)
    if expand_soft:
        expand_soft(softdir, rmcompressed = False)
    # extract gsm soft file metadata
    if extract_gsm:
        gsmsoftpath = os.path.join('recount-methylation-files','gsm_soft')
        os.makedirs(gsmsoftpath, exist_ok = True) # mkdir noclobber
        extract_gsm_soft(gse_soft_dir = softpath,
            gsm_soft_dir = gsmsoftpath,
            eq
            )
    # convert gsm soft to json
    try:
        gsmsoftlist = os.listdir(gsmsoftpath)
    except:
        print("Error finding GSM SOFT files. Returning...")
        return
    if conv_json:
        if gsmsoftlist:
            for gsmfile in gsmsoftlist:
                gsm_soft2json(gsm_soft_fn=gsmfile)
        else:
            print('No GSM SOFT files detected, skipping JSON conversion...')
    # get the JSON filelist
    try:
        jsonfilelist = os.listdir(gse_soft_dir)
    except:
        print("Error finding JSON files. Returning...")
        return
    if msrap_prepare:
        if jsonfilelist:
            msrap_prepare_json(gsm_json_filelist = jsonfilelist)
        else:
            print('No GSM JSON files detected. Skipping MSRA-p preparation...')
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

import preprocess
import os
import re

idir = os.path.join('recount-methylation-files','idats')
sdir = os.path.join('recount-methylation-files','gse_soft')

pfd = preprocess.get_mongofiles_for_preprocessing(idatsdir=idir,
softdir = sdir)

# filter on valid files


pfd['idats'][list(pfd['idats'].keys())[0]]

"""
