#!/usr/bin/env python3
import pymongo
import sys
import os
import datetime
import inspect

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
    gsmindex = list(set([record['gsmid'] for record in idatslist]))
    # grab unique gse ids
    softlist = list(client.recount_methylation.gse.soft.find())
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
    softrecordsfilt = []
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

