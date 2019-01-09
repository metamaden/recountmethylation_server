#!/usr/bin/env python3

""" update_rmdb.py
    Functions to update Recount Methylation MongoDB (or "RMDB"), a database of 
    file metadata for recount methylation. 
    Notes:
        * Job: This function is run as part of the GSE-based job definition for
            celery job queue. After attempting idat and soft file downloads,
            any newly detected files have their metadata stored as docs in RMDB.
    Functions:
        * update_rmdb: Update RMDB with any metadata for newly downloaded files.
"""

import pymongo
import datetime
import os
import sys
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import settings
settings.init()

def update_rmdb(ddidat, ddsoft, host=settings.rmdbhost, port=settings.rmdbport):
    """ Update recount-methylation database compilations with new documents.
        Arguments
            * ddidat : Download dictionary from dl_idats, as returned by 
                dl_idats().
            * ddsoft : Download dicitonary from dl_soft, as returned by 
                dl_soft(). 
        Returns
            * statusdict object (dictionary): Result list (1 = new doc added, 
                0 = no doc added) .
    """
    statusdict = {}
    client = pymongo.MongoClient('localhost', 27017)
    rmdb = client.recount_methylation
    if ddidat:
        lvals = []
        statusdict['ddidat'] = []
        gsmc = rmdb.gsm
        idatsc = gsmc.idats
        for gsmkey in list(ddidat.keys()):
            lvals = ddidat[gsmkey]
            for lval in lvals:
                if lval[-1]==True:
                    new_idatdoc = {"gsmid":lval[0],
                        "ftpaddress":lval[1],
                        "filepath":lval[2],
                        "exitstatus":lval[3],
                        "date":lval[4]
                        }
                    idatsc.insert_one(new_idatdoc)
                    statusdict['ddidat'].append(1)
                else:
                    statusdict['ddidat'].append(0)
    if ddsoft:
        svals = []
        statusdict['ddsoft'] = []
        gsec = rmdb.gse
        softc = gsec.soft
        for gsekey in list(ddsoft.keys()):
            svals = ddsoft[gsekey]
            if svals[-1]==True:
                ssoft = svals[1]
                new_softdoc = {"gseid":gsekey,
                        "ftpaddress":ssoft[1],
                        "filepath":ssoft[2],
                        "exitstatus":ssoft[3],
                        "date":ssoft[4]
                        }
                softc.insert_one(new_softdoc)
                statusdict['ddsoft'].append(1)
            else:
                statusdict['ddsoft'].append(0)
    return statusdict
