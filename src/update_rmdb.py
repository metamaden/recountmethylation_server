#!/usr/bin/env python3

import pymongo
import datetime

def update_rmdb(ddidat,ddsoft,host='localhost',port=27017):
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

"""

# clear the db
import pymongo 
client = pymongo.MongoClient('localhost', 27017)
rmdb = client.recount_methylation
list(rmdb.list_collections())
rmdb.drop_collection('gsm.idats')
rmdb.drop_collection('gse.soft')

# Sample update

import datetime
import sys
sys.path.insert(0,'/usr/local/lib/python3.7/site-packages')
import pymongo
from mongodb_test import update_rmdb

# from dl_idats
ddidat = {'GSM1505330': [['GSM1505330', 'geo/samples/GSM1505nnn/GSM1505330/suppl/', 
'connection success, valid num idats found'], 
['GSM1505330', True, 
'temp_dir/tmp50wvd8sj/GSM1505330.timestamp.GSM1505330_9376538060_R04C01_Grn.idat.gz', 
'226 Transfer complete', datetime.datetime(2014, 9, 16, 15, 15, 45), 'new_date'], 
['GSM1505330', True, 
'temp_dir/tmp50wvd8sj/GSM1505330.timestamp.GSM1505330_9376538060_R04C01_Red.idat.gz', 
'226 Transfer complete', datetime.datetime(2014, 9, 16, 15, 15, 45), 'new_date']]}

gsmkeys = list(ddidat.keys())
gsmkey = gsmkeys[0]
lvals = ddidat[gsmkey]
lval = lvals[1]
lval[1]==True


ddsoft = {'GSE109904': [
['GSE109904', 'geo/series/GSE109nnn/GSE109904/soft/', 'success'], 
['GSE109904', 'geo/series/GSE109nnn/GSE109904/soft/GSE109904_family.soft.gz', 
    'temp_dir/tmpzgz2ma1t/GSE109904.1543256167.GSE109904_family.soft.gz', 
    '226 Transfer complete', datetime.datetime(2018, 11, 21, 4, 29, 57), 
    'new_date'], 
    False]}

update_rmdb(ddidat=ddidat,ddsoft=ddsoft)

client = pymongo.MongoClient('localhost', 27017)
rmdb = client.recount_methylation
gsmc = rmdb.gsm
idatsc = gsmc.idats
list(idatsc.find())

gsec = rmdb.gse
softc = gsec.soft
list(softc.find())

softc.find({'gse' : gse},{'date' : 1})
softc.find({'gse' : gse},{'date' : 1})

"""