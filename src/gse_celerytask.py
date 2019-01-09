#!/usr/bin/env python3

""" gse_celerytask.py
    Task script for celery. This script defines a task for celerybeat, where 
    individual task managed by celeryd and worker. Configuration is also 
    accessed from celeryconfig.py, including SQLite backend database info.
    Notes:
        * Broker: For best results, ensure RabbitMQ broker and celery are both 
            running in background.
    Functions:
        * gse_task: Job task definition for celery queue. Returns sparse info. 
            and GSE ID, for access from backend db.
"""

import celery
from celery import Celery
import os
import sys
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from utilities import gettime_ntp, get_queryfilt_dict
from dl import soft_mongo_date, idat_mongo_date, dl_idat, dl_soft
from update_rmdb import update_rmdb
import settings
settings.init()

app = Celery()
app.config_from_object('celeryconfig')

@app.task
def gse_task(gse_id, gsefiltdict = get_queryfilt_dict(), 
    timestamp = gettime_ntp()):
    """ GSE based task for celery job queue.
        Arguments
            * gse_id (str) : A single valid GSE id
            * gsefiltdict (dict) : GSE filtered query object, as dictionary read
                using querydict().
            * timestamp (str) : NTP timestamp for versioning file downloads.
        Returns
            * rl (list) : List of download dictionaries and rmdb update 
                statuses.
    """
    if not timestamp:
        run_timestamp = gettime_ntp()
    else:
        run_timestamp = timestamp
    print('Beginning GSE task, id: '+gse_id)
    rl = []
    rl.append(gse_id)
    if gsefiltdict:
        print('File gsefiltdict provided, continuing...')
        gsmlist = gsefiltdict[gse_id]
        print("Detected N = "+str(len(gsmlist))+' GSM IDs...')
        rl.append(True)
        print("Beginning soft file download...")
        ddsoft = dl_soft(gse_list=[gse_id], timestamp=run_timestamp)
        rl.append(True)
        print('Beginning idat download...')
        ddidat = dl_idat(input_list=gsmlist, timestamp=run_timestamp)
        rl.append(True)
        print('updating rmdb...')
        updateobj = update_rmdb(ddidat=ddidat, ddsoft=ddsoft)
        rl.append(True)
        print('Task completed! Returning...')
        return rl
    else:
        print("Error: no gse query filt file provided. Returning...")
        rl.append(None)
        return rl
