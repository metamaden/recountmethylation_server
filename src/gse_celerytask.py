#!/usr/bin/env python3

""" gse_celerytask.py
    
    Authors: Sean Maden, Abhi Nellore
    
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

import celery, os, sys; from celery import Celery
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, get_queryfilt_dict
from dl import soft_mongo_date, idat_mongo_date, dl_idat, dl_soft
from update_rmdb import update_rmdb
import settings; settings.init()

app = Celery(); app.config_from_object('celeryconfig')

@app.task
def gse_task(gse_id, gsefiltdict = get_queryfilt_dict(), 
    timestamp = gettime_ntp()):
    """ GSE based task for celery job queue.
        
        Arguments
            * gse_id : A single valid GSE id (str).
            * gsefiltdict : GSE filtered query object, as dictionary read
                using querydict() (dict).
            * timestamp : NTP timestamp for versioning file downloads (str).
            
        Returns
            * rl, a list of download dictionaries and rmdb update statuses.
            
    """
    if not timestamp:
        run_timestamp = gettime_ntp()
    else:
        run_timestamp = timestamp
    print('Beginning GSE task, ID: '+gse_id); rl = []; rl.append(gse_id)
    if gsefiltdict:
        print('File gsefiltdict provided, continuing...')
        gsmlist = gsefiltdict[gse_id]
        print('Detected N = '+str(len(gsmlist))+' GSM IDs...')
        if len(gsmlist) > 0:
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
        else:
            print('No valid GSM IDs detected for study GSE ID ', gse_id, 
                  ', skipping...')
            rl.append(None)
        print('Task completed! Returning...')
        return rl
    else:
        print("Error: no GSE query filt file provided. Returning...")
        rl.append(None)
        return rl
