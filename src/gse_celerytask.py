#!/usr/bin/env python3

""" gse_celerytask.py
    Task script for celery.
    
    This script defines a task for celerybeat, where individual task managed by 
    celeryd and worker.

    Notes:
        * For best results, ensure RabbitMQ and celery are both running
"""

import celery
from celery import Celery
import os
import sys
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from serverutilities import gettime_ntp, get_queryfilt_dict
from dl import soft_mongo_date, idat_mongo_date, dl_idat, dl_soft
from update_rmdb import update_rmdb

app = Celery()
app.config_from_object('celeryconfig')

@app.task
def gse_task(gse_id, gsefiltdict = get_queryfilt_dict(), timestamp = gettime_ntp()):
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
    if gsefiltdict:
        print('File gsefiltdict provided, continuing...')
        gsmlist = gsefiltdict[gse_id]
        print("Detected N = "+str(len(gsmlist))+' GSM IDs...')
        rl.append(gsmlist)
        print("Beginning soft file download...")
        ddsoft = dl_soft(gse_list=[gse_id], timestamp = run_timestamp)
        rl.append(ddsoft)
        print('Beginning idat download...')
        ddidat = dl_idat(input_list = gsmlist, timestamp = run_timestamp)
        rl.append(ddidat)
        print('updating rmdb...')
        updateobj = update_rmdb(ddidat = ddidat, ddsoft = ddsoft)
        rl.append(updateobj)
        print('Task completed! Returning...')
        return rl
    else:
        print("Error: no gse query filt file provided. Returning...")
        return 0 

""" Notes and Tutorial

# testing server.py
import subprocess
import os
from gse_celerytask import gse_task

# start rabbitmq and celery
# brew services start rabbitmq # start the broker
# celery worker -A gse_tasks -l INFO

# run the queue
qstatlist = []
gselist = somelist
for gse in gse_list:
    qstatlist.append(gse_task.delay(gse))

"""