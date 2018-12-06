import celery
from celery import Celery
import glob
import os
import sys
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from dl import gettime_ntp, soft_mongo_date, idat_mongo_date, dl_idat, dl_soft
from update_rmdb import update_rmdb

""" gse_celerytask.py
    Task script for celery
	
    This script defines a task for celerybeat, where individual task managed by 
    celeryd and worker.

	Notes:
        * Check that RabbitMQ and celery are both running
"""

app = Celery()
app.config_from_object('config')

@app.task
def gse_task(gse_id, gsefiltdict):
    """ Define a GSE-based task for celery queue
        Arguments
            * gse_id : a single valid GSE id
            * gsefiltdict : gse filtered query as dictionary (see 
                'edirect_query.querydict()' function)
        Returns
            * rl (list) : list of download dictionaries and rmdb update statuses
    """
    print('Beginning gse_task for id: '+gse_id)
    rl = []
    if gsefiltdict:
        print('gsefiltdict provided, continuing...')
        # get gsms to pass to dl_idats
        gsmlist = gsefiltdict[gse_id]
        print('n = '+str(len(gsmlist))+' gsm ids detected...')
        # check for valid gsm ids here...
        rl.append(gsmlist)
        print('beginning soft file download...')
        ddsoft = dl_soft(gse_list=[gse_id])
        rl.append(ddsoft)
        print('beginning idat download...')
        ddidat = dl_idat(input_list=gsmlist)
        rl.append(ddidat)
        print('updating rmdb...')
        updateobj = update_rmdb(ddidat=ddidat,ddsoft=ddsoft)
        rl.append(updateobj)
        print('Task completed! Returning...')
        return rl
    else:
        print("Error: no gse query filt file provided. Returning...")
        return 0 

""" EXAMPLES

# testing server.py
import subprocess
import os
from gse_celerytask import gse_task

# start rabbitmq and celery
# brew services start rabbitmq # start the broker
# celery worker -A gse_tasks -l INFO

args1 = ["brew","services","start","rabbitmq"] # for macOSX
args2 = ["celery","worker","-A","gse_celerytask","-l","INFO"]
subprocess.check_command([args1],shell=False)
subprocess.check_command([args2],shell=False)

# on first run
# get new gse list from filtered equery

# on run >1st
# get new gse list from stored filtered equery file

# run the queue
qstatlist = []
gselist = somelist
for gse in gse_list:
    qstatlist.append(gse_task.delay(gse))

# 

"""