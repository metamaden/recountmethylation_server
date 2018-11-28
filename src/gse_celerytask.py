import celery
from celery import Celery
import glob
import os
import sys
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import dl
import update_rmdb
from update_rmdb import update_rmdb
import edirect_query
from edirect_query import querydict

""" This is the task script for celery
	This script defines a task to be executed by celerybeat
	tasks managed by celeryd and worker
	Notes:
		* RabbitMQ-server running (>rabbitmq-server)
		* Celery installed and running
		* Celery defaults to using RabbitMQ broker
"""

app = Celery()
app.config_from_object('config')

@app.task
def gse_task(gse_id, target='equery'):
    """ Define a GSE-based task for celery queue
        Arguments
            * gse_id : a single valid GSE id
            * target : dir or path to 'gsequery_filt.*' file
        Returns
            * rl (list) : list of download dictionaries and rmdb update statuses
    """
    rl = []
    # grab the latest gsequery_filt.* file from target
    gsefiltfiles = glob.glob('.'.join([os.path.join(target, 'gsequery_filt'),
        '*'])
    )
    if gsefiltfiles:
        if len(gsefiltfiles)>1:
            gsefiltfiles.sort(key=lambda x: int(x.split('.')[1]))
            gsefiltfile_mostrecent = gsefiltfiles[-1]
        else:
            gsefiltfile_mostrecent = gsefiltfiles[0]
        # load gse query filt object, get filtered gsm list
        gsefilt = querydict(query=os.path.join(target,gsefiltfile_mostrecent))
        gsmlist = gsefilt[gse_id]
        rl.append(gsmlist)
        ddsoft = dl_soft(gse_list=gse_id)
        rl.append(ddsoft)
        ddidat = dl_idat(input_list=gsmlist)
        rl.append(ddidat)
        updateobj = update_rmdb(ddidat=ddidat,ddsoft=ddsoft)
        rl.append(updateobj)
        return rl
    else:
        print("Error: no gsequery_filt.* files found in target. Returning...")
        return

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