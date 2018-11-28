#!/usr/bin/env python3

import subprocess
import glob
import sys
import os
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import edirect_query
from edirect_query import gsm_query, gse_query, querydict, gsequery_filter
from gse_celerytask import gse_task

def getlatest_filepath(filepath,filestr):
    """ Get path the latest version of a file, based on timestamp
        Arguments
            * filepath (str) : path to dir to search
            * filestr (str) : pattern of file to search
        Returns
            * latest_file_path (str) : path to latest version of file, OR
            * 0 : search turned up no files at location
    """
    filelist = glob.glob('.'.join([os.path.join(filepath, filestr), '*']))
    if filelist:
        if len(filelist) > 1:
            # sort on timestamp
            filelist.sort(key=lambda x: int(x.split('.')[1]))
            latest_file_path = filelist[-1]
        else:
            latest_file_path = filelist[0]
        return latest_file_path
    else:
        return 0       

def firsttime_run(filedir = 'recount-methylation-files'):
    """ Tasks performed on first time setup
        Arguments
            * filedir (str): dir name for db files 
        Returns
            * gseidlist (list): list of valid GSE IDs
    """
    equery_dest = os.path.join(filedir,'equery')
    temp = os.path.join(filedir,'temp')
    # run equeries 
    gse_query(dest = equery_dest, temp = temp)
    gsm_query(dest = equery_dest, temp = temp)
    # run filter
    gsequeryfile = glob.glob(os.path.join(equery_dest,'gse_edirectquery.*'))
    gsmqueryfile = glob.glob(os.path.join(equery_dest,'gsm_edirectquery.*'))
    gsequery_filter(gsequery = gsequeryfile,
        gsmquery = gsmqueryfile,
        target = equery_dest,
        splitdelim = ' '
        )
    gsefiltpath = glob.glob(os.path.join(equery_dest,'gsequery_filt.*'))
    gsefiltd = querydict(query = gsefiltpath)
    gseidlist = list(gsefilt.keys())
    print("GSE id list of len "+str(len(gseidlist))+" found. Returning...")
    return gseidlist

def scheduled_run(eqfilt_path=False, filedir = 'recount-methylation-files'):
    """ Tasks performed on regular schedule, after first setup
        Arguments
            * filedir : file to search for filtered equery files
            * eqfilt_path : path to GSE equery filtered file
        Returns
            * 0 (int) : If error encountered, or
            * gse_list (list) : list of valid GSE IDs
    """
    if eqfilt_path:
        gsefiltd = querydict(query=eqfilt_path,splitdelim=' ')
        if gsefiltd:
            gseidlist = list(gsefiltd.keys())
            return gseidlist
        else:
            print("Error processing file at provided path: "+eqfilt_path)
            return 0     
    else:
        print("No gse filt file dir provided, checking default location...")
        # check if gse query filt file exists
        eqpath = os.path.join(filedir,'equery')
        gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt')
        if not gsefilt_latest or gsefilt_latest == 0:
            print("No gse query filt file found, checking for gse and gsm "
                +"queries...")
            # run fresh gse/gsm queries if not found
            gsequery_latest = getlatest_filepath(eqpath,'gse_edirectquery')
            gsmquery_latest = getlatest_filepath(eqpath,'gsm_edirectquery')
            if not gsequery_latest or gsequery_latest == 0:
                gse_query(dest = eqpath,
                    temp=os.path.join(filedir,'temp')
                    )
            if not gsmquery_latest or gsmquery_latest == 0:
                gsm_query(dest = eqpath,
                    temp = os.path.join(filedir,'temp')
                    )
            # run new gse filt
            print("Running filter on gse query...")
            gsequery_filter(gsequery=gsequery_latest,
                    gsmquery = gsmquery_latest,
                    target = eqpath, splitdelim = ' '
                    )
            gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt')
        else:
            print("Gse filt file found. Continuing...")
    if gsefilt_latest:
        print("Getting dictionary from gse filt file...")
        gsefiltd = querydict(query = gsefilt_latest)
        gseidlist = list(gsefiltd.keys())
        # print(gseidlist)
        print("GSE id list of len "+str(len(gseidlist)) + " found. Returning..")
        return gseidlist
    else: 
        print("Error, attempt to produce gseidlist failed. Check for a valid "
            +"gse query filt file and gse/gsm query files.")
        return 0

def serverstart(rabbitmq_args, celery_args, mongodb_args, platform = "mac"):
    """ What to run on server start (setup dependencies in console)
        Arguments
            * rabbitmq_args (list) : subprocess commands to start RabbitMQ
            * celery_args (list) : subprocess commands to start Celery
            * mongodb_args (list) : subprocess commands to start MongoDB
            * platform (str) : Either 'mac', 'ux', or 'pc'. 
                If provided, overrides given lists with new lists
        Returns
            * rd (dict) : dictionary of subprocess objects
    """
    if platform:
        if platform=="mac":
            rabbitmq_args = ["brew","services","start","rabbitmq"]
            celery_args = ["celery","worker","-A","gse_celerytask","-l","INFO"]
            mongodb_args = ["mongod"]
        if platform=="nix":
            rabbitmq_args = []
            celery_args = []
            mongodb_args = []
        if platform=="pc":
            rabbitmq_args = []
            celery_args = []
            mongodb_args = []
        else: 
            print("invalid platform provided. Returning...")
            return 0
    # define process running (check active processes if possible)
    if not processrunning:
        celery_run = 0
        try: 
            celery_start = subprocess.check_call(celery_args, shell=False)
        except subprocess.CalledProcessError as e:
             celery_start = e
    else:
        celery_run = 1
        celery_start = 0
    # rabbitmq
    if not processrunning:
        rabbitmq_run = 0
        try:
            rabbitmq_start = subprocess.check_call(rabbitmq_args, shell=False)
        except subprocess.CalledProcessError as e:
            rabbitmq_start = e
    else:
        rabbitmq_run = 1
        rabbitmq_start = 0
    # mongodb
    if not processrunning:
        mongodb_run = 0
        try:
            mongodb_start = subprocess.check_call(mongodb_args, shell=False)
        except subprocess.CalledProcessError as e:
            mongodb_start = e
    else:
        mongodb_run = 1
        mongodb_start = 0
    rd = {'rabbitmq': [rabbitmq_run, rabbitmq_start],
            'celery': [rabbitmq_run, celery_start],
            'mongodb': [mongodb_run, mongodb_start]
            }
    return rd

def run_gsequeue(gse_list):
    """ Run task queue with GSE IDs
        Arguments
            * gse_list
        Returns
            * 
    """
    qstatlist = []
    # run tasks synchronously in celery
    for gse in gse_list:
        qstatlist.append(gse_task.apply_async(gse))
    return qstatlist

def on_restart():
    """ Handling server.py restart, inc. interruptions when processing queue
        Arguments
        Returns
    """

def main(files_dir='recount-methylation-files'):
    """ Script to run on call from cl.
        Arguments
        Returns
    """

if __name__ == "__main__":
    files_dir = 'recount-methylation-files'
    if os.path.exists(files_dir):
        print(files_dir+" found. Running scheduled_run...")
        gselist = scheduled_run()
        print(str(len(gselist)))
        print(str(gselist[0]))
        # queuerun = run_gsequeue(gse_list = scheduled_run())
        queuerun = run_gsequeue(gse_list = gselist[0])
    else:    
        print(files_dir+" not found. Creating dir and running firsttime_run...")
        os.makedirs(files_dir)
        queuerun = run_gsequeue(gse_list = firsttime_run())
        
        
