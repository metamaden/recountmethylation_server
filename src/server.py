#!/usr/bin/env python3

""" server.py

Authors: Sean Maden, Abhinav Nellore

Description:
    Server script to manage an instance of the recount-methylation database.

Overview:
    A recount-methylation instance consists of files (namely edirect query 
    results, experiment metadata in soft format, and methylation array intensity 
    data or 'idat' files) obtained from edirect queries and ftp-called downloads 
    from the Gene Expression Omnibus (GEO). RMDB is a recount-methylation Mongo 
    database that aggregates file metadata as documents, including experiement 
    (GSE) and sample (GSM) ids, ftp addresses and file paths to downloaded 
    files, and a datetime-formatted date corresponding to last file update. 
    Files are versioned using NTP timestamps in filenames. 

    For best results, we recommend users attempt an initial setup of their 
    recount-methylation instance using default generated directory trees and 
    filenames, and do not directly change locations or names of files initially 
    downloaded.

Server Processes:
    The server.py script manages process queues, error handling, and 
    coordination of recount-methylation. It currently uses Celery distributed 
    task queue to queue jobs synchronously. Jobs are brokered using RabbitMQ, 
    and queue details are backed up locally in a SQLite db. It is recommend you 
    consult the SQLite backend database for details about interruptions to 
    server operations.

Dependencies and Setup:
    1. Recount-methylation primarily uses Python 3 for download handling and 
        file management. R is used for SOFT-to-JSON conversion, and for 
        preprocessing arrays. MetaSRA-pipeline, which runs using Python 2, is 
        used for mapping experiment metadata to ENCODE ontology terms, and we 
        recommend installing a fork of the original repo (available here: 
        <https://github.com/metamaden/MetaSRA-pipeline>). 
    2. Clone the recount-methylation-server repo from GitHub (available here: 
        <>).
    3. To run recount-methylation server.py, follow all provided setup and 
        readme instructions. Also ensure the following resources are installed 
        and running:
            * Celery (http://www.celeryproject.org/)
            * RabbitMQ (https://www.rabbitmq.com/)
            * MongoDB (https://www.mongodb.com/)
            * SQLite (https://www.sqlite.org/)
"""

import subprocess
import glob
import sys
import os
import re
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
import edirect_query
from edirect_query import gsm_query, gse_query, gsequery_filter  
from utilities import gettime_ntp, getlatest_filepath, querydict
from utilities import get_queryfilt_dict  
import settings
settings.init()

def firsttime_run(filedir='recount-methylation-files', 
    run_timestamp=gettime_ntp()):
    """ firsttime_run
        On first setup, run new equeries and query filter.
        Arguments:
            * filedir (str): Dir name for db files. 
            * run_timestamp (str) : NTP timestamp or function to retrieve it.
        Returns:
            * gseidlist (list): List of valid GSE IDs.
    """
    print("Beginning first time server run...")
    equery_dest = settings.equerypath
    temppath = settings.temppath
    gse_query()
    gsm_query()
    gseqfile = getlatest_filepath(equery_dest,'gse_edirectquery') 
    gsmqfile = getlatest_filepath(equery_dest,'gsm_edirectquery')
    gsequery_filter()
    gsefiltpath = getlatest_filepath(equery_dest,'gsequery_filt')
    if gsefiltpath:
        gsefiltd = querydict(querypath=gsefiltpath,splitdelim=' ')
        gseidlist = list(gsefiltd.keys())
        print("GSE id list of len "+str(len(gseidlist))+" found. Returning...")
        return gseidlist
    else:
        print("Error retrieving gse query filtered file. Returning...")
        return None

def scheduled_run(eqfilt_path=False, run_timestamp=gettime_ntp()):
    """ scheduled_run
        Tasks performed on regular schedule, after first setup. For the job 
        queue, a list of GSE ids is returned. The id list is filtered on 
        existing GSE soft files to prioritize unrepresented experiments for 
        download. 
        Arguments:
            * eqfilt_path (str) : Filepath to edirect query filter file.
            * filedir (str) : Root name of files directory.
            * run_timestamp (str) : NTP timestamp or function to retrieve it.
        Returns:
            * gse_list (list) : list of valid GSE IDs, or None if error occurs 

    """
    try:
        gsefiltd = get_queryfilt_dict()
    except:
        print("No gse query filt file found, checking for gse and gsm "
            +"queries...")
        gsequery_latest = getlatest_filepath(filepath=eqpath,
            filestr='gse_edirectquery')
        if not gsequery_latest:
            gse_query()
        gsmquery_latest = getlatest_filepath(eqpath,'gsm_edirectquery')
        if not gsmquery_latest:
            gsm_query()
        print("Running filter on gse query...")
        gsequery_filter()
        gsefiltd = get_queryfilt_dict()
    # get list of gse ids from existing soft files
    gsesoftfiles = os.listdir(settings.gsesoftpath)
    print(str(gsesoftfiles))
    rxgse = re.compile('GSE[0-9]*')
    gseid_softexists = [str(rxgse.findall(softfn)[0]) for softfn in gsesoftfiles
        if rxgse.findall(softfn)
    ]
    if gsefiltd:
        gseid_listall = list(gsefiltd.keys())
        print("GSE id list of len "+str(len(gseid_listall)) + " found. Filtering..")
        if gseid_softexists and len(gseid_softexists)>0:
            gseid_filt = [gseid for gseid in gseid_listall
                if not gseid in gseid_softexists
            ]
        else:
            gseid_filt = gseid_listall
        print("After filtering on existing soft files, N = "+str(len(gseid_filt))
            +" GSE ids remain. Returning id list...")
        # if all gse ids represented, return all gse ids for brand new run
        if len(gseid_filt)==len(gseid_listall):
            gseid_filt = gseid_listall
        return gseid_filt
    else: 
        print("Error forming equery filt dictionary. Returning...")
        return None

if __name__ == "__main__":
    """ Recount-methylation sever server.py main
        Code addresses various contingencies precluding generation of GSE ID 
        list. Once list can be made, it is used to populate a new Celery queue.
    """
    print("Starting server.py...")
    from gse_celerytask import gse_task
    import argparse
    from random import shuffle
    gselist = [] # queue input, gse-based
    qstatlist = [] # job status object, also stored at sqlite db
    print("Getting timestamp...")
    run_timestamp = gettime_ntp() # pass this result to child functions
    
    # If GSE ID provided, immediately parse it for download
    parser = argparse.ArgumentParser(description='Arguments for server.py')
    parser.add_argument("--gseid", type=str, required=False,
        default=None, help='Option to enter valid GSE ID for immediate download.')
    args = parser.parse_args()
    # For the job queue, either from provided argument or automation
    if args.gseid:
        print("Provided GSE id detected. Processing...")
        gqd = get_queryfilt_dict()
        qstatlist.append(gse_task(gse_id = args.gseid, gsefiltdict=gqd, 
                timestamp = run_timestamp
                )
            )
    else:
        print("No GSE id(s) provided. Forming id list for job queue...")
        files_dir = settings.filesdir
        if os.path.exists(files_dir):
            print("Directory : "+files_dir+" found. Running scheduled_run...")
            gselist = scheduled_run(run_timestamp=run_timestamp)
        else:    
            print("Directory : "+files_dir+" not found. Creating filesdir and "
                +"running firsttime_run...")
            os.makedirs(files_dir, exist_ok=True)
            gselist = firsttime_run(run_timestamp=run_timestamp)
        if gselist:
            print("Shuffling GSE id list...")
            shuffle(gselist) # randomize GSE ID order
            print("Beginning job queue for GSE id list of "+str(len(gselist))
                +" samples...")
            gqd = get_queryfilt_dict() # one eqfilt call for all jobs this run
            for gse in gselist:
                qstatlist.append(gse_task(gse_id=gse, gsefiltdict=gqd, 
                    timestamp=run_timestamp
                    )
                )    
        else:
            print("Error: Valid gselist absent. Returning...")
            
