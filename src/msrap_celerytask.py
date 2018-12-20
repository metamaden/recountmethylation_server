#!/usr/bin/env python3

import celery
from celery import Celery
import os
import sys
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from utilities import gettime_ntp, get_queryfilt_dict, getlatest_filepath
from process_soft import run_metasrapipeline

app = Celery()

@app.task
def msrap_task(gseid, gsm_jsonfnlist, gsefiltd = get_queryfilt_dict(),
    filesdir = 'recount-methylation-files', timestamp = gettime_ntp(),
    gsmjsondir = 'gsm_json'):
    gsmlist = gsefiltd[gseid]
    gsmlist = [gsmid for gsmid in gsmlist if gsmid in gsmfn_jsonlist]
    for gsmid in gsmlist:
        # select the latest file to run
        jsonfpath = getlatest_filepath(filepath = os.path.join(filesdir,gsmjsondir),
            filestr=gsmid,embeddedpattern=True)
        jsonfn = os.path.basename(jsonfpath)
        if jsonfn and not jsonfn==0:
            run_metasrapipeline(json_flist=[jsonfn])




    