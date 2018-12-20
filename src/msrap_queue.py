#!/usr/bin/env python3
import os 
import sys
from msrap_celerytask import msrap_task
import argparse
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from utilities import gettime_ntp, get_queryfilt_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("--nchunk", type=int, required=True,
        default=5, help='Number of jobs/GSE IDs per queue chunk.')
    args = parser.parse_args()
    filesdir = 'recount-methylation-files'
    gsmjsondir = 'gsm_json'
    gsmjsonpath = os.path.join(filesdir,gsmjsondir)
    gsmjson_fileslist = os.listdir(gsmjsonpath)
    gsmjson_idlist = [os.splitext(gsmfn)[1] for gsmfn in gsmjson_fileslist]
    # get list of GSE IDs corresponding to any GSM JSON files
    gse_fgsm = [gseid for gseid in list(gsefiltdict.keys())
        if bool(set(gsefiltdict[gseid]) & set(gsmjson_idlist))
    ]
    if args.nchunk:
        jobs = msrap_task.chunks(gse_idlist, args.nchunk) # break the list into 30 chunks. Experiment with what number works best here.
        jobs.apply_async(gsm_jsonfnlist = gsmjson_idlist)
    else:
        print("Error determining nchunk! Returning...")


