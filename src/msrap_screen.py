#!/usr/bin/env python3

import re
import subprocess
import argparse
import sys
import os
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from process_soft import run_metasrapipeline
from utilities import gettime_ntp

""" msrap_screen.py
    Manages n screen calls for MetaSRA-pipeline, running the pipeline on
    subsets of the JSON files indicated
    Arguments
        * nscreen (int): number of screens in total
    Returns
        * null, generating n screen sessions for MetaSRA-pipeline as a side-effect.
"""


def list_all_gsms_toprocess(filesdir="recount-methylation-files",
    gsmjsondir="gsm_json",msrapdir="gsm_msrap_outfiles", jsonstr="json",
    msrapstr="msrapout",gsmstr="GSM"):
    """ Grab GSM IDs for which there are JSON files but not MSRAP files
        * returnfn (True/False, bool.): whether to return filenames
    """
    rxgsm = re.compile(".*"+gsmstr+".*")
    # get JSON GSM IDs
    rxjson = re.compile(".*"+jsonstr+"$")
    json_flist = os.listdir(os.path.join(filesdir,gsmjsondir))
    json_flist = list(filter(rxjson.match,jsonflist))
    json_gsmid = [jfn.split(".")[1] for jfn in json_flist]
    json_gsmid = list(filter(rxgsm.match,json_gsmid))
    # get MSRAP GSM IDs
    rxmsrap = re.compile(".*"+msrapstr+"$")
    msrap_flist = os.listdir(os.path.join(filesdir,msrapdir))
    msrap_flist = list(filter(rxjson.match,msrapflist))
    msrap_gsmid = [mfn.split(".")[1] for mfn in msrap_flist]
    msrap_gsmid = list(filter(rxgsm.match,msrap_gsmid))
    # get diffs list
    gsmid_toprocess = [gsmid for gsmid in json_gsmid
        if not gsmid in msrap_gsmid
    ]
    return gsmid_toprocess

def generate_gsm_fnsublist(startindex, endindex, idlist=list_all_gsms_toprocess(), 
    filesdir="recount-methylation-files", gsmjsondir="gsm_json", 
    jsonstr="json", gsmstr="GSM"):
    """ Generate a GSM files or ids sublist
        Arguments:
            * returnfn (True/False, bool.): whether to return filenames
        Returns:
            * gsm_fnlist (list) : list of gsm filenames from JSON dir 
    """
    try:
        idsublist = idlist[statindex:endindex]
        # get JSON GSM IDs
        rxjson = re.compile(".*"+jsonstr+"$")
        json_flist = os.listdir(os.path.join(filesdir,gsmjsondir))
        json_flist = list(filter(rxjson.match,jsonflist))
        json_filt_flist = []
        for gsmid in idsublist:
            rxgsmi = re.compile(".*"+gsmid+".*")
            gsmjsoni = list(filter(rxgsmi.match,json_flist))
            if len(gsmjsoni)>0:
                json_filt_flist.append(gsmjsoni)
        return return gsm_fnlist   
    except:
        print("Error getting sublist from provided indices and id list!")
        return 0

if __name__ == "__main__":
    parser.add_argument("--nscreens", type=int, required=True,
        default=5, help='Number of screens to deploy')
    args = parser.parse_args()
    # screen denom
    ns = args.nscreens
    # form list of gsm ids to process from json to msrap
    gsmid_list = get_gsmids_toprocess()
    # form indices
    istart = list(range(0,len(gsmid_list),ns))
    iend = istart[1::]
    msrap_screenscript_path = os.path.join('recount-methylation-server',
            'src','msrap_runscreen.py')
    for index, item in enumerate(iend): 
        cmdlist = ['screen',
            'python3',
            msrap_screenscript_path,
            '--startindex ', str(istart[index]),
            '--endindex ', str(iend[index])
            ]
        subprocess.call(cmdlist,shell=False)
    # final index range
    cmdlist = ['screen',
            'python3',
            msrap_screenscript_path,
            '--startindex ', str(istart[-1]),
            '--endindex ', str(len(gsmid_list)-1)
            ]

