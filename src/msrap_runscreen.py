import re
import subprocess
import argparse
import sys
import os
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from utilities import gettime_ntp

def run_msrap_screen(startindex,endindex, gsm_jsondir = 'gsm_json',
    filesdir = 'recount-methylation-files', msrap_fn = 'msrapout', 
    msrap_destdir = 'gsm_msrap_outfiles', msrap_dir = '.',
    timestamp = gettime_ntp(), gsmsoftindex = 1, gsmsoftindexdenom = 5,
    deployscreens = False):
    """ Run MetaSRA-pipeline on available GSM JSON files.
        
    """
    # get the gsm id list
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
    # get the fn list
    idsublist = gsmid_toprocess[startindex:endindex]
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
    gsm_json_fn_list = json_filt_flist
    # iterate over valid filenames at gsm json dir
    for gsm_json_fn in gsm_json_fn_list:
        outfn = os.path.splitext(gsm_json_fn)[0] # fn without extension
        gsmjson_readpath = os.path.join(gsm_jsonpath, gsm_json_fn)
        gsm_msrapout_writepath = os.path.join(msrap_destpath,
            ".".join([timestamp,outfn,msrap_fn]))
        try:
            cmdlist = ['python2',
                msrap_path,
                gsmjson_readpath,
                gsm_msrapout_writepath
                ]
            subprocess.call(cmdlist,shell=False)
            msrap_statlist.append(1)
        except subprocess.CalledProcessError as e:
            msrap_statlist.append(e)
    return msrap_statlist

if __name__ == "__main__":
    parser.add_argument("--startindex", type=int, required=True,
        default=1, help='')
    parser.add_argument("--endindex", type=int, required=True,
        default=50, help='')
    args = parser.parse_args()
    run_msrap_screen(args.startindex,args.endindex)






