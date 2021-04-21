#!/usr/bin/env python3

"""

Author: Sean Maden

Description:

    Exclude a vector of GSM IDs from this instance.

Functions:

    * eqd_gsm_exclude
 
"""

import subprocess, os, socket, struct, sys, time, tempfile, atexit, shutil
import glob, filecmp, re; from itertools import chain
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, querydict, getlatest_filepath
import settings; settings.init(); from server import firsttime_run
from utilities import get_queryfilt_dict, querydict, gettime_ntp
from edirect_query import gse_query_diffs

def eqd_gsm_exclude(equery_dest=settings.equerypath, filesdir=settings.filesdir,
    gsmv_fname="gsmv.txt", exclude_dpath=os.path.join("inst", "freeze_gsmv")):
    """ Exclude GSM IDs from edirecty query objects

    Arguments:
        * gsmv_fname: Name of the file to load. Should include only 
            space-separated sample/GSM IDs in a single line.
        * exclude_dpath: Path to directory containing the file gsmv_fname.

    Returns:
        * Returns the path to the new filtered file at settings.equerypath.

    """
    gsmv_fpath = os.path.join(exclude_dpath, gsmv_fname)
    if not os.path.exists(gsmv_fpath):
        print("Couldn't find sample ID file")
    gsmv_exclude = [line.rstrip('\n').split(" ") 
                        for line in open(gsmv_fpath)][0]
    # gsmv_exclude = [i for sublist in gsmv_exclude for i in sublist]
    eqpath = settings.equerypath
    gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt', 
            embeddedpattern=True, tslocindex=1, returntype='returnlist'
        )[0]
    print("Starting with latest detected filter file: "+gsefilt_latest)
    querylines = [line.rstrip('\n') for line in open(gsefilt_latest)]
    qlnew = []; print("Applying filter..."); numgsm_old = len(querylines)
    for line in querylines:
        line = line.split(" ")
        ldat = [gid for gid in line if not gid in gsmv_exclude]
        numgsm_new = len(ldat)
        if len(ldat) > 1:
            qlnew.append(ldat)
    print("After filter, retained " + str(len(qlnew)) + " studies.")
    nts = gettime_ntp()
    newfpath = os.path.join(eqpath, ".".join(["gsequery_filt",nts]))
    print("Writing new filter file: ", newfpath)
    with open(newfpath, "w") as wf:
        for line in qlnew:
            wf.write(" ".join(line) + "\n")
        return newfpath

if __name__ == "__main__":
    """ gsm_exclude

    Exclude a file containing a space-sparated list of GSM IDs (a.k.a. gsmv 
    file) prior to synchronizing a recountmethylation instance.

    Arguments:

    * fname: Name of the file containing the GSM IDs for exclusion.
    * dpath: Path to the directory containing the file with GSM IDs for 
        exclusion.

    Returns:

    None, modifies the contents of files targeted by an instance.

    """
    import argparse, os
    parser = argparse.ArgumentParser(description='Arguments for gsm_exclude.py')
    parser.add_argument("--fname", type=str, required=False, default="gsmv.txt", 
        help='File containing space-separated GSM IDs to exclude from instance')
    parser.add_argument("--dpath", type=str, required=False, 
        default=os.path.join("inst", "freeze_gsmv"), 
        help='Path to directory containg file with GSM IDs for exclusion.')
    args = parser.parse_args()
    gsmv_fpath = os.path.join(args.dpath, args.fname)
    if os.path.exists(gsmv_fpath):
        eqd_gsm_exclude(gsmv_fname = args.fname, exclude_dpath = args.dpath)
    else:
        print("Error, couldn't find gsmv file at path: " + str(gsmv_fpath))
