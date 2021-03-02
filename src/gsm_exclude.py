#!/usr/bin/env python3

"""

Author: Sean Maden

Description:
    Exclude a vector of GSM IDs from this instance.

Functions:
    * eqd_gsm_exclude
 
"""

import subprocess, os, socket, struct, sys, time, tempfile, atexit, shutil
import glob, filecmp, re
from itertools import chain
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, querydict, getlatest_filepath
import settings
settings.init()
from server import firsttime_run
from utilities import get_queryfilt_dict, querydict, gettime_ntp
from edirect_query import gse_query_diffs

def eqd_gsm_exclude(gsmv_fname = "gsm_exclude", 
    exclude_dpath = os.path.join("."),  filesdir = settings.filesdir,
    equery_dest = settings.equerypath):
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
    gsmv_exclude = [line.rstrip('\n') for line in open(gsmv_fpath)]
    eqpath = settings.equerypath
    gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt', 
            embeddedpattern=True, tslocindex=1, returntype='returnlist'
        )[0]
    print("Starting with latest detected filter file: "+gsefilt_latest)
    querylines = [line.rstrip('\n') for line in open(gsefilt_latest)]
    qlnew = []
    print("Applying filter...")
    for line in querylines:
        ldat = line.split(" ")
        numgsm_old = len(ldat)
        ldat = [gid for gid in ldat 
                    if gid[0:3] == "GSE" or 
                    gid not in gsmv_exclude]
        numgsm_new = len(ldat)
        if len(ldat) > 1:
            qlnew.append(ldat)
    print("After filter, retained " + len(qlnew) + " studies.")
    nts = gettime_ntp()
    newfpath = os.path.join(eqpath, ".".join(["gsequery_filt",nts]))
    print("Writing new filter file: ", newfpath)
    with open(newfpath, "w") as wf:
        for line in qlnew:
            wf.write(" ".join(line) + "\n")
    return newfpath

if __name__ == "__main__":
    """
    """
    eqd_gsm_exclude()
