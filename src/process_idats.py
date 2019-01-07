#!/usr/bin/env python3

import os
import sys
import re
import gzip
import shutil
# sys.path.insert(0, os.path.join("recount-methylation-server","src"))

""" process_idats.py
    Functions to preprocess idats before being read into minfi.
    Functions:
        * expand_idats: Expand compressed idat files.
"""

def expand_idats(filesdir = 'recount-methylation-files', idatsdir = 'idats'):
    """ Detect and expand available idat files.
        Arguments:
            * filesdir (str) : Root name of directory containing files.
            * idatsdir (str) : Name of directory containing idat files.
        Returns:
            * ridatd object (dictionary)
    """
    idatspath = os.path.join(filesdir,idatsdir)
    idats_fnlist = os.listdir(idatspath)
    rcompressed1 = re.compile(".*idat.gz$")
    idats_fnlist_filt = list(filter(rcompressed1.match, idats_fnlist)) 
    ridatd = {} # return dictionary
    for compidat in idats_fnlist_filt:
        idat_fn = os.path.splitext(compidat)[0]
        statuslist = []
        ridatd[compidat] = []
        with gzip.open(os.path.join(idatspath, compidat), 'rb') as f_in:
            with open(os.path.join(idatspath, idat_fn), 'wb') as f_out:
                try:
                    shutil.copyfileobj(f_in, f_out)
                    ridatd[compidat].append(1)
                except shutil.Error as se:
                    ridatd[compidat].append(se)
    return ridatd
