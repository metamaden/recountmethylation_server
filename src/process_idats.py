#!/usr/bin/env python3

""" process_idats.py

    Authors: Sean Maden, Abhi Nellore
    
    Functions to preprocess idats before being read into minfi.
    Functions:
        * expand_idats: Expand compressed idat files.
"""

import os
import sys
import re
import gzip
import shutil
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
import settings
settings.init()

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

def cleanup_idats(gsmfpathdict):
    """
    """
    eqd = get_queryfilt_dict()
    gsmvalidlist = list(set([gsmid for gselist in list(eqd.values()) 
        for gsmid in gselist
    ]))
    gsmvalid_fpathlist = {key:value for (key,value) in gsmfpathdict.items() 
        if key in gsmvalidlist
    }
    for gsmindex, gsmid in enumerate(gsmvalid_fpathlist, 1):
        if not gsmid in gsmvalidlist:
            
qstr = "esearch -db gds -query 'GPL21145[ACCN] and idat[suppFile] and gsm[ETYP]' | efetch -format docsum | xtract -pattern DocumentSummary -element Id Accession > gsmid"
output=subprocess.check_output(qstr, shell=True)

qstr = "esearch -db gds -query 'GPL21145 [ACCN] and idat [suppFile] and gsm [ETYP]'"
output=subprocess.check_output(qstr, shell=True)

dldict['gsmquery'].append(dlfilename)
    subp_strlist1 = ["esearch","-db","gds","-query",
    "'"+settings.platformid+"[ACCN] AND idat[suppFile] AND gsm[ETYP]'"
    ]
    subp_strlist2 = ["efetch","-format","docsum"]
    subp_strlist3 = ["xtract","-pattern","DocumentSummary",
        "-element","Id Accession",">",
        os.path.join(temp_make,dlfilename)
        ]
    args = " | ".join([" ".join(subp_strlist1),
        " ".join(subp_strlist2),
        " ".join(subp_strlist3)])



for item in idatfiles:
    if not item.split(".")[0] in of:
        print("File not valid. removing " + item)
        os.remove(os.path.join("recount-methylation-files", "idats", item))
    else:
        print(item + " has valid gsm id. Continuing...")



