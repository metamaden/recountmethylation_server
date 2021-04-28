#!/usr/bin/env python3

""" rsheet.py

    Author: Sean Maden

    Description:
        Get data for a report on a recountmethylation instance.
    
    Functions:
        * idats_report: Get report stats for IDAT files.
        * soft_report: Get report stats for SOFT files.
        * equery_fractions: Get fractions of targeted GSM and GSE IDs currently
            represented among available instance files.

"""

import sys, os, re
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, get_queryfilt_dict
import settings; settings.init()

def idats_report(strmatchl = [".*idat.gz$", ".*idat$", ".*hlink.*"], 
    ipath = settings.idatspath):
    """ Get report stats for IDATs

    Arguments:

    Returns:
        * ddidat, dictionary of IDAT file stats.
    """
    eqd = get_queryfilt_dict()
    gsmv = list(set([v for sublist in eqd.values() for v in sublist]))
    idatsv = os.listdir(settings.idatspath)
    ddidat = {}
    ugsmv = list(set([i.split(".")[0] for i in idatsv]))
    ddidat["unique.gsmv"] = ugsmv
    ddidat["unique.gsm"] = len(ugsmv)
    fract_num = len(set([gsm for gsm in gsmv if gsm in ugsmv]))
    fract_denom = len(set(gsmv))
    ddidat["eqd.fract.gsm"] = fract_num/fract_denom
    for t in strmatchl:
        lm = list(filter(re.compile(t).match, idatsv))
        ddidat[t] = len(lm)
    return ddidat

def soft_report(ddidat, gsesoftpath = settings.gsesoftpath, 
    gsmsoftpath = settings.gsmsoftpath, 
    gsmjsonfiltpath = settings.gsmjsonfiltpath,
    strmatchl_gsesoft = ".*family.soft$", strmatchl_gsmsoft = ".*soft$",
    strmatchl_gsejsonfilt = ".*\\.json\\.filt$"):
    """ Get report stats for SOFT files

    Arguments:
        * ddidat: Results object returned by `idats_report()`.
        * gsesoftpath: Path to GSE SOFT files.
        * gsmsoftpath: Path to GSM SOFT files.
        * gsmjsonfiltpath: Path to filtered JSON files.
        * strmatchl_gsesoft: List of regex patterns for GSE SOFT files.
        * strmatchl_gsmsoft: List of regex patterns for GSM SOFT files.
        * strmatchl_gsejsonfilt: List of regex patterns for GSM filtered 
            JSON files.

    Returns:
        * ddsoft, dictionary of SOFT file stats for report.

    """
    eqd = get_queryfilt_dict()
    gsev = list(set([k for k in eqd.keys()]))
    gsmv = list(set([v for sublist in eqd.values() for v in sublist]))
    gsm_idatv = ddidat["unique.gsmv"] # gsm ids vector from idats_report()
    gsesoftv = list(set(os.listdir(gsesoftpath)))
    gsmsoftv = list(set(os.listdir(gsmsoftpath)))
    gsmjsonfiltv = list(set(os.listdir(gsmjsonfiltpath)))
    ddsoft = {}
    lm = list(filter(re.compile(strmatchl_gsesoft).match, gsesoftv))
    ddsoft["gsesoft " + strmatchl_gsesoft] = len(lm)
    lm = list(filter(re.compile(strmatchl_gsesoft).match, gsmsoftv))
    ddsoft["gsmsoft " + strmatchl_gsesoft] = len(lm)    
    lm = list(filter(re.compile(strmatchl_gsejsonfilt).match, gsmsoftv))
    ddsoft["gsmjsonfilt " + strmatchl_gsesoft] = len(lm)
    ddsoft["gsm_soft_and_idat"] = len(set([gsm for gsm in gsmsoftv
        if gsm in gsm_idatv]))
    gsesoftv_gsev = [fn.split(".")[0] for fn in gsesoftv]
    gsenum = len(list(set([gse for gse in gsev if gse in gsesoftv_gsev])))
    gsedenom = len(list(set(gsev)))
    ddsoft["eqd.fract.gse"] = gsenum/gsedenom
    gsmsoftv_gsmv = [fn.split(".")[0] for fn in gsmsoftv]
    gsmnum = len(list(set([gsmv for gsm in gsmv if gsm in gsmsoftv_gsmv])))
    gsedenom = len(list(set(gsmv)))
    ddsoft["eqd.fract.gsm.soft"] = gsmnum/gsedenom
    gsmboth_gsmv = [gsm for gsm in gsmsoftv_gsmv if gsm in gsm_idatv]
    gsmbothnum = len(list(set(gsmboth_gsmv)))
    gsmbothdenom = len(list(set(gsmv)))
    ddsoft["eqd.fract.gsm.both"] = gsmnum/gsedenom
    return ddsoft

if __name__ == "__main__":
    """ Report on recountmethylation instance files
    """
    ddidat = idats_report()
    ddsoft = soft_report()
    ddreport = {"ddidat" : ddidat, "ddsoft" : ddsoft}
    ddreport