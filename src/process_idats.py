#!/usr/bin/env python3

""" process_idats.py

    Authors: Sean Maden, Abhi Nellore
    
    Functions to preprocess idats before being read into minfi.
    
    Functions:
        * expand_idats: Expand compressed idat files.

"""

import os, sys, re, gzip, shutil; from random import shuffle
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
import settings; settings.init()

def expand_idats(idatspath = settings.idatspath, compext = ".*idat.gz$", 
    expext = ".*idat$"):
    """ expand_idats

        Detect and expand available idat files.
        
        Arguments:
        * idatspath : Path to instance directory containing downloaded 
            IDATs (valid file path).
        * compext : Regular expression pattern for extension of compressed
            IDAT files (string, regex pattern).
        * expext : Regular expression pattern for extension of expanded 
            IDAT files (string, regex pattern).

        Returns:
        * ridatd dictionary containing expanded IDAT info.
    """
    idats_fnlist = os.listdir(idatspath)
    rexpanded1 = re.compile(expext); rcompressed1 = re.compile(compext)
    idats_fnlist_filt1 = list(filter(rexpanded1.match, idats_fnlist))
    idats_fnlist_filt1 = [fn.split(".")[0] for fn in idats_fnlist_filt1]
    idats_fnlist_filt2 = list(filter(rcompressed1.match, idats_fnlist))
    idats_fnlist_filt = [fn for fn in idats_fnlist_filt2 
                            if not fn.split(".")[0] in idats_fnlist_filt1]
    shuffle(idats_fnlist_filt); ridatd = {}; nfile = 1
    print("Expanding "+str(len(idats_fnlist_filt))+" compressed IDATs...")
    for compidat in idats_fnlist_filt:
        print("Working on file "+compidat+", number "+str(nfile))
        idat_fn = os.path.splitext(compidat)[0];statuslist=[]; ridatd[compidat] = []
        with gzip.open(os.path.join(idatspath, compidat), 'rb') as f_in:
            with open(os.path.join(idatspath, idat_fn), 'wb') as f_out:
                try:
                    shutil.copyfileobj(f_in, f_out); ridatd[compidat].append(1)
                except:
                    ridatd[compidat].append(shutil.Error)
        print("Finished with file "+compidat+", number "+str(nfile));nfile+=1
    return ridatd

def new_idat_hlinks(gsmid, ts, igrn_fn, ired_fn):
    """ new_idat_hlinks

        Make new hlink files and return new hlink paths.

        Arguments
        * gsmid: Valid sample GSM ID (str).
        * ts: Timestamp (int or str).
        * igrn_fn: Filename of green intensity data file (IDAT, str).
        * ired_fn: Filename of red intensity data file (IDAT, str).
        
        Returns
        * rlist (list): List of path to new grn [0] and red [1] hlinked idats
    """
    lidat=[i for i in os.listdir(settings.idatspath) if "idat" in i and not "idat.gz" in i]
    lidatgsm=[i.split(".")[0] for i in lidat];setidatgsm=[gsm for gsm in set(lidatgsm)]
    nhlink=1; timestamp = gettime_ntp()
    for gsm in setidatgsm:
        lgsm = [i for i in lidat if gsm in i and not "hlink" in i]
        lgsmhlink = [i for i in lgsm if "hlink" in lgsm]
        if len(lgsmhlink) == 0:
            if len([i for i in lgsm if "Grn.idat" in i]) >= 1 & len([i for i in lgsm if "Red.idat" in i]) >= 1:
                print("Found paired IDATS...")
                igrn_fn = [i for i in lgsm if "Grn.idat" in i][0]
                ired_fn = [i for i in lgsm if "Red.idat" in i][0]
                print("Generating new hlink filenames...")
                gnewhlinkfn = '.'.join([gsm, str(ts), 'hlink', 
                    '.'.join(igrn_fn.split('.')[2:])])
                rnewhlinkfn = '.'.join([gsm, str(ts), 'hlink', 
                '.'.join(ired_fn.split('.')[2:])])
                grn_hlink_path = os.path.join(settings.idatspath, gnewhlinkfn)
                red_hlink_path = os.path.join(settings.idatspath, rnewhlinkfn)
                if not os.path.exists(grn_hlink_path) and not os.path.exists(red_hlink_path):
                    grn_hlink = os.link(os.path.join(settings.idatspath,igrn_fn), grn_hlink_path)
                    red_hlink = os.link(os.path.join(settings.idatspath,ired_fn), red_hlink_path)
                    print("Made new hlinks " + gnewhlinkfn + " and " + rnewhlinkfn)
                else:
                    print("Hlinks already exist for sample " + gsm)
            else:
                print("Couldn't find paired IDATs for gsm "+gsm)
        else:
            print("Found existing hlink files for gsm "+gsm)
        print("Finished with sample "+gsm+" number "+str(nhlink))
        nhlink+=1
    return None

if __name__ == "__main__":
    """ Process downloaded IDATs

    Process downloaded IDATs for a recountmethylation instance. Prepares IDATs
    for compilation.

    """
    expand_idats()
    new_idat_hlinks()
    expand_idats()