#!/usr/bin/env python3

""" edirect_query.py
    
    Authors: Sean Maden, Abhi Nellore
    
    Get files of GSE and GSM ids from GEO via edirect query (NCBI entrez 
    utlities function). IDs correspond to valid HM450k array experiments and 
    samples.
    
    Notes:
        * Scheduling: New edirect queries should be scheduled periodically to 
            check for latest experiment/sample/file info and to detect novel
            uploaded experiments/samples/files.
        * Filters: The edirect queries (GSE, GSM, and filtered query file) work
            together to form a filter on valid sample and experiment ids whose 
            files are to be downloaded and preprocessed. 
    
    Functions:
        * gse_query_diffs: Quickly detect and return differences between two 
            edirect query files.
        * gsm_query: Get valid HM450k sample (GSM) IDs, based on presence of 
            HM450k platform and availability of raw idat array files in 
            supplement.
        * gse_query: Get valid HM450k experiment (GSE) IDs, based on 
            specification of the correct platform accession (GPL13534) in the 
            experiment annotation. Note, many experiments combine multiple 
            platform ids, so it is necessary to filter out non-HM450k array 
            samples from experiment GSM ID lists.
        * gsequery_filter: Generate a new edirect query filter file, containing
            valid GSE and GSM IDs for HM450k array experiments/samples with
            idats available in sample supplemental files.
"""

import subprocess, os, socket, struct, sys, time, tempfile, atexit, shutil
import glob, filecmp; from itertools import chain
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, querydict, getlatest_filepath
import settings
settings.init()

def gse_query_diffs(query1, query2, rstat=False):
    """ gse_query_diffs
        
        Compares two GSE query results, returning query file diffs or boolean.
        
        Arguments:
            * query1 (str) : first edirect query, filename
            * query2 (str) : second edirect query, filename
            * rstat (True/False, bool.) : whether to return boolean only, 
                or else return list (default)
        
        Returns:
            * boolean (T/F) or query diffs (list of GSE IDs). Boolean is 'True' 
                if query objects are the same, 'False' otherwise
    """
    difflist = []
    qd1 = querydict(query1) # eg. for first value: qd1[list(qd1.keys())[0]]
    qd2 = querydict(query2)
    # if gse doesn't exist in qd2
    for key in qd1:
        evalkey = ''
        if key in list(qd2.keys()):
            if not qd1[key]==qd2[key]:
                for item in qd1[key]:
                    if not item in qd2[key]:
                        evalkey = False
                for item in qd2[key]:
                    if not item in qd1[key]:
                        evalkey = False
        else:
            evalkey = None
        if not evalkey:
            difflist.append(key)
    for key in qd2:
        if not key in list(qd1.keys()):
            difflist.append(key)
    if rstat:
        if len(difflist)>0:
            return False
        else:
            return True
    else:
        return difflist

def gsm_query(validate=True, timestamp=gettime_ntp()):
    """ gsm_query
        Get GSM level query object, from edirect query.
        Arguments:
            * validate (True/False, bool.) : whether to validate the file after 
                ownload.
            * timestamp (str) : NTP timestamp or function to retrieve it.
        Returns: 
            * Error (str) or download object (dictionary). 
    """
    # timestamp = str(gettime_ntp())
    eqdestpath = settings.equerypath
    temppath = settings.temppath
    os.makedirs(eqdestpath, exist_ok=True)
    os.makedirs(temppath, exist_ok=True)
    temp_make = tempfile.mkdtemp(dir=temppath)
    atexit.register(shutil.rmtree, temp_make)
    dldict = {}
    dldict['gsmquery'] = []
    dlfilename = ".".join(['gsm_edirectquery',timestamp])
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
    output=subprocess.check_output(args, shell=True)
    dldict['gsmquery'].append(output) 
    if validate:
        gsmquery_filewritten = os.path.join(temp_make,dlfilename)
        gsmquery_old = glob.glob('.'.join([os.path.join(eqdestpath, 'gsm_edirectquery'), '*',]))
        if gsmquery_old:
            if len(gsmquery_old)>1:
                gsmquery_old.sort(key=lambda x: int(x.split('.')[1]))
                gsmquery_old_mostrecent = gsmquery_old[-1]
            else:
                gsmquery_old_mostrecent = gsmquery_old[0]
            # filecmp should work (equesry file order preserved on reps)
            if filecmp.cmp(gsmquery_old_mostrecent,gsmquery_filewritten):
                print("Downloaded gsm query file same as most recent stored."+
                    " Removing..."
                    )
                os.remove(gsmquery_filewritten)
                dldict['gsmquery'].append(False)
            else:
                print("Downloaded file is new, moving to dest...")
                shutil.move(gsmquery_filewritten, os.path.join(
                                eqdestpath, os.path.basename(gsmquery_filewritten))
                            )
                dldict['gsmquery'].append(True)
        else:
            print("Downloaded file is new, moving...")
            shutil.move(gsmquery_filewritten, os.path.join(
                eqdestpath, os.path.basename(gsmquery_filewritten))
                )
            dldict['gsmquery'].append(True)
    return dldict

def gse_query(validate=True, timestamp=gettime_ntp()):
    """ gse_query
        
        Get GSE level query object from edirect query.
        
        Arguments:
            * validate (True/False, bool) : Whether to validate the file after 
                download.
            * timestamp (str) : NTP timestamp or function to retrieve it.
        
        Returns: 
            * Error (str) or download object (dictionary).
    """
    eqdestpath = settings.equerypath
    os.makedirs(eqdestpath, exist_ok=True)
    temppath = settings.temppath
    os.makedirs(temppath, exist_ok=True)
    temp_make = tempfile.mkdtemp(dir=temppath)
    atexit.register(shutil.rmtree, temp_make)
    dldict = {}
    dldict['gsequery'] = []
    dlfilename = ".".join(['gse_edirectquery',timestamp])
    dldict['gsequery'].append(dlfilename)
    subp_strlist1 = ["esearch","-db","gds","-query",
    "'"+settings.platformid+"[ACCN] AND idat[suppFile] AND gse[ETYP]'"
    ]
    subp_strlist2 = ["efetch","-format","docsum"]
    subp_strlist3 = ["xtract","-pattern","DocumentSummary",
        "-element","Id Accession",">",
        os.path.join(temp_make,dlfilename)
        ]
    args = " | ".join([" ".join(subp_strlist1),
        " ".join(subp_strlist2),
        " ".join(subp_strlist3)])
    output=subprocess.check_output(args, shell=True)
    dldict['gsequery'].append(output)
    if validate:
        gsequery_filewritten = os.path.join(temp_make,dlfilename)
        gsequery_old = glob.glob('.'.join([os.path.join(eqdestpath, 'gse_edirectquery'), '*',]))
        if gsequery_old:
            if len(gsequery_old)>1:
                gsequery_old.sort(key=lambda x: int(x.split('.')[1]))
                gsequery_old_mostrecent = gsequery_old[-1]
            else:
                gsequery_old_mostrecent = gsequery_old[0]
            # get diffs manually (edirect can return id's in different order)
            diffs = gse_query_diffs(query1=gsequery_old_mostrecent,
                query2=gsequery_filewritten,
                rstat=True)
            if diffs:
                print("Downloaded gse query file same as most recent stored."+
                    " Removing..."
                    )
                os.remove(gsequery_filewritten)
                dldict['gsequery'].append(False)
            else:
                print("Downloaded file is new, moving to dest...")
                shutil.move(gsequery_filewritten, os.path.join(
                                eqdestpath, os.path.basename(gsequery_filewritten))
                            )
                dldict['gsequery'].append(True)
        else:
            print("Downloaded file is new, moving...")
            shutil.move(gsequery_filewritten, os.path.join(
                eqdestpath, os.path.basename(gsequery_filewritten))
                )
            dldict['gsequery'].append(True)
    return dldict

def gsequery_filter(splitdelim='\t', timestamp=gettime_ntp()):
    """ gsequery_filter
        
        Prepare an edirect query file. Filter a GSE query file on its GSM 
            membership. 
        
        Arguments:
            * splitdelim (str) : Delimiter to split ids in querydict() call.
            * timestamp (str) : NTP timestamp or function to retrieve it.
        
        Returns:
            * gsequeryfiltered (list): Filtered GSE query object (list), writes
                filtered query file as side effect.
    """
    eqpath = settings.equerypath
    gsequerystr = settings.gsequerystr
    gsmquerystr = settings.gsmquerystr
    # get GSM list from gsm query file
    gsmqueryf_latestpath = getlatest_filepath(filepath=eqpath,
            filestr=gsmquerystr, embeddedpattern=True, tslocindex=1, 
            returntype='returnlist'
        )
    if gsmqueryf_latestpath:
        print("Latest gsmquery file detected: "+str(gsmqueryf_latestpath))
    else:
        print("Error detecting latest gsmquery file! Returning...")
        return
    gsmlines = [line.rstrip('\n') for line in open(gsmqueryf_latestpath[0])]
    gsmlist = [line.split('\t')[1::][0] for line in gsmlines] 
    # get GSE dictionary object
    gsequeryf_latestpath = getlatest_filepath(filepath=eqpath, filestr=gsequerystr,
            embeddedpattern=True, tslocindex=1, returntype='returnlist'
        )
    if gsequeryf_latestpath:
        print("Latest gsequery file detected: "+str(gsequeryf_latestpath))
    else:
        print("Error detecting latest gsequery file! Returning...")
        return
    gsed_obj = querydict(querypath=gsequeryf_latestpath[0], splitdelim='\t')
    gsefiltl = []
    for gsekey in list(gsed_obj.keys()):
        samplelist_original = gsed_obj[gsekey]
        samplelist_filt = [sample for sample in samplelist_original
            if sample in gsmlist
        ]
        if samplelist_filt and len(samplelist_filt)>0:
            gsefiltl.append(' '.join([gsekey,' '.join(samplelist_filt)]))      
    print('writing filt file...')
    if eqpath:
        filtfn = ".".join(["gsequery_filt",timestamp])
        with open(os.path.join(eqpath, filtfn), 'w') as filtfile:
            for item in gsefiltl:
                filtfile.write("%s\n" % item)
    return gsefiltl

if __name__ == "__main__":
    """ Run a new EDirect query

    Query the GEO DataSets API for valid data run using the target platform.

    """
    print("Beginning EDirect query...")
    equery_dest = settings.equerypath; temppath = settings.temppath
    gse_query(); gsm_query(); gsequery_filter()