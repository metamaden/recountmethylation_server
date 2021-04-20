#!/usr/bin/env python3

""" preprocess_mdat2.py
    Scripts for preprocesssing idats. 
    Sequence of operations (and relevant function or functions):
    1. Compile a status dictionary (see: 'scan_gsmstatdict')
    2. Run preprocessing using wrapper function (see: 'preprocess_mdat_wrapper')
    3. Append new data to compilations (see: 'append_compilations')
"""

import sys
import os
import inspect
import re
import shutil
import subprocess
import pickle
from datetime import datetime
import time
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
sys.path.insert(0, os.path.join("recount-methylation-analysis","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
from utilities import monitor_processes
import settings
settings.init()

def checkmdatpaths():
    """ checkdirs
        Verify mdat directories exist.
    """
    os.makedirs(settings.mdatdirpath, exist_ok=True)
    os.makedirs(settings.rawbetapath, exist_ok=True)
    os.makedirs(settings.rawgrnpath, exist_ok=True)
    os.makedirs(settings.rawredpath, exist_ok=True)
    os.makedirs(settings.rawpvalspath, exist_ok=True)
    os.makedirs(settings.noobbetapath, exist_ok=True)
    os.makedirs(settings.mdatlogspath, exist_ok=True)
    # adds M and UM tables
    os.makedirs(settings.rawmpath, exist_ok=True)
    os.makedirs(settings.rawumpath, exist_ok=True)
    os.makedirs(settings.normmpath, exist_ok=True)
    os.makedirs(settings.normumpath, exist_ok=True)

def getbn(maxbn=1000):
    """ getbn
        Returns list of valid hlink basenames retrieved from idats dir.
        Arguments:
        * maxbn (int): Upper limit to basenames returned.
        Returns:
        * bnlist (list): List of valid basenames.
    """
    # unit testing with local files
    idatspath = settings.idatspath
    # filter for valid hlink idat filenames
    flistfilt = os.listdir(idatspath)
    flistfilt = [fn for fn in flistfilt
        if 'hlink' in fn
        and fn[-4:]=='idat'
    ]
    bnlist = []
    gsmlist = list(set([fn.split('.')[0] for fn in flistfilt])) # unique gsms
    for index, gsm in enumerate(gsmlist):
        #print("starting GSM : "+gsm)
        gsmfnlist = [fn for fn in flistfilt
            if fn.split('.')[0]==gsm
        ]
        gsmgrn = [fn for fn in gsmfnlist
                if re.search('.*_Grn\\.idat$',fn)
            ]
        gsmred = [fn for fn in gsmfnlist
            if re.search('.*_Red\\.idat$',fn)
        ]
        # filt on valid gsms having grn and red idats
        if gsmgrn and gsmred:
            #print("found red and grn match. continuing...")
            gsmgrn = gsmgrn[0]
            gsmred = gsmred[0]
            if gsmgrn[:-9]==gsmred[:-9]:
                #print("Found valid basename for GSM id. continuing...")
                bnlist.append(gsmgrn[:-9])
            #print("appending new red and grn matches to bnlist, new len = "
            #    +str(len(bnlist)))
        print("Done with " + str(round(index/len(gsmlist),5))
                + " percent of GSMs; with bnlist length = "
                + str(len(bnlist)), end = "\r"
            )
        if len(bnlist) >= maxbn:
            print("Compiled basenames to max number. Returning...")
            return bnlist
    return bnlist

def preprocess_mdat(bnlistpass, timelim=40, nsampproc=10, nprocmax=4, statint=2):
    """ preprocess_mdat
        Preprocess mdat files via background subprocesses, monitoring, and 
        logging.
        Arguments
            * bnlistpass (list) : List of valid basenames
            * timelim (int) : Time limit for running processes, in minutes.
            * nsampproc (int): Number of samples per process launched.
            * nprocmax (int): Total processes to launch
            * statint (int): Seconds to wait before monitor status updates.
        Returns
            * None, produces status log in stdout and new logfile as side effect
    """
    # form the array of bn lists
    print("Forming basenames array...")
    bnscreenarray = [] # array of bn lists for batch processing
    n = nsampproc
    bnscreenarray = [' '.join(bnlistpass[i * n:(i + 1) * n]) for i in 
        range((len(bnlistpass) + n - 1) // n)
    ]
    bnscreenarray = bnscreenarray[0:nprocmax]
    print("Finished forming basenames array of length = "+str(len(bnscreenarray)))
    # new screen deployment
    print("Getting timestamp...")
    timestamp = gettime_ntp() 
    process_list = [] # process list for status monitoring and stderr
    # getting string limit
    print("Getting args maxstr...")
    argmaxstr = int(str(subprocess.check_output(['getconf', 
            'ARG_MAX'
        ])).replace("b'",'').replace("\\n'",""))
    print("Detected argmaxstr of "+str(argmaxstr)+". Continuing...")
    print("Launching background subprocesses...")
    for bi, bnstr in enumerate(bnscreenarray, 0):
        cmd = ['Rscript', settings.mdatscriptpath, timestamp, str(bi), 
            '"'+bnstr+'"'
        ]
        cmdcharlen = len(''.join(cmd))
        print("Formed cmd str of len = " + str(cmdcharlen)
            +", checking args str limit...")
        if cmdcharlen <= argmaxstr:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE
                )
            process_list.append(proc)
            print("Launched background subprocess and appended poll to "
                +"statuslist. Continuing...")
        else:
            print("Error! Char length of cmd exceeds system limit for args. "
                +"Try modifying argument 'nsampproc'. Continuing...")
    # process monitoring start
    monitor_processes(process_list=process_list, logpath=settings.mdatlogspath)
    print("Completed preprocessing. Returning...")
    return None

def preprocess_mdat_wrapper(filttype='notrun', timelim=1000, nprocmax=4, 
    nsampproc=10, statint=2):
    """ preprocess_mdat_wrapper
        Run preprocessing batches as child processes, processively. Each batch 
        parallelizes among an array sublist passed to it by the wrapper.
        This function automatically detects valid basenames, forms basenames 
        into sublists, and passes these sublists to preprocess_mdat for 
        parallelization within the sublist.
        Arguments:
        * filttype: How to pre-filter basename files to be run. Options:
            - 'notrun' : any samples lacking all mdata
            - [dattype] : any samples lacking a specific datatype, one of either
                'detp', 'rawgrn', 'rawred', 'rawbeta', 'mraw', 'umraw', 
                'noobbeta', 'mnorm', or 'umnorm'.
        * timelim (int): Max time (minutes) to monitor processes
        * nprocmax (int): Max processes for child process to deploy
        * nsampproc (int): Max samples per child process
        Returns: 
        * None. Creates logs for preprocessing batches as side-effect.
    """
    # verify paths
    print("Verifying mdat paths...")
    checkmdatpaths()
    # prefilter 
    print("Loading GSM status pickle dictionary...")
    gsmstatdict = pickle.load( open( settings.gsmstatpicklepath, "rb" ) )
    print("Forming basenames from status object...")
    longbnlist = [k.replace('"','') for k in list(gsmstatdict.keys())]
    gsmfiltlist = []
    lbnfilt = []
    print("Applying filter to basenames...")
    if filttype == 'notrun':
        gsmfiltlist = [k.replace('"','') for k in list(gsmstatdict.keys())
            if gsmstatdict[k] == []
        ]
        lbnfilt = [bn for bn in longbnlist if bn in gsmfiltlist]
    else:
        gsmfiltlist = [k.replace('"','') for k in list(gsmstatdict.keys())
            if not filttype in gsmstatdict[k]
        ]
        lbnfilt = [bn for bn in longbnlist if bn in gsmfiltlist]
    if lbnfilt:
        print("Generated filtered bn list of length "+str(len(lbnfilt)))
    else:
        print("There was a problem forming the filtered basenames list. "
                +"Please check entry for argument 'filttype' is valid."
            )
        return None
    # Form array of sub lists, to process in batches
    print("Forming array of sublists for preprocessing...")
    longbnarray = []
    totalbn = len(lbnfilt)
    bnpass = nprocmax*nsampproc # bn per process passed
    # form array of bn sublists corresponding to each process launched
    n = bnpass
    longbnarray = [' '.join(lbnfilt[i * n:(i + 1) * n]) for i in 
        range((totalbn + n - 1) // n)
    ]
    print("Finished forming large bn array of length = " 
            + str(len(longbnarray)) + " to pass for preprocessing. Continuing.."
        )
    print("Starting batch processing of basenames sublists...")
    # launch processes successively, with internal parallelization by 'nprocmax' 
    for si, sublist in enumerate(longbnarray, 0):
        print("Working on sublist "+str(si))
        preprocess_mdat(bnlistpass=sublist.split(' '), timelim=timelim, 
                nsampproc=nsampproc, nprocmax=nprocmax, statint=statint
            )
        print("Finished processing sublist "+str(si))
    print("Finished processing all basenames sublists. Returning...")
    return None

def check_and_form_compilations():
    """ check_and_form_compilations
        Detect available compilations, and form new compilations from subset
        tables if valid compilation not available for filetype.
        Arguments:
        * None
        Returns:
        * None, checks and/or forms compilation files as side effect.
    """
    # check compilations dir
    os.makedirs(settings.compilationspath, exist_ok=True)
    # get valid compilations filenames
    compflist = os.listdir(settings.compilationspath)
    ctypelist = ["detp", "rawgrn", "rawred", "rawbeta", "noobbeta", "mraw", 
        "umraw", "mnorm", "umnorm"
    ]
    ctypepaths = [settings.rawpvalspath, settings.rawgrnpath, 
        settings.rawredpath, settings.rawbetapath, settings.noobbetapath,
        settings.rawmpath, settings.rawumpath, settings.normmpath, 
        settings.normumpath
    ]
    print("Loading gsm status dictionary...")
    gsmstatdict = pickle.load(open( settings.gsmstatpicklepath, "rb"))
    print("Beginning iterations on compilation filetypes...")
    for index, ctype in enumerate(ctypelist, 0):
        print("Starting check for compilation of type " + ctype + "...")
        fpattern = '.*\\.'+ctype+'\\.compilation\\.mdat$'
        matchlist = list(filter(re.compile(fpattern).match, compflist))
        if not matchlist:
            print("No compilation found of type " + ctype + ". Continuing...")
            pflist = os.listdir(ctypepaths[index])
            pfn = list(filter(re.compile('.*\\.mdat$').match, pflist))[0]
            if pfn:
                print("Found valid preprocessed file of type " + ctype 
                        + ". Continuing..."
                    )
                ts = pfn.split('.')[0]
                cfn = '.'.join([str(ts), ctype, 'compilation', 'mdat'])
                shutil.copyfile(os.path.join(ctypepaths[index], pfn),
                        os.path.join(settings.compilationspath, cfn)
                    )
                if os.path.exists(os.path.join(settings.compilationspath, cfn)):
                    print("Successfully created new compilations file of type "
                            + ctype + ". Updating gsm status dictionary..."
                        )
                    with open(os.path.join(settings.compilationspath, cfn), "r") as oc:
                        for li, line in enumerate(oc, 0):
                            if line[0:4]=='"GSM':
                                itemkey = line.split(' ')[0].replace('"','')   
                                if itemkey in list(gsmstatdict.keys()):
                                    if not ctype in gsmstatdict[itemkey]:
                                        gsmstatdict[itemkey].append(ctype)
                                    else:
                                        gsmstatdict[itemkey] = [ctype]
                                print("Finished with cfile line "+str(li), 
                                    end="\r")
                    continue
                else:
                    print("There was an error copying the preprocessed file for type "
                            + ctype + ". Continuing..."
                        )
                    continue
            else:
                print("Could not find valid file of type " + ctype
                        + "Please run preprocessing prior to assembling "
                        + "compilations."
                    )
                continue
    print("Finished checking compilations. Returning...")
    return None

def scan_gsmstatdict(usersheet=True, maxbn=40000,
        gsmstatdictpath=settings.gsmstatpicklepath):
    """ scan_gsmstatdict
        Make a new GSM status dictionary, or update an existing dictionary with
        latest sample data from compilations files.
        Arguments:
        * usersheet (Bool.): Whether to load sample basenames from latest 
            detected rsheet. If 'False', detect basenames de novo with 
            "getbn()".
        * maxbn (int): Max basenames allowed when forming new status dictionary.
        * gsmstatdictpath (path/str): Path from which to read status dictionary.
        Returns: 
        * None, or status dictionary object if loadobj
    """
    if not os.path.exists(gsmstatdictpath):
        basenames = []
        if usersheet:
            rslatest = getlatest_filepath(filepath=settings.sheetspath,
                filestr="rsheet", embeddedpattern=True, tslocindex=0,
                returntype='returnlist'
            )
            if rslatest:
                rslpath = rslatest[0]
                print("Detected latest rsheet. Reading sample ids...")
                with open(rslpath, "r") as rso:
                    for linect, line in enumerate(rso, 1):
                        if line[0:3]=='GSM':
                            basenames.append(line.split(' ')[7].replace('\n',''))
                        print("Finished reading line num "+str(linect), 
                            end="\r")
                print("Finished reading rsheet. Continuing...")
        else:
            # form the new status dictionary
            print("Getting basenames with 'getbn()'...")
            basenames = getbn(maxbn=maxbn)
        if not basenames:
            print("Error obtaining basenames. Returning...")
            return None
        else:
            print("Finished retrieving n = "+str(len(basenames))
                +" basenames. Forming dictionary...")
            gsmstatdict = {bn:[] for bn in basenames}
            pickle_out = open(gsmstatdictpath, "wb")
            pickle.dump(gsmstatdict, pickle_out)
            pickle_out.close()
    # check path for existing file
    if os.path.exists(gsmstatdictpath):
        print("Detected sample status dictionary. Updating...")
        tasktype = "update dictionary"
        gsmstatdict = pickle.load(open(gsmstatdictpath, "rb"))
        cflist = os.listdir(settings.compilationspath)
        cflist = [cfn for cfn in cflist if 'compilation' in cfn]
        for cfn in cflist:
            print("Starting on cfn "+str(cfn))
            cftype = cfn.split('.')[1]
            print("Detected compilation type "+str(cftype))
            cfnpath = os.path.join(settings.compilationspath, cfn)
            with open(cfnpath, "r") as opencfn:
                for li, line in enumerate(opencfn, 1):
                    if line.split(' ')[0][0:4]=='"GSM':
                        gsmfname = line.split(' ')[0].replace('"','')
                        if gsmfname in list(gsmstatdict.keys()):
                            if not cftype in gsmstatdict[gsmfname]:
                                gsmstatdict[gsmfname].append(cftype)
                        else:
                            gsmstatdict[gsmfname] = [cftype]
                        print("Finished reading line num "+str(li), end="\r")
            print("Finished reading lines from cfn. "
                    +"Saving updated sample status dictionary.")
            pickle_out = open(gsmstatdictpath, "wb")
            pickle.dump(gsmstatdict, pickle_out)
            pickle_out.close()
            print("Finished saving updated dictionary. Continuing...")
    else:
        print("Error, could not detect gsm status dictionary at settings path. "
            +"Returning...")
    return None

def append_compilations():
    """ append_compilations
        Append new data to compilations from subset tables, updating the sample 
        status dictionary in the process.
        Arguments:
        * None
        Returns:
        * None, appends new data to compilations as side effect.
    """
    gsmstatdict = pickle.load(open(settings.gsmstatpicklepath, "rb"))
    if gsmstatdict:
        print("Successfully loaded GSM status pickle dictionary. Continuing...")
    else:
        print("Error loading GSM status pickle dictionary. Returning...")
        return None
    stdirlist = [settings.rawbetapath, settings.rawpvalspath, settings.rawredpath,
        settings.rawgrnpath, settings.noobbetapath, settings.rawmpath,
        settings.rawumpath, settings.normmpath, settings.normumpath
    ]
    # check and form new compilations
    check_and_form_compilations()
    cfiles = os.listdir(settings.compilationspath) # compilation files
    cfiles = [cfn for cfn in cfiles if 'compilation' in cfn]
    print("Found n = "+str(len(cfiles))+" compilation filenames. Continuing...")
    # iterate over subset table directory paths
    for stdirpath in stdirlist:
        print("Working on subset table dirpath: "+str(stdirpath))
        stlist = os.listdir(stdirpath)
        stlist = [stfn for stfn in stlist if 'mdat' in stfn]
        ctype = os.path.basename(stdirpath)
        cfname = [cf for cf in cfiles if ctype in cf][0]
        if cfname:
            print("Detected compilation file of type: "+str(cfname))
        else:
            print("Error extracting cfname for type "+str(ctype)+", breaking..")
            break
        print("reading in compilation file colnames...")
        with open(os.path.join(settings.compilationspath, cfname), "r") as rcf:
            cnamecomp = rcf.readline() # compilation file colnames
        for st in stlist:
            print("Working on subset table file "+str(st))
            with open(os.path.join(stdirpath, st), "r") as rst:
                stcnames = rst.readline() # subset table colnames
            # check colnames
            if stcnames == cnamecomp:
                print("Subset and comp table colnames identical. Continuing...")
                flist = []
                # open subset table for reading new data
                with open(os.path.join(stdirpath, st), "r") as rst:
                    linect = 0
                    for line in rst:
                        flist.append(line)
                        print("Finished appending line "+str(linect), 
                                end="\r")
                        linect += 1
                print("Finished appending n = "+str(len(flist))+" lines. "
                        +"Checking and appending new data..")
                # open compilation file for appending new data
                with open(os.path.join(settings.compilationspath, cfname), "a") as rcf:
                    for itemi, item in enumerate(flist, 0):
                        if item.split(' ')[0][0:4]=='"GSM':
                            # check if sample data is new, and update statdict
                            itemkey = item.split(' ')[0].replace('"','')   
                            if itemkey in list(gsmstatdict.keys()):
                                if not ctype in gsmstatdict[itemkey]:
                                    rcf.write(item)
                                    gsmstatdict[itemkey].append(ctype)
                                    continue
                            else:
                                rcf.write(item)
                                gsmstatdict[itemkey] = [ctype]
                                continue
                        print("Finished with item "+str(itemi), end = "\r")
            else:
                print("Error, colnames don't match between compilation and st."
                        +" Breaking...")
                break
            print("Finished with st "+str(st)+". Updating status dictionary...")
            # save the status dictionary between subset table iterations
            pickle_out = open(settings.gsmstatpicklepath, "wb")
            pickle.dump(gsmstatdict, pickle_out)
            pickle_out.close()
            print("Finished updating sample status dictionary. Continuing...")
            continue      
        print("Finished processing all subset tables at path. Continuing...")
    print("Finished processing all subset table paths. Returning...")    
    return None

def purgemdatdir(type="subsettables"):
    """ purgemdatdir
        Purge preprocessing files from mdata dir in perparation for new 
        analysis.
        Argumernts:
        * type (str): What filetype(s) to purge (options: 'subsettables', 'all', 
            'compilations')
        Returns:
        * None, purges files as side effect
    """
    ml = os.listdir(settings.mdatdirpath)
    validdirlist = ['compilations', 'rawbeta', 'detp', 'mraw', 'umraw', 
        'rawgrn', 'rawred', 'noobbeta', 'mnorm', 'umnorm', 'logs'
    ]
    mlfilt = [dn for dn in ml if dn in validdirlist]
    if type=="subsettables":
        mlfilt = [dn for dn in mlfilt if not "compilations" in dn]
    if type=="compilations":
        mlfilt=["compilations"]
    print("After type detection, working on list of "+str(len(mlfilt))
        +" files.")
    for dn in mlfilt:
        print("Removing contents of dir: "+str(dn))
        cmd = ' '.join(['rm', '-rf', os.path.join(settings.mdatdirpath, dn, 
            "*")])
        print(cmd)
        subprocess.call(cmd, shell=True)
    # rm gsm status dict
    subprocess.Popen(' '.join(['rm', settings.gsmstatpicklepath]), shell=True)
    print("Finished purging mdata dir. Returning...")
    return None






