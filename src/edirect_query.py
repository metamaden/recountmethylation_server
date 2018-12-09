import subprocess
import os
import socket
import struct
import sys
import time
import tempfile
import atexit
import shutil
import glob
import filecmp
from itertools import chain
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from dl import gettime_ntp

def querydict(query, splitdelim = '\t'):
    """ Process edirect query result into dictionary of GSEs
        Arguments
            * query (file) : filepath of edirect results 
            (format: 1 GSE per newline)
        Returns
            * query (dictionary) : dictionary (keys = GSEs, vals = GSM lists)
    """
    querylines = [line.rstrip('\n') for line in open(query)]
    querylist = []
    querydict = {}
    for line in querylines:
        querylist.append(line.split(splitdelim))
    for idlist in querylist:
        gsekey = list(filter(lambda x:'GSE' in x, idlist))[0]
        gsmlist = list(filter(lambda x:'GSM' in x, idlist))
        querydict[gsekey] = gsmlist
    return querydict

def gse_query_diffs(query1,query2,rstat=False):
    """ Compares two GSE query results, returning diffs or boolean
        Arguments
            * query1 : first edirect query, filename
            * query2 : second edirect query, filename
            * rstat : whether to return boolean only, else return list
        Returns
            * boolean (T/F) or diffs object (list of GSEs with diffs)
            Boolean is 'True' if query objects are the same, 'False' otherwise
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
                        evalkey = 'notalike'
                for item in qd2[key]:
                    if not item in qd1[key]:
                        evalkey = 'notalike'
        else:
            evalkey = 'notpresent'
        if (evalkey == 'notalike') or (evalkey == 'notpresent'):
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

def gsm_query(dest = 'equery', temp = 'temp', validate = True, 
    timestamp = str(gettime_ntp())):
    """ Get GSM-level query object from Edirect
        Arguments
            * dest : destination directory for query object
            * temp : temporary dir becfore validation
            * validate : (T/F) whether to validate the file after download
        Retursn 
            * err or successful download 
    """
    # timestamp = str(gettime_ntp())
    os.makedirs(dest, exist_ok=True)
    os.makedirs(temp, exist_ok=True)
    temp_make = tempfile.mkdtemp(dir=temp)
    atexit.register(shutil.rmtree, temp_make)
    dldict = {}
    dldict['gsmquery'] = []
    dlfilename = ".".join(['gsm_edirectquery',timestamp])
    dldict['gsmquery'].append(dlfilename)
    subp_strlist1 = ["esearch","-db","gds","-query",
    "'GPL13534[ACCN] AND idat[suppFile] AND gsm[ETYP]'"
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
        gsmquery_old = glob.glob('.'.join([os.path.join(dest, 'gsm_edirectquery'), '*',]))
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
                                dest, os.path.basename(gsmquery_filewritten))
                            )
                dldict['gsmquery'].append(True)
        else:
            print("Downloaded file is new, moving...")
            shutil.move(gsmquery_filewritten, os.path.join(
                dest, os.path.basename(gsmquery_filewritten))
                )
            dldict['gsmquery'].append(True)
    return dldict

def gse_query(dest = 'equery', temp = 'temp', validate = True, 
    timestamp = str(gettime_ntp())):
    """ Get GSE-level query object from Edirect
        Arguments
            * dest : destination directory for query object
            * temp : temporary dir becfore validation
            * validate : (T/F) whether to validate the file after download
        Retursn 
            * err or successful download 
    """
    # timestamp = str(gettime_ntp())
    os.makedirs(dest, exist_ok=True)
    os.makedirs(temp, exist_ok=True)
    temp_make = tempfile.mkdtemp(dir=temp)
    atexit.register(shutil.rmtree, temp_make)
    dldict = {}
    dldict['gsequery'] = []
    dlfilename = ".".join(['gse_edirectquery',timestamp])
    dldict['gsequery'].append(dlfilename)
    
    subp_strlist1 = ["esearch","-db","gds","-query",
    "'GPL13534[ACCN] AND idat[suppFile] AND gse[ETYP]'"
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
        gsequery_old = glob.glob('.'.join([os.path.join(dest, 'gse_edirectquery'), '*',]))
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
                                dest, os.path.basename(gsequery_filewritten))
                            )
                dldict['gsequery'].append(True)
        else:
            print("Downloaded file is new, moving...")
            shutil.move(gsequery_filewritten, os.path.join(
                dest, os.path.basename(gsequery_filewritten))
                )
            dldict['gsequery'].append(True)
    return dldict

def gsequery_filter(gsequery='gsequery',gsmquery='gsmquery',target='equery',
    splitdelim = '\t', timestamp = str(gettime_ntp())):
    """ Prepare an edirect query file
        Filter GSE query on GSM query membership
        Arguments
            * gsequery (str): path to GSE query file from edirect.sh
            * gsmqeury (str): path to GSM query file from edirect.sh
            * splitdelim (str) : delimiter to split ids in querydict
        Returns
            * gsequeryfiltered (list): filtered GSE query file
    """
    # timestamp = str(gettime_ntp())
    gsed = querydict(query = gsequery,splitdelim = '\t')
    gsmlines = [line.rstrip('\n') for line in open(gsmquery)]
    gsmlist = [] # list of valid gsms with idat suppfiles
    for line in gsmlines:
        gsmlist.append(line.split('\t')[1::][0])
    gsefiltl = []
    for gsekey in list(gsed.keys()):
        ival = gsed[gsekey]
        # print(gsmkey)
        if len(ival)>0:
            for gsm in ival:
                if not gsm in gsmlist:
                    ival.remove(gsm)
            if len(gsed[gsekey])>0:
                gsefiltl.append(' '.join([gsekey,' '.join(ival)]))
    print('writing filt file...')
    if target:
        towrite = ".".join(["gsequery_filt",timestamp])
        with open(os.path.join(target,towrite), 'w') as filtfile:
            for item in gsefiltl:
                filtfile.write("%s\n" % item)
    # to read in the gsequery_filt file:
    # gsefiltd = querydict(query = 'gsequery_filt.t',splitdelim = ' ')
    return gsefiltl

"""

import subprocess
import os
import socket
import struct
import sys
import time
import tempfile
import atexit
import shutil
import glob
import filecmp

from edirect_query import gettime_ntp, gsm_query, gse_query, querydict, gse_query_diffs

gsm_query()

gse_query()

gse_query_diffs(query1='./equery/gse_edirectquery.timestamp2',
    query2='./temp/tmpavadzdzh/gse_edirectquery.sometime',
    rstat=True)


#from edirect_query import equery_gse
#form edirect_query import equery_gsm

dest_dir = 'equery'
temp_dir = 'temp'
dest = dest_dir
temp = temp_dir

os.makedirs(dest, exist_ok=True)
os.makedirs(temp, exist_ok=True)
temp_make = tempfile.mkdtemp(dir=temp)
atexit.register(shutil.rmtree, temp_make)

# dlfilename = ".".join(['gse_edirectquery',timestamp])
dlfilename = "gsefile"

#---------------
# pipes testing
#---------------
# refs
# 1. <https://security.openstack.org/guidelines/dg_avoid-shell-true.html>


subp_strlist1 = ["esearch","-db","gds","-query",
    "'GPL13534[ACCN] AND idat[suppFile] AND gse[ETYP]'"
    ]
subp_strlist2 = ["efetch","-format","docsum"]
subp_strlist3 = ["xtract","-pattern","DocumentSummary",
    "-element","Id Accession",">",
    os.path.join(temp_make,dlfilename)
    ]

pipe1 = subprocess.Popen(subp_strlist1, stdout=subprocess.PIPE, shell=False,
    stderr=subprocess.PIPE
    )
pipe2 = subprocess.Popen(subp_strlist2, stdin=pipe1.stdout, 
    stdout=subprocess.PIPE, shell=False, stderr=subprocess.PIPE
    )
pipe3 = subprocess.Popen(subp_strlist3, stdin=pipe2.stdout, 
    stdout=subprocess.PIPE, shell=False, stderr=subprocess.PIPE
    )
pipe1.stdout.close()
output = pipe3.communicate()[0]

# try 2
subp_strlist1 = ["esearch","-db","gds","-query",
    "'GPL13534[ACCN] AND idat[suppFile] AND gse[ETYP]'"
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


equery_gse(dest=dest_dir,temp=temp_dir,vlaidate=True)
equery_gse(dest=dest_dir,temp=temp_dir,vlaidate=True)

#---------
# validate test
#---------
dest = 'equery'
if validate:
    gsequery_filewritten = os.path.join(temp_make,dlfilename)
    gseqeury_old = glob.glob('.'.join([os.path.join(dest, 'gse_edirectquery'), '*',]))
    if gsequery_old:
        if len(gsequery_old)>1:
            gseqeury_old.sort(key=lambda x: int(x.split('.')[1]))
            gsequery_old_mostrecent = gseqeury_old[-1]
        else:
            gsequery_old_mostrecent = gseqeury_old
        if filecmp.cmp(gsequery_filewritten, gsequery_old_mostrecent):
            print("Downloaded gse query file same as most recent stored."+
                " Removing..."
                )
            os.remove(gsequery_filewritten)
            dldict[index][1] = False
        else:
            print("Downloaded file is new, moving to dest...")
            shutil.move(gsequery_filewritten, os.path.join(
                            dest, os.path.basename(gsequery_filewritten))
                        )
            dldict[index][1] = True
    else:
        print("Downloaded file is new, moving...")
        shutil.move(gsequery_filewritten, os.path.join(
            dest, os.path.basename(gsequery_filewritten))
            )
        dldict[index][1] = True




if validate:
        print("Validating downloaded files...")
        for gsm_id, file_written, index in files_written:
            print("file written is "+file_written)
            gsms = glob.glob('.'.join([
                        os.path.join(dest_dir, gsm_id), '*', '.'.join(
                                        file_written.split('.')[2:]
                                )
                    ]))
            gsmstr = "; ".join(str(e) for e in gsms)
            print("gsms found : "+gsmstr)
            if gsms:
                if len(gsms)>1:
                    gsms.sort(key=lambda x: int(x.split('.')[1]))
                    most_recent = gsms[-1]
                else:
                    most_recent = gsms[0]
                print('most recent file :'+most_recent)
                if filecmp.cmp(most_recent, file_written):
                    print("Downloaded file is same as recent file. Removing...")
                    os.remove(file_written)
                    # If filename is false, we found it was the same
                    dldict[gsm_id][index][1] = False
                else:
                    print("Downloaded file is new, moving to dest_dir...")
                    shutil.move(file_written, os.path.join(
                            dest_dir, os.path.basename(file_written))
                        )
                    dldict[gsm_id][index][1] = True
            else:
                print("Downloaded file is new, moving...")
                shutil.move(file_written, os.path.join(
                            dest_dir, os.path.basename(file_written))
                        )
                dldict[gsm_id][index][1] = True


# example from SO
from subprocess import Popen
from subprocess import PIPE
p1 = Popen(["dmesg"], stdout=PIPE)
p2 = Popen(["grep", "hda"], stdin=p1.stdout, stdout=PIPE)
p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
output = p2.communicate()[0]




















"""