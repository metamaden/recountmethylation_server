import os
import re
import gzip
import shutil
import subprocess

def gettime_ntp(addr='time.nist.gov'):
    """ Get NTP Timestamp,
        code from:
        <https://stackoverflow.com/questions/39466780/simple-sntp-python-script>
        
        Arguments
            * addr : valid NTP address ('0.uk.pool.ntp.org','time.nist.gov' etc)
        Returns
            * timestamp : NTP seconds timestamp of type 'int'
    """
    TIME1970 = 2208988800
    client = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    data = '\x1b' + 47 * '\0'
    client.sendto(data.encode('utf-8'), (addr, 123))
    data, address = client.recvfrom(1024)
    t = struct.unpack('!12I', data)[10] - TIME1970
    return t

def expand_soft(soft_dir='soft',rmcompressed=True):
    """ Automatically expand detected soft files
        Arguments
            * soft_dir : directory to search for soft files
            * rmcompressed : remove compressed files successfully expanded
        Returns
            * list : compressed filenames and statuslist
    """
    r1 = re.compile(".*soft.*")
    dirlist = os.listdir(soft_dir)
    soft_list = list(filter(r1.match, dirlist)) 
    r2 = re.compile(".*\.gz$")
    softl_compressed = list(filter(r2.match,soft_list))
    if len(softl_compressed)>0:
        statuslist = []
        for softcompfile in softl_compressed:
            with gzip.open(os.path.join(soft_dir,softcompfile), 'rb') as f_in:
                with open(os.path.join(soft_dir,softcompfile[:-3]), 'wb') as f_out:
                    try:
                        shutil.copyfileobj(f_in, f_out)
                        statuslist.append(True)
                    except shutil.Error as se:
                        statuslist.append(se)
        statusindices = [i for i, x in enumerate(statuslist) if x == True]
        rmsuccess = [softl_compressed[i] for i in statusindices]
        if rmcompressed & len(rmsuccess)>0:
            for compfilename in rmsuccess:
                os.remove(os.path.join(soft_dir,compfilename))
    else: 
        print("Error: no compressed soft files found at specified soft_dir.")
        return 
    return [softl_compressed,statuslist]

def extract_gsm_soft(gse_soft_dir='gse_soft',gsm_soft_dir='gsm_soft',
    temp_dir='temp_dir',edir_filt,validate=True):
    """ Extract GSM soft file sections from GSE soft file
        Arguments 
            * gse_soft_dir : directory to search for gse soft files
            * gsm_soft_dir : dir to write new gsm soft files
        Returns
            * null, generates gsm soft files as side effect
    """
    timestamp = str(gettime_ntp())
    r1 = re.compile(".*soft.*")
    dirlist = os.listdir(soft_dir)
    soft_list = list(filter(r1.match, dirlist))
    for soft in soft_list:
        openindex = []
        closeindex = []
        rxopen = re.compile('!Sample_title')
        rxclose = re.compile('!sample_table_end')
        lsoft = []
        with open(softfile) as file:
            for num, line in enumerate(file, 0):
                if rxclose.search(line):
                    closeindex.append(num)
                if rxopen.search(line):
                    openindex.append(num)
                lsoft.append(line)
        rxgsm = re.compile('GSM[0-9]*')
        for num, openi in enumerate(openindex,0):
            # read gsm lines
            gsmsoftlines = lsoft[openi:closeindex[num]]
            gsm_newfile_path = os.path.join(
                gsm_soft_dir,
                timestamp+str(rxgsm.findall(gsmsoft[1])[0])+'.soft'
                )
            with open(gsm_newfile_path,"w+") as gsmfile:
                for line in gsmsoftlines:
                    gsmfile.write(line)
                gsmfile.close()
    return

def gsm_soft2json(gsm_soft_fn,gsm_json_destdir='gsm_soft_json'):
    """ Convert GSM soft file to JSON format
        Call R script to process GSM soft (XML) to JSON format
        Arguments
            * gsm_soft_fn : name of GSM soft file to process
            * gsm_json_destdir : name of dest dir for new JSON files
        Returns
            * null, or error
    """
    try:
        cmdlist = ['Rscript','soft2json.R',gsm_soft_fn,gsm_json_destdir]
        subprocess.call(cmdlist,shell=False)
        #subprocess.check_call(cmdlist,executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        raise e

""" Examples and tests

import os
import re
import gzip
import shutil
import subprocess

from process_soft import expand_soft
from process_soft import process_soft

process_soft()


soft_dir = 'soft'
softfile = "GSE109904.1542514817.GSE109904_family.soft"
subprocess.check_call(['process_gsm_from_soft.sh',os.path.join(soft_dir,softfile)], executable='/bin/bash')


gsm_soft2json('./gsm_soft/GSM3177405.soft')

"""
