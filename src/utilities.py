#!/usr/bin/env python3

""" utilities.py
    Commonly referenced functions for recount methylation server.
    Functions:
        * gettime_ntp: Return an NTP timestamp as a string, for file versioning.
        * getlatest_filepath: Access the path to the latest version of a file in
            a provided files directory. References NTP timestamp in filename to 
            determine latest available file.
        * querydict: Form a dictionary object by reading in an edirect file.
        * get_queryfilt_dict: Retrieve latest edirect query filtered file.
"""

import os
import glob
import socket
import struct
import sys
import time
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import settings
settings.init()

def gettime_ntp(addr='time.nist.gov'):
    """ gettime_ntp
        Get NTP Timestamp for file versioning.
        Arguments
            * addr (str) : valid NTP address (e.g. '0.uk.pool.ntp.org',
                'time.nist.gov' etc)
        Returns
            * timestamp (str) : NTP seconds timestamp, converted to string.
    """
    TIME1970 = 2208988800
    client = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    data = '\x1b' + 47 * '\0'
    client.sendto(data.encode('utf-8'), (addr, 123))
    data, address = client.recvfrom(1024)
    t = struct.unpack('!12I', data)[10] - TIME1970
    return str(t)

def getlatest_filepath(filepath, filestr, embeddedpattern=False, tslocindex=1,
    returntype='returnstr'):
    """ getlatest_filepath
        Get path the latest version of a file, based on its timestamp. Can 
        return >1 files sharing latest timestamp as type list or str.
        Arguments:
            * filepath (str) : Path to directory to search.
            * filestr (str) : Pattern of filename filter for search.
            * embeddedpattern (T/F, bool) : Whether filestr pattern is embedded
                in filename (assumes pattern at start of name otherwise)
            * tslocindex (int) : Relative location index of timestamp in fn.
            * returntype (str) : Return file path(s) as str (use 'returnstr') or 
                list (use 'returnlist')
        Returns:
            * latest_file_path (str) or status (0, int) : Path to latest version 
                of file, or else 0 if search turned up no files at location
    """
    # print("start getlatest")
    if embeddedpattern:
        embedpattern = str('*' + filestr + '*')
        pathstr = str(os.path.join(filepath, embedpattern))
        filelist = glob.glob(pathstr)
    else:
        pathpattern = '.'.join([os.path.join(filepath, filestr), '*'])
        filelist = glob.glob(pathpattern)
    if filelist:
        flfilt = []
        # filter filelist on possible int/valid timestamp
        for fp in filelist:
            try:
                int(os.path.basename(fp).split('.')[tslocindex])
                flfilt.append(fp)
            except ValueError:
                break
        if flfilt and not len(flfilt)==0: 
            if len(flfilt) > 1:
                lfr = []
                # sort on timestamp
                flfilt.sort(key=lambda x: int(os.path.basename(x).split('.')[tslocindex]))
                # last list item
                lastitem = flfilt[-1]
                latestts = lastitem.split('.')[tslocindex]
                lfr = [flitem for flitem in flfilt
                    if flitem.split('.')[tslocindex] == latestts
                ]
            else:
                lfr = [flfilt[0]]
            # print("getlatest end")
            if returntype=='returnstr':
                return ' '.join(i for i in lfr)
            if returntype=='returnlist':
                return lfr
            else:
                return None
        else:
            return None
    else:
        return None 

def querydict(querypath, splitdelim='\t'):
    """ querydict
        Read edirect query results file into a dictionary object.
        Arguments:
            * querypath (str) : Filepath of edirect results (format: 1 GSE per 
                newline).
        Returns:
            * query (dictionary) : Dictionary (format : keys = GSEs, vals = GSM 
                lists)
    """
    querylines = [line.rstrip('\n') for line in open(querypath)]
    querylist = []
    querydict = {}
    for line in querylines:
        querylist.append(line.split(splitdelim))
    for idlist in querylist:
        gsekey = list(filter(lambda x:'GSE' in x, idlist))[0]
        gsmlist = list(filter(lambda x:'GSM' in x, idlist))
        querydict[gsekey] = gsmlist
    return querydict

def get_queryfilt_dict(filesdir='recount-methylation-files',
    eqtarget='equery'):
    """ get_queryfilt_dict
        Return the latest filtered GSE query file as a dictionary.
        Arguments:
            * filesdir: Root name of directory containing database files.
            * eqtarget: Name of equery files destination directory.
        Returns:
            * gsefiltd (dict): GSE filtered query, as a dictionary object.
    """
    # eqpath = os.path.join('recount-methylation-files','equery')
    eqpath = settings.equerypath
    gsefilt_latest = getlatest_filepath(eqpath,'gsequery_filt', 
            embeddedpattern=True, tslocindex=1, returntype='returnlist'
        )
    if gsefilt_latest and len(gsefilt_latest)==1:
        gsefiltd = querydict(querypath=gsefilt_latest[0], splitdelim=' ')
        return gsefiltd
    else:
        print("Error: could not retrieve latest equery filt filepath! "
            +"Are there more than one latest file at the search directory?")
        return
