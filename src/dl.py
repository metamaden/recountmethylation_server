#!/usr/bin/env python3

""" dl.py

    Authors: Sean Maden, Abhi Nellore
    
    Functions to manage download of idats and experiment soft files from GEO.
    Notes:
        * Validation for each download function checks new files against 
            existing files in the corresponding destination files directory.
        * Downloads are launched from the job definition for the celery job 
            queue manager. Each queued job is based around a valid GSE id.
    Functions:
        * soft_mongo_date: grab latest update date for a soft file from the 
            Recount Methylation Mongo db (or 'RMDB'). 
        * idat_mongo_date: grab latest update date for an idat file from RMDB.
        * dl_idat: Download and validate idat files.
        * dl_soft: Download and validate soft files.
"""

import ftplib
import datetime
import os
import sys
import subprocess
import glob
import fnmatch
import filecmp
import pymongo
import time
import tempfile
import shutil
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, getlatest_filepath
import settings
settings.init()

def soft_mongo_date(gse, filename, client):
    """ soft_mongo_date
        Get date(s) from mongo soft subcollection.
        Arguments:
            * gse (str) : Valid GSE ID.
            * filename (str) : Name of a valid soft file.
            * client (conn.) : A client connection to RMDB.
        Returns:
            * mongo_date_list : list of resultant date(s) from query, or empty 
                list if no docs detected
    """
    mongo_date_list = []
    rmdb = client.recount_methylation
    gsec = rmdb.gse
    softc = gsec.soft  
    mongo_date_list.append(
        d['date'] for d in softc.find(
        {'gse' : gse},
        {'date' : 1}))
    return mongo_date_list

def idat_mongo_date(gsm_id,filename,client):
    """ idat_mongo_date
        Get date(s) from mongo idat subcollections.
        Arguments:
            * gsm_id (str) : A valid sample GSM ID.
            * filename (str) : Name of valid array idat file, inc. 'grn' or 
                'red'.
            * client (conn.) : A client connection to RMDB.
        Returns:
            * mongo_date_list (list) : list of resultant date(s) from query, or 
                empty list if no docs detected
    """
    mongo_date_list = []
    rmdb = client.recount_methylation
    gsmc = rmdb.gsm
    idatc = gsmc.idats  
    if 'grn' in filename:
        mongo_date_list.append(
            d['date'] for d in idatc.grn.find(
            {'gsm' : gsm_id},
            {'date' : 1}))
    else:
        mongo_date_list.append(
            d['date'] for d in idatc.red.find(
            {'gsm' : gsm_id},
            {'date' : 1}))
    return mongo_date_list

def dl_idat(input_list, retries_connection=3, retries_files=3, interval_con=.1, 
    interval_file=.01, validate=True, timestamp=gettime_ntp()):
    """ dl_idat
        Download idats, reading in either list of GSM IDs or ftp addresses.
        Arguments
            * input list (list, required) : A list of valid GSM IDs.
            * retries_connection (int) : Number of ftp connection retries 
                allowed.
            # retries_files : Number of retry attempts allowed for sample file
                downloads.
            * interval_con (float) : Time (in seconds) to sleep before retrying 
                a database connection.
            * interval_file (float) : Time (in seconds) to sleep before retrying 
                a file connection. 
            * validate (Bool.): Validate new files against existing idats?
            * timestamp (str) : An NTP timestamp for versioning.
        Returns 
            * dldict (dictionary) : Records, dates, and exit statuses of ftp 
                calls, OR error string over connection issues. Downloads and 
                moves new and validated files as side effect. 
    """
    idatspath = settings.idatspath
    temppath = settings.temppath
    os.makedirs(idatspath, exist_ok=True)
    os.makedirs(temppath, exist_ok=True)
    temp_dir_make = tempfile.mkdtemp(dir=temppath)
    item = input_list[0]
    if not item.startswith('GSM'):
        raise RuntimeError("GSM IDs must begin with \"GSM\".")
    ftptoken_login = 'ftp.ncbi.nlm.nih.gov'
    retries_left_connection = retries_connection
    while retries_left_connection:
        print('trying ftp connection')
        try:
            ftp = ftplib.FTP(ftptoken_login)
            loginstat = ftp.login()
            print('connection successful, continuing...')
            break
        except ftplib.all_errors as e:
            if retries_left_connection:
                retries_left_connection -= 1
                print('continuing with connection retries left = '
                    +str(retries_left_connection))
                time.sleep(interval_con)
                continue
            else:
                print('connection retries exhausted, returning...')
                return str(e)
    # mongodb connection
    client = pymongo.MongoClient(settings.rmdbhost, settings.rmdbport) 
    dldict = {}
    files_written = []
    for gsm_id in input_list:
        print('Starting GSM: '+gsm_id)
        dldict[gsm_id] = []
        id_ftptokens = [
                'ftp.ncbi.nlm.nih.gov', 'geo', 'samples',
                gsm_id[:-3] + 'nnn', gsm_id, 'suppl'
            ]
        id_ftpadd = '/'.join(id_ftptokens[1::])+'/'
        filenames = []
        retries_left_files = retries_files
        try:
            filenames = ftp.nlst(id_ftpadd)
            if len(filenames)>0:
                filestr = '; '.join(str(e) for e in filenames)
                print("files found: "+filestr)
                dldict[gsm_id].append([gsm_id,
                    id_ftpadd,
                    "connection success, valid num idats found"]
                    )
                print("Idat filenames detected for "+gsm_id+", continuing...") 
                for file in filenames:
                    print("Beginning iteration for file: "+file)
                    filedate = ""
                    filedate_estat = ""
                    filedl_estat = ""
                    file_tokens = file.split('/')
                    try:
                        filedate = ftp.sendcmd("MDTM /" + '/'.join(file_tokens))
                        filedate = datetime.datetime.strptime(filedate[4:],
                            "%Y%m%d%H%M%S")
                        mongo_date = idat_mongo_date(gsm_id,file,client)
                        if filedate in mongo_date:
                            filedate_estat = "same_as_local_date"
                            dldict[gsm_id].append(
                                [gsm_id,
                                file,
                                filedate,
                                filedate_estat]
                                )
                            print('Online date same as local date. Breaking..')
                            break
                        else:
                            filedate_estat = "new_date"
                            to_write = os.path.join(
                                    temp_dir_make,
                                    '.'.join([gsm_id, str(timestamp), 
                                    file_tokens[-1]])
                                )
                            file_ftpadd = '/'.join(file_tokens[:-1])
                            file_ftpadd = file_ftpadd+'/'+file_tokens[-1:][0]
                            print('Attempting file download, for file: '
                                +file)
                            try:
                                with open(to_write, 'wb') as output_stream:
                                    filedl_estat = ftp.retrbinary(
                                            "RETR /"+file_ftpadd,
                                            output_stream.write
                                        )
                                dldict[gsm_id].append(
                                        [gsm_id,
                                        file_ftpadd,
                                        to_write,
                                        filedl_estat,
                                        filedate,
                                        filedate_estat]
                                    )
                                if '226 Transfer complete' in filedl_estat:
                                    files_written.append(
                                            (gsm_id, to_write, 
                                                len(dldict[gsm_id]) - 1)
                                        )
                                print("File successfully downloaded. "
                                    +"Continuing...")
                                continue
                            except ftplib.all_errors as efiledl:
                                if retries_left_files:
                                    retries_left_files -= 1
                                    print('ftp file dl error, retries left = '
                                    +str(retries_left_files))
                                    time.sleep(interval_file)
                                    continue
                                else:
                                    print('File retries exhausted. Breaking...')
                                    filedl_estat = str(efiledl)
                                    dldict[gsm_id].append(
                                            [gsm_id,
                                            file_ftpadd,
                                            to_write,
                                            filedl_estat,
                                            filedate,
                                            filedate_estat]
                                        )
                                    break
                            break
                        break    
                    except ftplib.all_errors as efiledate:
                        if retries_left_files:
                            retries_left_files -= 1
                            print('ftplib file date error, retries left = '
                            +str(retries_left_files))
                            time.sleep(interval_file)
                            continue
                        else:
                            print('File retries exhausted. Breaking...')
                            filedate_estat = str(efiledate)
                            filedate = "not_available"
                            dldict[gsm_id].append(
                                    [gsm_id,
                                    file,
                                    filedate,
                                    filedate_estat]
                                    )
                            break
                    continue 
            else:
                dldict[gsm_id].append([gsm_id,
                    "no files at ftp address"]
                    )
                break
        except ftplib.error_temp as eid:
            if retries_left_files:
                retries_left_files -= 1
                print('ftplib filenames error, retries left = '
                    +str(retries_left_files))
                time.sleep(interval_file)
                continue
            else:
                print('File retries exhausted. Breaking...')
                dldict[gsm_id].append([gsm_id, id_ftpadd, str(eid)])
                break
    if validate:
        print("Validating downloaded files...")
        for gsm_id, file_written, index in files_written:
            print("file written is "+file_written)
            filestr = os.path.basename(file_written).split('.')[2::]
            filestr = str('.'.join(filestr))
            print('filestr written : '+filestr)
            print('dir to search latest: '+idatspath)
            gsmidat_latest = getlatest_filepath(idatspath, filestr,
                    embeddedpattern=True, returntype='returnlist', tslocindex=1
                )
            print('gsm latest: '+str(gsmidat_latest))
            if gsmidat_latest:
                gsmidat_latest = gsmidat_latest[0]
                print('cmp result: '
                        +str(filecmp.cmp(gsmidat_latest, file_written))
                    )
                if filecmp.cmp(gsmidat_latest, file_written):
                    print("Downloaded file is same as recent file. Removing...")
                    os.remove(file_written)
                    # If filename is false, we found it was the same
                    dldict[gsm_id][index].append(False)
                else:
                    print("Downloaded file is new, moving to idatspath...")
                    shutil.move(file_written, os.path.join(
                            idatspath, os.path.basename(file_written))
                        )
                    dldict[gsm_id][index].append(True)
                    dldict[gsm_id][index][2] = os.path.join(
                                idatspath, os.path.basename(file_written)
                            )
            else:
                print("Downloaded file is new, moving...")
                shutil.move(file_written, os.path.join(
                            idatspath, os.path.basename(file_written))
                        )
                dldict[gsm_id][index].append(True)
                dldict[gsm_id][index][2] = os.path.join(
                                idatspath, os.path.basename(file_written)
                            )
        shutil.rmtree(temp_dir_make)
    return dldict

def dl_soft(gse_list=[], retries_connection=3, retries_files=3, interval_con=.1, 
    interval_file=.01, validate=True, timestamp=gettime_ntp()):
    """ dl_soft
        Download GSE soft file(s). Accepts either a list of GSM IDs or ftp 
        addresses.
        Arguments:
            * gse_list (list, required) : A list of valid GSE id(s).
            * retries_connection (int) : Number of ftp connection retries 
                allowed. 
            * retries_files : Number of retry attempts allowed for sample file
                downloads. 
            * interval_con (float) : Time (in seconds) to sleep before retrying 
                a database connection. 
            * interval_file (float) : Time (in seconds) to sleep before retrying 
                a file connection. 
            * validate (Bool.): Validate new files against existing idats?
            * timestamp (str) : An NTP timestamp for versioning.     
        Returns: 
            * Dictionary showing records, dates, and exit statuses of ftp calls
                OR error string over connection issues
    """
    gsesoftpath = settings.gsesoftpath
    temppath = settings.temppath
    os.makedirs(gsesoftpath, exist_ok=True)
    os.makedirs(temppath, exist_ok=True)
    temp_dir_make = tempfile.mkdtemp(dir=temppath)
    item = gse_list[0]
    if not item.startswith('GSE'):
        raise RuntimeError("GSE IDs must begin with \"GSE\".")
    ftptoken_login = 'ftp.ncbi.nlm.nih.gov'
    retries_left_connection = retries_connection
    while retries_left_connection:
        print('trying ftp connection')
        try:
            ftp = ftplib.FTP(ftptoken_login)
            loginstat = ftp.login()
            print('connection successful, continuing...')
            break
        except ftplib.all_errors as e:
            if retries_left_connection:
                retries_left_connection -= 1
                print('continuing with connection retries left = '
                    +str(retries_left_connection))
                time.sleep(interval_con)
                continue
            else:
                print('connection retries exhausted, returning...')
                return str(e)
    # mongodb connection
    client = pymongo.MongoClient(settings.rmdbhost, settings.rmdbport) 
    dldict = {}
    print('beginning iterations over gse list...')
    for gse in gse_list:
        print('beginning download for gse: '+gse)
        retries_left_files = retries_files 
        dldict[gse] = []
        files_written = []
        filenames = []
        # tokens for soft file ftp address
        id_ftptokens = [
                'ftp.ncbi.nlm.nih.gov', 'geo', 'series',
                gse[:-3] + 'nnn', gse, 'soft'
            ]
        id_ftpadd = '/'.join(id_ftptokens[1::])+'/'
        while retries_left_files:
            try:
                filenames = ftp.nlst(id_ftpadd)
                # filter for only soft file names
                file = list(filter(lambda x:'family.soft' in x,filenames))[0]
                dldict[gse].append([gse,
                    id_ftpadd,
                    "success"]
                    )
                filedate = ""
                filedate_estat = ""
                filedl_estat = ""
                file_tokens = file.split('/')
                try:
                    print('getting date from '+'/'.join(file_tokens))
                    filedate = ftp.sendcmd("MDTM /" + '/'.join(file_tokens))
                    filedate = datetime.datetime.strptime(filedate[4:],
                        "%Y%m%d%H%M%S")
                    mongo_date = soft_mongo_date(gse,file,client)
                    if filedate in mongo_date:
                        print('online  date same as local date,'
                            +'breaking...')
                        filedate_estat = "same_as_local_date"
                        dldict[gse].append(
                            [gse,
                            file,
                            filedate,
                            filedate_estat]
                            )
                        break
                    else:
                        print('new online date found, continuing...')
                        filedate_estat = "new_date"
                        to_write = os.path.join(
                                temp_dir_make,
                                '.'.join([gse, timestamp, 
                                file_tokens[-1]])
                            )
                        file_ftpadd = '/'.join(file_tokens[:-1])
                        file_ftpadd = file_ftpadd+'/'+file_tokens[-1:][0]
                        try:
                            print('downloading soft from '+file_ftpadd)
                            with open(to_write, 'wb') as output_stream:
                                filedl_estat = ftp.retrbinary(
                                        "RETR /"+file_ftpadd,
                                        output_stream.write
                                    )
                            dldict[gse].append(
                                    [gse,
                                    file_ftpadd,
                                    to_write,
                                    filedl_estat,
                                    filedate,
                                    filedate_estat]
                                )
                            if '226 Transfer complete' in filedl_estat:
                                files_written.append(
                                        (gse, to_write, len(dldict[gse]) - 1)
                                    )
                            print('total files written = '
                                +str(len(files_written)))
                            print('soft transfer successful for '
                                    +to_write+', breaking...')
                            break
                        except ftplib.all_errors as efiledl:
                            print('file download error from '+file_ftpadd)
                            if retries_left_files:
                                retries_left_files -= 1
                                print('continuing with file retries left ='
                                    +str(retries_left_files))
                                time.sleep(interval_file)
                                continue
                            else:
                                print('file retries exhausted, breaking..')
                                filedl_estat = str(efiledl)
                                dldict[gse].append(
                                        [gse,
                                        file_ftpadd,
                                        to_write,
                                        filedl_estat,
                                        filedate,
                                        filedate_estat]
                                    )
                                break
                except ftplib.all_errors as efiledate:
                    print('error getting date from '+'/'.join(file_tokens))
                    if retries_left_files:
                        retries_left_files -= 1
                        print('continuing with file retries left = '
                            +str(retries_left_files))
                        time.sleep(interval_file)
                        continue
                    else:
                        print('file retries exhausted, breaking..')
                        filedate_estat = str(efiledate)
                        filedate = "not_available"
                        dldict[gse].append(
                                [gse,
                                file,
                                filedate,
                                filedate_estat]
                                )
                        break
            except ftplib.error_temp as eid:
                print('error making ftp connection to '+id_ftpadd)
                if retries_left_files:
                    retries_left_connection -= 1
                    print('ftplib error encountered, file retries left = '
                        +str(retries_left_files))
                    time.sleep(interval_file)
                    continue
                else:
                    print('file retries exhausted, breaking..')
                    dldict[gse].append([gse, id_ftpadd, str(eid)])
                    break
    if validate:
        print('commencing file validation...')
        for gse, new_filepath, index in files_written:
            filestr = os.path.basename(new_filepath).split('.')[0]
            gsesoft_latest = getlatest_filepath(gsesoftpath,filestr)
            if gsesoft_latest and not gsesoft_latest == 0:
                if filecmp.cmp(gsesoft_latest, new_filepath):
                    print('identical file found in dest_dir, removing...')
                    dldict[gse].append(False)
                    os.remove(new_filepath)
                else:
                    print('new file detected in temp_dir, moving to '
                        +'dest_dir...')
                    dldict[gse].append(True)
                    dldict[gse][index][2] = os.path.join(
                                gsesoftpath, os.path.basename(new_filepath)
                            )
                    shutil.move(new_filepath, os.path.join(
                            dest_dir, os.path.basename(new_filepath))
                        )
            else:
                print('new file detected in temp_dir, moving to dest_dir..')
                dldict[gse].append(True)
                dldict[gse][index][2] = os.path.join(
                                gsesoftpath, os.path.basename(new_filepath)
                            )
                shutil.move(new_filepath, os.path.join(
                            gsesoftpath, os.path.basename(new_filepath))
                        )
            continue
        shutil.rmtree(temp_dir_make)
    return dldict
