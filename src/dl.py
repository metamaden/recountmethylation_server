import ftplib
import datetime
import os
import ftplib
import datetime
import os
import subprocess
import glob
import filecmp
import pymongo
import socket
import struct
import sys
import time
import tempfile
import atexit
import shutil

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

def soft_mongo_date(gse,filename,client):
    """ Get date(s) from mongo soft subcollection

        Arguments
            * gse : valid GSE ID
            * filename : name of soft file
            * client : mongodb client connection

        Returns
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
    """ Get date(s) from mongo idat subcollections

        Arguments
            * gsm_id : valid sample GSM ID
            * filename : name of array idat file, inc. 'grn' or 'red'
            * client : mongodb client connection

        Returns
            * mongo_date_list : list of resultant date(s) from query, or empty 
                                list if no docs detected
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

def dl_idat(input_list, dest_dir='idats', temp_dir='temp', 
    retries_connection=3, retries_files=3, interval=.1, validate=True):
    """ Download idats, 
        Reads in either list of GSM IDs or ftp addresses
    
        Arguments
            * input list : list of GSM IDs
            * dest_dir : target directory for new downloads
            * temp_dir : temporary directory to store files
            * retries_connection: num ftp connection retries
            # retries_files : num retry attempts on sample files
            * interval: time (in sec) before first retry
            * validate: compare most recently downloaded version with previous
                version of file after download and if new file is same, delete
                it and make note in dictionary
        Returns 
            * Dictionary showing records, dates, and exit statuses of ftp calls
            OR error string over connection issues
    """
    # timestamp = str(gettime_ntp())
    timestamp = 'timestamp'
    os.makedirs(dest_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    temp_dir_make = tempfile.mkdtemp(dir=temp_dir)
    atexit.register(shutil.rmtree, temp_dir_make)
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
                time.sleep(interval)
                continue
            else:
                print('connection retries exhausted, returning...')
                return str(e)
    client = pymongo.MongoClient('localhost', 27017) # mongodb connection
    dldict = {}
    for gsm_id in input_list:
        print('Starting GSM: '+gsm_id)
        dldict[gsm_id] = []
        id_ftptokens = [
                'ftp.ncbi.nlm.nih.gov', 'geo', 'samples',
                gsm_id[:-3] + 'nnn', gsm_id, 'suppl'
            ]
        id_ftpadd = '/'.join(id_ftptokens[1::])+'/'
        files_written = []
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
                                    time.sleep(interval)
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
                            time.sleep(interval)
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
                time.sleep(interval)
                continue
            else:
                print('File retries exhausted. Breaking...')
                dldict[gsm_id].append([gsm_id, id_ftpadd, str(eid)])
                break
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
                    dldict[gsm_id][index].append(False)
                else:
                    print("Downloaded file is new, moving to dest_dir...")
                    shutil.move(file_written, os.path.join(
                            dest_dir, os.path.basename(file_written))
                        )
                    dldict[gsm_id][index].append(True)
                    dldict[gsm_id][index][2] = os.path.join(
                                dest_dir, os.path.basename(file_written)
                            )
            else:
                print("Downloaded file is new, moving...")
                shutil.move(file_written, os.path.join(
                            dest_dir, os.path.basename(file_written))
                        )
                dldict[gsm_id][index].append(True)
                dldict[gsm_id][index][2] = os.path.join(
                                dest_dir, os.path.basename(file_written)
                            )
    return dldict

def dl_soft(gse_list, dest_dir='gse_soft', temp_dir='temp', retries_connection=3,
    retries_files=3, interval=.1, validate=True):
    """ Download idats, 
        Reads in either list of GSM IDs or ftp addresses
    
        Arguments
            * input list : list of GSM IDs
            * dest_dir : target directory for new downloads
            * temp_dir : temporary directory to store files
            * retries_connection : num retries for connection, on ftp err.
            * retries_files : num retries for given file, on ftp err.
            * interval: time (in s) before first retry
            * validate: compare most recently downloaded version with previous
                version of file after download and if new file is same, delete
                it and make note in dictionary
        Returns 
            * Dictionary showing records, dates, and exit statuses of ftp calls
            OR error string over connection issues
    """
    timestamp = str(gettime_ntp())
    os.makedirs(dest_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    temp_dir_make = tempfile.mkdtemp(dir=temp_dir)
    atexit.register(shutil.rmtree, temp_dir_make)
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
                time.sleep(interval)
                continue
            else:
                print('connection retries exhausted, returning...')
                return str(e)
    client = pymongo.MongoClient('localhost', 27017) # mongodb connection
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
                                        (gse, to_write,len(dldict[gse]) - 1)
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
                                    +retries_left_files)
                                time.sleep(interval)
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
                            +retries_left_files)
                        time.sleep(interval)
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
                        +retries_left_files)
                    time.sleep(interval)
                    continue
                else:
                    print('file retries exhausted, breaking..')
                    dldict[gse].append([gse, id_ftpadd, str(eid)])
                    break
    if validate:
        print('commencing file validation...')
        for file in files_written:
            gse = file[0]
            file_written = file[1]
            print(gse)
            print(file_written)
            gses = glob.glob('.'.join([
                        os.path.join(dest_dir, gse), '*', '.'.join(
                                        file_written.split('.')[2:]
                                )
                    ]))
            if gses:
                if len(gses)>1:
                    most_recent = gses.sort(
                        key=lambda x: int(x.split('.')[1])
                    )
                    most_recent = most_recent[-1]
                else:
                    most_recent = gses[0]
                if filecmp.cmp(most_recent, file_written):
                    print('identical file found in dest_dir, removing...')
                    os.remove(file_written)
                    # If filename is false, we found it was the same
                    dldict[gse].append(False)
                else:
                    print('new file detected in temp_dir, moving to '
                        +'dest_dir...')
                    dldict[gse].append(True)
                    dldict[gse][index][2] = os.path.join(
                                dest_dir, os.path.basename(file_written)
                            )
                    shutil.move(file_written, os.path.join(
                            dest_dir, os.path.basename(file_written))
                        )
            else:
                print('new file detected in temp_dir, moving to dest_dir..')
                dldict[gse].append(True)
                dldict[gse][index][2] = os.path.join(
                                dest_dir, os.path.basename(file_written)
                            )
                shutil.move(file_written, os.path.join(
                            dest_dir, os.path.basename(file_written))
                        )
            continue
    return dldict

def update_rmdb(ddidat,ddsoft,host='localhost',port=27017):
    """ Update recount-methylation database with new docs
        Arguments
            * ddidat : download dictionary from dl_idats
            * ddsoft : download dicitonary from dl_soft
        Returns
            * statusdict : result list (1 = new doc added, 0 = no doc added) 
    """
    statusdict = {}
    client = pymongo.MongoClient('localhost', 27017)
    rmdb = client.recount_methylation
    if ddidat:
        statusdict['ddidat'] = ""
        gsmc = rmdb.gsm
        idatsc = gsmc.idats
        for gsmkey in list(ddidat.keys()):
            lvals = ddidat[gsmkey]
            for lval in lvals:
                if lval[-1]==True:
                    new_idatdoc = {"gsmid":lval[0],
                        "ftpaddress":lval[1],
                        "filepath":lval[2],
                        "exitstatus":lval[3],
                        "date":lval[4]
                        }
                    idatsc.insert_one(new_idatdoc)
                    statusdict['ddidat'] = 1
                else:
                    statusdict['ddidat'] = 0
    if ddsoft:
        statusdict['ddsoft'] = ""
        gsec = rmdb.gse
        softc = gsec.soft
        for gsekey in list(ddsoft.keys()):
            lvals = ddsoft[gsekey]
            if lvals[-1]==True:
                lsoft = lvals[1]
                new_softdoc = {"gseid":lval[0],
                        "ftpaddress":lval[1],
                        "filepath":lval[2],
                        "exitstatus":lval[3],
                        "date":lval[4]
                        }
                softc.insert_one(new_softdoc)
                statusdict['ddsoft'] = 1
            else:
                statusdict['ddsoft'] = 0
    return statusdict