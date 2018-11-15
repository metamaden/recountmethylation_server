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

def gettime_ntp(addr='0.uk.pool.ntp.org'):
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

def dl_idat(input_list, dest_dir='idats',temp_dir='temp_dir', retries=3, interval=.1, validate=True):
    """ Download idats, 
        Reads in either list of GSM IDs or ftp addresses
    
        Arguments
            * input list : list of GSM IDs
            * dest_dir : target directory for new downloads
            * temp_dir : temporary directory to store files
            * retries: number of times to retry a given file to download
            * interval: time (in s) before first retry
            * validate: compare most recently downloaded version with previous
                version of file after download and if new file is same, delete
                it and make note in dictionary
            ADD db connector to be able to make query in comment

        Returns 
            * Dictionary showing records, dates, and exit statuses of ftp calls
            OR error string over connection issues
    """
    timestamp = str(gettime_ntp())
    os.makedirs(dest_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    temp_dir_make = tempfile.mkdtemp(dir=temp_dir)
    atexit.register(shutil.rmtree, temp_dir_make)
    item = input_list[0]
    if not item.startswith('GSM'):
        raise RuntimeError("GSM IDs must begin with \"GSM\".")
    ftptoken_login = 'ftp.ncbi.nlm.nih.gov'
    try:
        ftp = ftplib.FTP(ftptoken_login)
        loginstat = ftp.login()
    except ftplib.all_errors as e:
        return str(e)
    client = pymongo.MongoClient('localhost', 27017) # mongodb connection
    dldict = {}
    for gsm_id in input_list:
        dldict[gsm_id] = []
        id_ftptokens = [
                'ftp.ncbi.nlm.nih.gov', 'geo', 'samples',
                gsm_id[:-3] + 'nnn', gsm_id, 'suppl'
            ]
        id_ftpadd = '/'.join(id_ftptokens[1::])+'/'
        retries_left = retries
        files_written = []
        filenames = []
        while retries_left:
            try:
                filenames = ftp.nlst(id_ftpadd)
                dldict[gsm_id].append([gsm_id,
                    id_ftpadd,
                    "success"]
                    ) 
                for file in filenames:
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
                            break
                        else:
                            filedate_estat = "new_date"
                            to_write = '.'.join([gsm_id, str(timestamp), 
                            file_tokens[-1]])
                            file_ftpadd = '/'.join(file_tokens[:-1])
                            file_ftpadd = file_ftpadd+'/'+file_tokens[-1:][0]
                            try:
                                with open(
                                        os.path.join(temp_dir_make, to_write), 'wb'
                                    ) as output_stream:
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
                                break
                            except ftplib.all_errors as efiledl:
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
                    except ftplib.all_errors as efiledate:
                        filedate_estat = str(efiledate)
                        filedate = "not_available"
                        dldict[gsm_id].append(
                                [gsm_id,
                                file,
                                filedate,
                                filedate_estat]
                                )
                        break
            except ftplib.error_temp as eid:
                if retries_left:
                    retries_left -= 1
                    print('ftplib error encountered, retries left = '
                        +retries_left)
                    time.sleep(interval)
                    continue
                else:
                    dldict[gsm_id].append([gsm_id, id_ftpadd, str(eid)])
                    break
        if validate:
            for gsm_id, file_written, index in files_written:
                gsms = glob.glob('.'.join([
                            os.path.join(dest_dir, gsm_id), '*', '.'.join(
                                            file_written.split('.')[2:]
                                    )
                        ]))
                most_recent = gsms.sort(
                        key=lambda x: int(x.split('.')[1])
                    )[-2:]
                if len(most_recent) == 2:
                    if filecmp.cmp(most_recent[0], most_recent[1]):
                        os.remove(most_recent[1])
                        # If filename is false, we found it was the same
                        dldict[gsm_id][index][1] = False
                else:
                    shutil.copyfile(
                        os.path.join(temp_dir_make,most_recent[0]), 
                        os.path.join(dest_dir,most_recent[0])
                        )
    return dldict


"""Example section

# data entries
gsm_id = 'GSM1505330'
file = 'geo/samples/GSM1505nnn/GSM1505330/suppl/
GSM1505330_9376538060_R04C01_Grn.idat.gz'
timestamp = '1541719636'
to_write = 'GSM1505330.1541719636.GSM1505330_9376538060_R04C01_Grn.idat.gz'
date1 = datetime.datetime(2014, 9, 16, 15, 15, 45)
date2 = datetime.datetime(2015, 9, 16, 15, 15, 45)

# mongo options and query
client = pymongo.MongoClient('localhost', 27017)
# make the db and collections
rmdb = client["recount_methylation"]
gsm = rmdb["gsm"]
idats = gsm["idats"]
grn = idats["grn"]
red = idats["red"]
client.list_database_names()
rmdb.list_collection_names() # returns: ['gsm.idats.grn']
# add a single doc
doc_add1 = {"gsm": 'GSM1505330',
"filename": "GSM1505330_9376538060_R04C01_Grn.idat","date": date1
}
doc_add2 = {"gsm": 'GSM1505330',
"filename": "GSM1505330_9376538060_R04C01_Grn.idat","date": date2
}
add1 = grn.insert_one(doc_add1).inserted_id
add2 = grn.insert_one(doc_add2).inserted_id

grn.find_one({'gsm' : gsm_id},{'date' : 1})
[print(d) for d in grn.find({'gsm' : gsm_id},{'date' : 1})]
for file in grn.find({'gsm' : gsm_id},{'date' : 1}):
    print(file)
datelist = []
[datelist.append(d['date']) for d in grn.find({'gsm' : gsm_id},{'date' : 1})]
dbdate = max(datelist)

datelist2 = []
[datelist2.append(d['date']) for d in red.find({'gsm' : gsm_id},{'date' : 1})]
[print(d) for d in red.find({'gsm' : gsm_id},{'date' : 1})]

dbdate = idat_mongo_query(gsm_id,file,client)

"""