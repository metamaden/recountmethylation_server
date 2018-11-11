import ftplib
import datetime
import re
import pymongo
import socket
import struct
import sys
import time

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

def mongo_docs():
    """ Get entity dates from mongo db
        Requires:
            * local 'recount_methylation', connection info
        Returns
            * dates object (list): dates for entities in mongo db
    """
    timestamp = gettime_ntp() # current timestamp
    mongo_docs = [timestamp]
    client = pymongo.MongoClient('localhost', 27017)
    rmdb = client.recount_methylation
    grnidats = rmdb.gsm.idats.grn
    redidats = rmdb.gsm.idats.red
    soft = rmdb.gse.soft
    grncursor = grnidats.find({})
    for doc in grncursor:
          mongo_docs.append(doc)
    redcursor = redidats.find({})
    for doc in redcursor:
          mongo_docs.append(doc)
    softcursor = soft.find({})
    for doc in softcursor:
          mongo_docs.append(doc)
    return mongo_docs

def update_date(ftp_addresses, retries=1):
    """ Gets last updated date of entity (GSE soft or GSM idat) from FTP URL

        ftp_addresses: where to get update dates -- all should have same
            domain

        Return value: ftp error string or list of datetime objects with
            last updated date/timestamps corresponding respectively to ftp
            addresses. Times are UTC.
        Note: can add num. retries/attempts
    """
    firstftp = ftp_addresses[0]
    if not firstftp.startswith('ftp://'):
        raise RuntimeError("FTP address must begin with \"ftp://\".")
    firstftp_tokens = firstftp.split('/')
    try: 
        ftp = ftplib.FTP(firstftp_tokens[2])
        ftp.login()
    except ftplib.all_errors as e:
        return str(e)
    datetimes_to_return = []
    for currentftp in ftp_addresses:
        currentftp_tokens = currentftp.split('/')
        try:
            datetime_to_parse = ftp.sendcmd(
                    "MDTM /" + '/'.join(currentftp_tokens[3:])
                )
            datetimes_to_return.append(datetime.datetime.strptime(
                        datetime_to_parse[4:], "%Y%m%d%H%M%S")
                    )
        except ftplib.all_errors as e:
            datetimes_to_return.append(str(e))
    return datetimes_to_return

def diffs_mongo():
    """ Get date differences between mongo db and online 
        Returns
        * diffs object (list): entities and dates for diffs
    """
    md = mongo_docs() # mongo doc query with timestamp
    # list ftp addresses and past to get date
    mdftp = []
    for doc in md:
        try:
            if all (key in doc for key in ("ftp","date")):
                mdftp.append(doc)
            else:
                continue
        except:
            continue
    timestamp = gettime_ntp()
    diffs = [timestamp]
    ftplookup = [doc['ftp'] for doc in mdftp]
    ftp_onlinedate = update_date(ftplookup)
    # simplest check for err, assume valid ftplist len >1
    if len(ftp_onlinedate) == len(ftplookup):
        index = 0
        for doc in mdftp:
            doc['update_date'] = ftp_onlinedate[index]
            index += 1
        for doc in mdftp:
            if not doc['date'] == doc['update_date']:
                diffs.append(doc)
    else:
        diffs.append(ftp_onlinedate)
        diffs.append("'ftp_onlinedate' length inequality error")
    return diffs

""" Notes and tutorial section

client = pymongo.MongoClient('localhost', 27017)
rmdb = client.recount_methylation
grnidats = rmdb.gsm.idats.grn
redidats = rmdb.gsm.idats.red
soft = rmdb.gse.soft

# to clear db contents
# client.recount_methylation.command("dropDatabase")
# rmdb.drop_collection("gsm.idats.grn")
# rmdb.drop_collection("gsm.idats.red")
# rmdb.drop_collection("gse.soft")
# rmdb.list_collection_names() # should return empty list []

note: collections only created after a doc has been added to one
soft = rmdb.gse.soft
grn = rmdb.gsm.idats.grn
red = rmdb.gsm.idats.red

# add one doc to instantiate a collection
doc_add1 = {"gse": ""}
soft.insert_one(doc_add1).inserted_id
grn.insert_one(doc_add1).inserted_id
red.insert_one(doc_add1).inserted_id

rmdb.list_collection_names() 
# returns: ['gsm.idats.grn', 'gse.soft', 'gsm.idats.red']

doc_add2 = {"gsm": 'GSM1505330',
"filename": "GSM1505330_9376538060_R04C01_Grn.idat.gz",
"date": datetime.datetime(2014, 9, 16, 15, 15, 45),
"ftp":'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1505nnn/GSM1505330/suppl/
    GSM1505330_9376538060_R04C01_Grn.idat.gz'
}
grn.insert_one(doc_add2).inserted_id

"""



