import ftplib
import datetime
import os
import subprocess
import glob
import filecmp

def gettime_ntp(addr='time.nist.gov'):
    # http://code.activestate.com/recipes/117211-simple-very-sntp-client/
    import socket
    import struct
    import sys
    import time
    TIME1970 = 2208988800L      # Thanks to F.Lundh
    client = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    data = '\x1b' + 47 * '\0'
    client.sendto(data, (addr, 123))
    data, address = client.recvfrom( 1024 )
    if data:
        t = struct.unpack( '!12I', data )[10]
        t -= TIME1970
    return t

def dl_idat(input_list, dldir, retries=3, interval=.1, validate=True):
    """ Download idats, 
        Reads in either list of GSM IDs or ftp addresses
    
        Arguments
            * input list : list of GSM IDs
            * dldir : target directory for new downloads
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
    item = input_list[0]
    if not item.startswith('GSM'):
        raise RuntimeError("GSM IDs must begin with \"GSM\".")
    ftptoken_login = 'ftp.ncbi.nlm.nih.gov'
    os.makedirs(dldir, exist_ok=True)
    try:
        ftp = ftplib.FTP(ftptoken_login)
        loginstat = ftp.login()
    except ftplib.all_errors as e:
        return str(e)
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
                        # Check if filedate is the same as in db, and if it is
                        # break out, or indicate in return value that this was
                        # the case.
                        filedate_estat = "success"
                    except ftplib.all_errors as efiledate:
                        filedate_estat = str(efiledate)
                        filedate = "not_available"
                    to_write = '.'.join([gsm_id, timestamp, file_tokens[-1]])
                    file_ftpadd = '/'.join(file_tokens[:-1])
                    file_ftpadd = file_ftpadd+'/'+file_tokens[-1:][0]
                    try:
                        with open(
                                os.path.join(dldir, to_write), 'wb'
                            ) as output_stream:
                            filedl_estat = ftp.retrbinary(
                                    "RETR /"+file_ftpadd,
                                    output_stream.write
                                )
                    except ftplib.all_errors as efiledl:
                        filedl_estat = str(efiledl)
                    dldict[gsm_id].append(
                        [gsm_id,
                        to_write,
                        file_ftpadd,
                        filedate,
                        filedate_estat,
                        filedl_estat]
                    )
                    if '226 Transfer complete' in filedl_estat:
                        files_written.append(
                                (gsm_id, to_write, len(dldict[gsm_id]) - 1)
                            )
            except ftplib.error_temp as eid:
                if retries_left:
                    retries_left -= 1
                    time.sleep(interval)
                    continue
                dldict[gsm_id].append([gsm_id, id_ftpadd, str(eid)])
            else:
                break
        if validate:
            for gsm_id, file_written, index in files_written:
                gsms = glob.glob('.'.join([
                            os.path.join(dldir, gsm_id), '*', '.'.join(
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
    return dldict
