import ftplib
import datetime
import os
import subprocess

def dl_idat(input_list, dldir, retries=3, interval=.1):
    """ Download idats, 
        Reads in either list of GSM IDs or ftp addresses
    
        Arguments
            * input list : list of GSM IDs
            * dldir : target directory for new downloads
            * retries: number of times to retry a given file to download
            * interval: time (in s) before first retry

        Returns 
            * Dictionary showing records, dates, and exit statuses of ftp calls
            OR error string over connection issues
    """
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
                        filedate_estat = "success"
                    except ftplib.all_errors as efiledate:
                        filedate_estat = str(efiledate)
                        filedate = "not_available"
                    to_write = file_tokens[-1]
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
                        filename_towrite,
                        file_ftpadd,
                        filedate,
                        filedate_estat,
                        filedl_estat]
                    )
            except ftplib.error_temp as eid:
                if retries_left:
                    retries_left -= 1
                    time.sleep(interval)
                    continue
                dldict[gsm_id].append([gsm_id, id_ftpadd, str(eid)])
    return dldict
