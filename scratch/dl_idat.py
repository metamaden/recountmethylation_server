import ftplib
import datetime
import os
import subprocess

def dl_idat(input_list,dldir):
    """
        Download idats, 
        Reads in either list of GSM IDs or ftp addresses
    
        Arguments
            * input list : list of GSM IDs
            * dldir : target directory for new downloads

        Returns 
            * Dictionary showing records, dates, and exit statuses of ftp calls
    """
    item = input_list[0]
    if not item.startswith('GSM'):
        raise RuntimeError("GSM IDs must begin with \"GSM\".")
    ftptoken_login = 'ftp.ncbi.nlm.nih.gov'
    if not os.path.exists(dldir):
        subprocess.call('mkdir '+dldir, shell=True)
    try:
        ftp = ftplib.FTP(ftptoken_login)
        loginstat = ftp.login()
    except ftplib.all_errors as e:
        return str(e)
    
    dldict = {}
    for id in input_list:
        dldict[id] = []
        id_ftptokens = ['ftp.ncbi.nlm.nih.gov','geo','samples',
        id[:-3]+'nnn',id,'suppl']
        id_ftpadd = '/'.join(id_ftptokens[1::])+'/'
        try:
            filenames = ftp.nlst(id_ftpadd)
            dldict[id].append([id,
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
                filename_towrite = file_tokens[-1]
                file_ftpadd = '/'.join(file_tokens[:-1])
                file_ftpadd = file_ftpadd+'/'+file_tokens[-1:][0]
                try:
                    filedl_estat = ftp.retrbinary("RETR /"+file_ftpadd,
                        open(os.path.join(dldir,
                            filename_towrite),'wb').write
                        )
                except ftplib.all_errors as efiledl:
                    filedl_estat = str(efiledl)
                dldict[id].append(
                    [id,
                    filename_towrite,
                    file_ftpadd,
                    filedate,
                    filedate_estat,
                    filedl_estat]
                    )
        except ftplib.error_temp as eid:
            dldict[id].append([id,id_ftpadd,str(eid)])
    return dldict
