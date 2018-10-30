import ftplib
import datetime
import re

def update_date(ftp_addresses):
    """ Gets last updated date of entity (GSE soft or GSM idat) from FTP URL

        ftp_addresses: where to get update dates -- all should have same
            domain

        Return value: ftp error string or list of datetime objects with
            last updated date/timestamps corresponding respectively to ftp
            addresses. Times are UTC.
    """
    if not ftp_address.startswith('ftp://'):
        raise RuntimeError("FTP address must begin with \"ftp://\".")
    ftp_tokens = ftp_address.split('/')
    try:
        ftp = ftplib.FTP(ftp_tokens[2])
        ftp.login()
    except ftplib.all_errors as e:
        return str(e)
    datetimes_to_return = []
    for ftp_address in ftp_addresses:
        try:
            datetime_to_parse = ftp.sendcmd(
                    "MDTM /" + '/'.join(ftp_tokens[3:])
                )
        except ftplib.all_errors as e:
            datetimes_to_return.append(str(e))
        else:
            datetimes_to_return.append(datetime.datetime.strptime(
                        datetime_to_parse[4:], "%Y%m%d%H%M%S")
                    )
    return datetimes_to_return
    