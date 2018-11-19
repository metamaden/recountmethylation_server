import filecmp
from itertools import chain

def querydict(query):
    """ Process edirect query result into dictionary of GSEs
        Arguments
            * query (file) : file of edirect results (format: 1 GSE per newline)
        Returns
            * query (dictionary) : dictionary (keys = GSEs, vals = GSM lists)
    """
    querylines = [line.rstrip('\n') for line in open(query)]
    querylist = []
    querydict = {}
    for line in querylines:
        querylist.append([line.split('\t')[1::]])
    for gselist in querylist:
        gselist = list(chain.from_iterable(gselist))
        gsekey = list(filter(lambda x:'GSE' in x, gselist))[0]
        gsmlist = list(filter(lambda x:'GSM' in x, gselist))
        querydict[gsekey] = gsmlist
    return querydict

def equery_filter(gsequery='gsequery',gsmquery='gsmquery'):
    """ Prepare an edirect query file
        Filter GSE query on GSM query membership
        Arguments
            * gsequery (str): filename for GSE query file from edirect.sh
            * gsmqeury (str): filename for GSM query file from edirect.sh
        Returns
            * gsequeryfiltered (list): filtered GSE query file
    """
    gsed = querydict(query=gsequery)
    gsmlines = [line.rstrip('\n') for line in open(gsmquery)]
    gsmlist = [] # list of valid gsms with idat suppfiles
    for line in gsmlines:
        gsmlist.append(line.split('\t')[1::][0])
    gsedfilt = gsed
    for key in gsed.keys():
        ival = gsed[key]
        for gsm in ival:
            if not gsm in gsmlist:
                gsedfilt[key] = gsedfilt[key].remove(gsm)
        return gsedfilt


def edirect_diffs(query1,query2):
    """ Compares two edirect query results, returning diffs
        Arguments
            * query1 : first edirect query, filename
            * query2 : second edirect query, filename
        Returns
            * diffs object (list): list of GSEs with diffs
    """
    difflist = []
    qd1 = querydict(query1) # eg. for first value: qd1[list(qd1.keys())[0]]
    qd2 = querydict(query2)
    # if gse doesn't exist in qd2
    for key in qd1:
        evalkey = ''
        if key in list(qd2.keys()):
            if not qd1[key]==qd2[key]:
                for item in qd1[key]:
                    if not item in qd2[key]:
                        evalkey = 'notalike'
                for item in qd2[key]:
                    if not item in qd1[key]:
                        evalkey = 'notalike'
        else:
            evalkey = 'notpresent'
        if (evalkey == 'notalike') or (evalkey == 'notpresent'):
            difflist.append(key)
    for key in qd2:
        if not key in list(qd1.keys()):
            difflist.append(key)
    return difflist


""" examples

query1 = 'newquery'
query2 = 'querydiff'
edirect_diffs(query1,query2)

import edirect_query
from itertools import chain

gse_filename='gsequery'
gsm_filename='gsmquery'

from equery_utilities import querydict
from equery_utilities import equery_filter

gsed1 = equery_filter()
gsed2 = querydict(query='gsequery')

gsed1[list(gsed1.keys())[0]]
gsed2[list(gsed1.keys())[0]]
testl = []
for key in gsed1.keys():
    if gsed1[key]==gsed2[key]:
        testl.append("same")
    else:
        testl.append("not_alike")

# filter the gse query on valid samples
gsed = edirect_query.querydict(query=gse_filename)
gsmquery = gsm_filename
gsmlines = [line.rstrip('\n') for line in open(gsmquery)]
gsmlist = []
for line in gsmlines:
    gsmlist.append(line.split('\t')[1::][0])
for key in gsed.keys():
    ival = gsed[key]
    for gsm in ival:
        if not gsm in gsmlist:
            gsed[key] = gsed[key].remove(gsm)
    
"""