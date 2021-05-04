#!/usr/bin/env python3

""" retry.py

Author: Sean Maden

Flexible functions to repeatedly attempt tasks for `recountmethylation_server`.

"""

import subprocess, glob, sys, os, re
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
settings.init()

def retry(script, complist, timelimit=180):
    """ retry
    
    Repeatedly attempt processes in `recountmethylation_server`. The indicated
    script should return a list of objects, to be compared with object provided
    in the `complist` argument. The process halts if either the time limit is
    reached or the percent provided by complistperc is achieved.

    Arguments:
    * script: Path to script for repeated attempts. Note, this must return 
        a list for comparison to the complist object.
    * complist: Comparison list for job status evaluation.
    * timelimit: Time limit in minutes (numeric, default 180min).

    Returns:

    Null, provides status updates over run.

    """

if __name__ == "__main__":
    print("Starting server.py..."); import subprocess, glob, sys, os, re
    sys.path.insert(0, os.path.join("recountmethylation_server","src"))
    import edirect_query, settings, argparse; settings.init()
    from edirect_query import gsm_query, gse_query, gsequery_filter  
    from utilities import gettime_ntp, getlatest_filepath, querydict
    from utilities import get_queryfilt_dict
    from gse_celerytask import gse_task; from random import shuffle
    gselist = [] # queue input, gse-based
    qstatlist = [] # job status object, also stored at sqlite db
    print("Getting timestamp...")
    run_timestamp = gettime_ntp() # pass this result to child functions
    # Parse the specified GSE ID.
    parser = argparse.ArgumentParser(description='Arguments for server.py')
    parser.add_argument("--gseid", type=str, required=False, default=None, 
        help='Option to enter valid GSE ID for immediate download.')
    args = parser.parse_args()

