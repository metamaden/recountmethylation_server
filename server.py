#!/usr/bin/env python3
"""
server.py
Sean Maden

Script coordinating server functionality for Slow Revolt (codename for a
generalized framework for maintaining a self-updating processed version of
archival data). This script is currently tailored for the recount-methylation
project.

Script:
 * uses bash utilities scripts in parallel
 * polls for updates to E-direct (NCBI tool) query returned, on designated
   interval
 * situational response to discrepancies detected in query
 * runs metadata mapping pipeline using utilities and MetaSRA-pipeline
 * populates a MongoDB instace with updates to the database

[insert license info here]
"""

import multiprocessing
import subprocess
import os
import sys
import time
import datetime
import dateparser # use to parse GEO record timestamps

# For polling, investigate subclassing threading.Thread to run polling on a separate thread
# This poll should add to the multiprocessing pool's jobs
# To do this successfully, you should probably declare your multiprocessing pool
# outside the run_idat_jobs function.
# In polling, you need to check if studies have been updated too
# So the polling class (and it is a class) is updating the mongodb with new
# metadata for studies on the server, at the same time diffing to detect
# changes, and also adding to the queue for processing

def get_diffs(cdt_path,ndt_path):
	"""
	
	Retrieve sample diffs for for job queue update
	Get status of data and metadata downloads
	cdt_path : filepath to current dirtree 'cdt' 
	ndt_path : filepath to new dirtree 'ndt'
	Currently expect same format for cdt and ndt objects
	todo: Handle options to automate cdt/ndt creation if files don't exist

	"""
	sqadd_idat=[]
	sqadd_metadata=[]
	cdt_list = []
	ndt_list = []

	cdt_lines = [line.rstrip('\n') for line in open(cdt_path)]
	for el in cdt_lines:
		cdt_list += [el.split('\t')[1::]]
	ndt_lines = [line.rstrip('\n') for line in open(ndt_path)]
	for el in ndt_lines:
		ndt_list += [el.split('\t')[1::]]

	# check cdt logfiles
	for e in cdt_list:
		if os.path.exists(os.path.join('GSE', str(e[0]))):
			for s in e[1::]:
				if not os.path.exists(os.path.join('GSE', str(e[0]), str(s),str("{0}_idat-transfer_SUCCESS".format(s)))):
					sqadd_idat += s
				if not os.path.exists(os.path.join('GSE', str(e[0]), str(s),str("{0}_msrap-metadata_SUCCESS".format(s)))):
					sqadd_metadata += s
		else:
			sqadd_idat += e[1::]
			sqadd_metadata += e[1::]

	# get diffs between cdt and ndt
	for ein in ndt_list:
		enamei = ein[0]
		flc = [e for sublist in cdt_list for e in sublist] # flat list
		if not enamei in flc and not enamei in sqadd_idat:
			sqadd_idat += ein[1::]
		if not enamei in flc and not enamei in sqadd_metadata:
			sqadd_metadata += ein[1::]
		for s in ein[1::]:
			if not s in flc and not s in sqadd_idat:
				sqadd_idat += s
			if not s in flc and not s in sqadd_metadata:
				sqadd_metadata += s

	print("Found {0} idat and {1} metadata diffs.".format(len(sqadd_idat),len(sqadd_metadata)))
	sqadd_dt = datetime.date.fromtimestamp(time.time())
	sqadd_list = [sqadd_idat,sqadd_metadata,sqadd_dt]
	return sqadd_list

def run_idat_jobs(temp_dir, num_processes=50):
	""" 

	Run new idat download jobs in parallel
	Uses Pools and process size to run batch download jobs

	"""
	gsm_dl_file = open("./tempfiles/gsmArray", "r")
    with open(os.path.join(temp_dir, 'gsmArray')) as gsm_array_stream:
        dl_list = gsm_array_stream.read().split()

	# make list of lists of len num_processes up to end of inputted array
	n_last_slice = len(dl_list) % num_processes # remainder is len of last joblist item
	n_list_slice = len(dl_list) // num_processes # use floor, increment up to lastjob
	joblist = []
	# iterate over items in dl_list, appending ids as lists to run as batches
	if n_list_slice > 0:
		for si in range(0, n_list_slice * num_processes - 1, num_processes):
			dlslice = dl_list[si:si+num_processes]
			print('interval start: ' + si, file=sys.stderr)
			joblist.append(dlslice)
	joblist.append(dl_list[(len(dl_list)-(1+n_last_slice)):(len(dl_list)-1)])
	
	# parallelize idat downloads
	pool = multiprocessing.Pool(processes=num_processes)
	# iterate over job batches, parallelizing downloads within batch
	for ji in joblist[0:(len(joblist)-1)]:
		# ji a list of dl jobs incremented as specified in num_processes
        # review https://stackoverflow.com/questions/19924104/python-multiprocessing-handling-child-errors-in-parent
        # to robustify if child processes fail; handle failure.
        # if it fails, add to end of queue. keep log. 
        # if K failures are encountered (make this eg a command line 
        # parameter), log failure and proceed
		print("Running jobs: ", ji, file=sys.stderr)
		pool.imap_unordered(lambda x: subprocess.check_call(
                ['dl_gsm_idat.sh', gsm_name],
                    executable='/bin/bash', ji
            ))
		print("Finished jobs: ", ji, file=sys.stderr)
		
#------
# MAIN
#------
def main(working_dir, temp_dir):
	"""
	
	Script to run on call from cl.

	"""
    import tempfile
    temp_dir = tempfile.mkdtemp(dir=temp_dir)
    import atexit
    import shutil
    atexit.register(shutil.rmtree, temp_dir)
    tally_file = os.path.join(temp_dir, 'all_gse_newtally')
	try:
        subprocess.check_call(['remeth_server_utilities.sh'], 
        	executable='/bin/bash')
    except CalledProcessError:
        # Handle failure
        raise
	# parallelize idat downloads
	run_idat_jobs()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
                    description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter
                )
    parser.add_argument('--working-dir', '-w', type=str, required=True,
                        help='directory in which to store all downloads and '
                             'processed outputs')
    parser.add_argument('--temp-dir', '-t', type=str, required=True,
                        help='directory in which to store temporary files')
    parser.add_argument('--db', type=str, required=True,
                        help='Mongo database file')
    args = parser.parse_args()
    main(args.working_dir, args.temp_dir)
