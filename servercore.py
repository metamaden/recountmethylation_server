#!/usr/bin/env python3

"""

Info:
Script to coordinate server functionality for Slow Revolt: Recount Methylation pilot, instance
Utilizes bash utilities scripts in parallel
Polls for updates to E-direct query returned, on designated interval
Situational response to discrepancies detected in query
Runs metadata mapping pipeline using utilities and MetaSRA-pipeline
Populates a MongoDB instace with updates to the database

Author:
Sean Maden

"""

#--------------
# DEPENDENCIES
#--------------
# for parallelizing jobs
import multiprocessing # for parallelizing jobs
# for passing shell arguments from python
import subprocess # for making shell calls and processing shell scripts
import shlex # for passing formatted commands to shell

#---------
# CLASSES
#---------

#-----------
# FUNCTIONS
#-----------

def get_gsm_idat(gsm_name):
	"""

	Function to call bash function+script to download an individual GSM's data files
	Need to import subprocess and shlex first
	Argument is sample GSM ID, which is entered by pools in parallelization

	"""
	subprocess.call(shlex.split('bash dl_gsm_idat.sh '+gsm_name))

def run_idat_jobs(nproc=50):
	"""

	Run new idat download jobs in parallel
	Uses Pools and process size to run batch download jobs

	"""

	gsm_dl_file = open("./tempfiles/gsmArray", "r")
	dl_list=gsm_dl_file.read().split()

	# make list of lists of len nproc up to end of inputted array
	n_last_slice=len(dl_list)%nproc # remainder is len of last joblist item
	n_list_slice=len(dl_list)//nproc # use floor, increment up to lastjob
	joblist=[]
	# iterate over items in dl_list, appending ids as lists to run as batches
	if n_list_slice>0:	
		slice_interval=range(0,(n_list_slice*nproc-1),nproc)
		for si in slice_interval:
			dlslice=dl_list[si:si+nproc]
			print(si)
			joblist.append(dlslice)
	joblist.append(dl_list[(len(dl_list)-(1+n_last_slice)):(len(dl_list)-1)])
	
	# parallelize idat downloads
	pool = multiprocessing.Pool(processes=nproc)
	# iterate over job batches, parallelizing downloads within batch
	for ji in joblist[0:(len(joblist)-1)]:
		# ji a list of dl jobs incremented as specified in nproc
		print("Running jobs :",ji)
		pool.imap_unordered(get_gsm_idat,ji)
		print("Finished jobs:",ji)
		
#------
# MAIN
#------
def main():
	"""
	
	Script to run on call from cl.

	"""

	# pass this on a duration interval (eg. '1 day', '5 days', '1 week', etc.)
	subprocess.call(shlex.split("mkdir -p tempfiles")) # tempfiles dir contains current edir queryfile
	open('./tempfiles/all_gse_newtally', 'w')
	#subprocess.call(shlex.split("echo '' > ./tempfiles/all_gse_newtally")) # instantiate file for new query data
	subprocess.call(shlex.split('bash remeth_server_utilities.sh')) # defines and runs bash functions
	# parallelize idat downloads
	run_idat_jobs()

if __name__ == "__main__":
    main() # call primary script




