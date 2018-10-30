#!/usr/bin/env python3

"""

Notes on Mongodb format
* one 'recount_methylation' db contains:
** 2 collections, one for 'gse' and one for 'gsm'
** each record in gse contains:
*** gse soft file doc
** each record in 'gsm' includes docs:
*** idat records (filenames, lastmod dates)
*** metadata records (filename, lastmod date)


"""

import csv
import re
import multiprocessing
import subprocess
import os
import sys
import time
import datetime
import dateparser
import pymongo




def main():
	# Setup Mongodb
	# subprocess.check_call(['brew services start mongodb'],executable='/bin/bash')
	myclient = pymongo.MongoClient("mongodb://localhost:27017/")
	# Check for remethdb and collections 
	if not 'recount_methylation' in myclient.list_database_names():
		remethdb = myclient["recount_methylation"]
	if not 'gse' remethdb.list_collection_names():
		mycolgse = remethdb["gse"]
	if not 'gsm' remethdb.list_collection_names():
		mycolgsm = remethdb["gsm"] # sample collection
		mycolidat = mycolgsm["idat"] # idat collection
		mycolgrn = mycolidat["grn"] # green channel idat collection
		mycolred = mycolidat["red"] # green channel idat collection
		mycolmd = mycolgsm["metadata"] # metadata collection

		# read in datidat files
	tableidat = []
	sdk = ['_id','GSE','filename','ftpaddress','exitstatus','date']
	for file in os.listdir(os.path.join('datidat'))[1::]:
		filename = file
		with open(os.path.join('datidat',filename), encoding = "ISO-8859-1") as tsvfile:
			readtsv = csv.reader(tsvfile, delimiter=' ')
			for row in readtsv:
				sdd = {}
				itemct = 0
				for item in row:
					sdd[sdk[itemct]] = item
					itemct += 1
				tableidat += [sdd]

		# modify idat info
	grnnew = [i for i in tableidat if re.search('Grn.idat',i['filename'])] # add check for existing ids, use sep. insert and update docs
	colgrnid = []
	for x in mycolgrn.find({},{ "_id": 1}):
	  colgrnid += [x]
	colgrnidlist = [doc['_id'] for doc in colgrnid]
	grnupdate = [i for i in grnadd if i['_id'] in colgrnidlist]
	grnadd = [i for i in grnadd if not i['_id'] in colgrnidlist]
	if len(grnadd) > 0:
		if len(grnadd) > 1:
			try:
				mycolgrn.insert_many(grnadd) # raises bulk write error - perhaps check if ids redundant before adding
			except BulkWriteError as bwe:
				print(bwe.details)
		else:
			try:
				mycolgrn.insert_one(grnadd[0])
			except BulkWriteError as bwe:
				print(bwe.details)

	"""

	On update
	Get diffs and update
	Iterate like this:
	if len(grnupdate) > 0:
		if len(grnupdate) > 1:
			mycolgrn.update_many(grnupdate) # raises bulk write error - perhaps check if ids redundant before adding
		else:
			try:
				mycolgrn.update_one(grnupdate[0])
			except BulkWriteError as bwe:
				print(bwe.details)

	"""

		"""

		Example query:
	myquery = { "filename": grnadd[0]['filename'] }
	mydoc = mycolgrn.find(myquery)
	for x in mydoc:
		print(x)

	myquery = {'_id':{}}
	mydoc = mycolgrn.find(myquery)
	for x in mydoc:
		print(x)

		"""

	# update all docs meeting criteria of query
myquery = { "address": { "$regex": "^S" } }
newvalues = { "$set": { "name": "Minnie" } }
x = mycol.update_many(myquery, newvalues)


if __name__ == "__main__":
	main()
