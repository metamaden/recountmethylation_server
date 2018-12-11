#!/usr/bin/env bash

# start_server.sh 
# Script to start vital server processes and run server.py
# 1. start celery
# 2. start mongodb

# start mongodb
mongodpath="./usr/local/bin/mongodb-linux-x86_64-4.0.4/bin/mongod"
dbpath="./home/metamaden/data/db"
cmd=$mongodpath+" --dbpath "+$dbpath
eval $cmd

# start celery task queue manager
celery worker -A ./recount-methylation-server/src/gse_celerytask -l INFO

# run server.py
python3 ./recount-methylation-server/src/server.py