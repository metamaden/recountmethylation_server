#!/usr/bin/env bash

# start_server.sh 
# Script to start vital server processes and run server.py.
# To be run from root dir for recount-methylation, above 
# recount-methylation-files and recount-methylation-server.
# Requisite processes for recount-methylation-server to run:
#       * mongodb: database for file metadata
#       * broker: messaging program for celery
#       * celery queue management software for python

# start mongodb
conpath="/home/metamaden/usr/local/bin/mongodb-linux-x86_64-4.0.4/bin/mongod"
dbpath="/home/metamaden/data/db"
dbconcmd=$conpath" --dbpath "$dbpath
eval $dbconcmd & # execute as background process

# start broker
brokerpath="/home/metamaden/usr/local/src/..."
eval $brokerpath &

# Start celery task queue manager
cd './recount-methylation-server/src'
celery worker -A gse_celerytask -l INFO & # execute as background process
cd -

# Start the server
python3 ./recount-methylation-server/src/server.py & # execute as bg process

# Notes and How-To's
# kill $! # kills last job
#