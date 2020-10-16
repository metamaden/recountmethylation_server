#!/usr/bin/env bash

# start_server.sh 
# 
# Authors: Sean Maden, Abhi Nellore
# 
# Script to start vital server processes and run server.py.
# To be run from root dir for recount-methylation, above 
# recount-methylation-files and recount-methylation-server.
#
# Resources (state)
# * MongoDB (running): Supports file metadata DB (aka. 'RMDB').
# * RabbitMQ (running): Broker messaging program for celery.
# * Celery (running): Queue management software for Python.
# * SQLite (installed): Used for backend job queue reporting.

# 1. Start MongoDB
conpath="/home/metamaden/usr/local/bin/mongodb-linux-x86_64-4.0.4/bin/mongod"
dbpath="/home/metamaden/data/db"
dbconcmd=$conpath" --dbpath "$dbpath
eval $dbconcmd & # execute as background process

# 2. Start Broker
rabbitmq-server start

# 3. Start celery task queue manager
cd './recount-methylation-server/src'
celery worker -A gse_celerytask -l INFO & # execute as background process
cd -

# 4. Start Recount Methylation server
python3 ./recount-methylation-server/src/server.py & # executes backgrd. process
