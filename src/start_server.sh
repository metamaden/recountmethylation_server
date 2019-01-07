#!/usr/bin/env bash

# start_server.sh 
# Script to start vital server processes and run server.py.
# To be run from root dir for recount-methylation, above 
# recount-methylation-files and recount-methylation-server.
#
# Requisite processes for recount-methylation-server to run properly:
# * mongodb: database for file metadata ('RMDB')
# * broker: messaging program for celery (e.g. RabbitMQ)
# * celery: queue management software for python

# 1. Start MongoDB
conpath="/home/metamaden/usr/local/bin/mongodb-linux-x86_64-4.0.4/bin/mongod"
dbpath="/home/metamaden/data/db"
dbconcmd=$conpath" --dbpath "$dbpath
eval $dbconcmd & # execute as background process

# 2. Start Broker
# brokerpath="/home/metamaden/usr/local/src/..."
# eval $brokerpath &
rabbitmq-server start

# 3. Start celery task queue manager
cd './recount-methylation-server/src'
celery worker -A gse_celerytask -l INFO & # execute as background process
cd -

# 4. Start Recount Methylation server
python3 ./recount-methylation-server/src/server.py & # executes backgrd. process
