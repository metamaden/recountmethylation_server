#!/usr/bin/env bash

# Installs, clones, and sets up dependencies for Recount Methylation server

# 1. Clone Recount Methylation Server and MetaSRA-pipeline fork
git clone https://github.com/metamaden/recount-methylation-server
# python3 dependencies
pip3 install pymongo
pip3 install celery

# 2. Setup MetaSRA-pipeline fork (via. metamaden or eventually pdxgx repo)
# Notes: 
# MetaSRA-pipeline uses Python2. 
# It has several dependencies, including: numpy, scipy, scikit-learn, 
# setuptools, and marisa-trie. 
# This fork modifies the code slightly, including support for an argument to 
# provide a filepath to store the mapped metadata, which is otherwise simply 
# output to console.
pip2 install numpy 
pip2 install scipy 
pip2 install scikit-learn
pip2 install setuptools
pip2 install marisa-trie
pip2 install dill
pip2 install nltk
# nltk requires punkt resource
python 2 
import nltk
nltk.download('punkt')
quit()
# Clone and setup MetaSRA-pipeline
git clone https://github.com/metamaden/MetaSRA-pipeline
cd MetaSRA-pipeline
git clone https://github.com/metamaden/bktree
cd setup_map_sra_to_ontology/
./setup.sh
# Notes on running MetaSRA-pipeline:
# setup should run through each script
# obo and lex terms should download
# a long list of term linkages should be printed out
# term linkages should then be inspected and checked
# setup then should finish quietly

# 3. Install celery

# 4. Install RabbitMQ

# 5. Install Mongo