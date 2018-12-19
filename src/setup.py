#!/usr/bin/env python3

import ftplib
import datetime
import os
import subprocess
import glob
import filecmp
import re
import gzip
import socket
import struct
import time
import tempfile
import atexit
import shutil
from itertools import chain
from celery import Celery
import sys
# sys.path.insert(0, os.path.join("recount-methylation-server","src"))
import pymongo

# recount methylation scripts

# clone vital repositories
# git clone https://github.com/metamaden/MetaSRA-pipeline
# cd MetaSRA-pipeline
# git clone https://github.com/metamaden/bktree
# cd -