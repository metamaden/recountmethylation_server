#!/usr/bin/env python3

""" settings.py
    Sets global variables for recount methylation, namely for filepaths and 
    platformid. Settings global objects are accessible by importing settings 
    then running 'settings.init()' and declaring a global object (e.g. 
    'settings.filesdir').
    Notes:
        * It is highly recommended you do not change any settings unnecessarily,
            including files directory names and paths.
        * You will need to change certain sections in settings.py, namely under
            the '[resource path]' section, depending on your system.
        * Consult all setup and README documents for details and troubleshooting
            instructions.
"""

import os

def init():
    # [platform]
    global platformid
    platformid = 'GPL13534'

    # [filedirectories]
    global filesdir
    global tempdir
    global gsesoftdir
    global gsmsoftdir
    global idatsdir
    global equerydir
    global gsmjsondir
    global gsmmsrapoutdir
    global sheetsdir
    global temppath
    global gsesoftpath
    global gsmsoftpath
    global idatspath
    global equerypath
    global gsmjsonpath
    global gsmmsrapoutpath
    global sheetspath
    filesdir = 'recount-methylation-files'
    tempdir = 'temp'
    gsesoftdir = 'gse_soft'
    gsmsoftdir = 'gsm_soft'
    idatsdir = 'idats'
    equerydir = 'equery'
    gsmjsondir = 'gsm_json'
    gsmmsrapoutdir = 'gsm_msrap_outfiles'
    sheetsdir = 'sheetfiles'
    temppath = os.path.join(filesdir, tempdir)
    gsesoftpath = os.path.join(filesdir, gsesoftdir)
    gsmsoftpath = os.path.join(filesdir, gsmsoftdir)
    idatspath = os.path.join(filesdir, idatsdir)
    equerypath = os.path.join(filesdir, equerydir)
    gsmjsonpath = os.path.join(filesdir, gsmjsondir)
    gsmmsrapoutpath = os.path.join(filesdir, gsmmsrapoutdir)
    sheetspath = os.path.join(filesdir, sheetsdir)

    # [serverdirectories]
    global serverdir
    global srcdir
    global serversrcpath
    global psoftscriptname
    global soft2jsonscriptname
    global psoftscriptpath
    global s2jscriptpath
    serverdir = 'recount-methylation-server'
    srcdir = 'src'
    psoftscriptname = 'process_soft.py'
    soft2jsonscriptname = 'soft2json.R'
    serversrcpath = os.path.join(serverdir, srcdir)
    psoftscriptpath = os.path.join(serversrcpath, psoftscriptname)
    s2jscriptpath = os.path.join(serversrcpath, soft2jsonscriptname)

    # [metasrapipeline_directories]
    global pipelinedir
    global msrapdir
    global runscriptname
    global msrapfnstem
    global msrappath
    global msraprunscriptpath
    pipelinedir = '.'
    msrapdir = 'MetaSRA-pipeline'
    runscriptname = 'run_pipeline.py'
    msrapfnstem = 'msrapout'
    msrappath = os.path.join(pipelinedir, msrapdir)
    msraprunscriptpath = os.path.join(msrappath, runscriptname)

    # [rmdb connection]
    global rmdbhost
    global rmdbport
    rmdbhost = 'localhost'
    rmdbport = 27017

    # [resource paths]
    global mongoconnpath
    global mongodbpath
    global runmongocmd
    global runrabbitmqcmd
    global runcelerycmd
    mongoconnpath = '/home/metamaden/usr/local/bin/mongodb-linux-x86_64-4.0.4/bin/mongod'
    mongodbpath = '/home/metamaden/data/db'
    runmongocmd = ' '.join(['eval',mongoconnpath,'--dbpath',mongodbpath])
    runrabbitmqcmd = 'rabbitmq-server start'
    runcelerycmd = ' | '.join(['cd '+serversrcpath,
            'celery worker -A gse_celerytask -l INFO &',
            'cd -'
        ]
    )

    # [launch server.py]
    global launchserverpycmd
    launchserverpycmd = 'python3 ./recount-methylation-server/src/server.py'

    # [misc]
    global gsepatt
    global softallpatt
    global compsoftpatt
    global expsoftpatt
    global jsonfnpattern
    global msrapoutfnpattern
    global gsequerystr
    global gsmquerystr
    gsepatt = '^GSE.*'
    softallpatt = '.*\.soft.*'
    compsoftpatt = '.*\.gz$'
    expsoftpatt = '.*\.soft$'
    jsonfnpattern = '.*json$'
    msrapoutfnpattern = '.*\.msrapout$'
    gsequerystr = 'gse_edirectquery' 
    gsmquerystr = 'gsm_edirectquery'

    # global eqfiltdict
    # eqfiltdict = utilities.get_queryfilt_dict()

