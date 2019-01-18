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
    global temppath
    global gsesoftpath
    global gsmsoftpath
    global idatspath
    global equerypath
    global gsmjsonpath
    global gsmmsrapoutpath
    filesdir = 'recount-methylation-files'
    tempdir = 'temp'
    gsesoftdir = 'gse_soft'
    gsmsoftdir = 'gsm_soft'
    idatsdir = 'idats'
    equerydir = 'equery'
    gsmjsondir = 'gsm_json'
    gsmmsrapoutdir = 'gsm_msrap_outfiles'
    temppath = os.path.join(filesdir, tempdir)
    gsesoftpath = os.path.join(filesdir, gsesoftdir)
    gsmsoftpath = os.path.join(filesdir, gsmsoftdir)
    idatspath = os.path.join(filesdir, idatsdir)
    equerypath = os.path.join(filesdir, equerydir)
    gsmjsonpath = os.path.join(filesdir, gsmjsondir)
    gsmmsrapoutpath = os.path.join(filesdir, gsmmsrapoutdir)

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

    # [analysis files, directories, and values]
    global analysisdir
    global analysisfilesdir
    global analysissrcdir
    global analysissrcpath
    global analysisfilespath
    global sheetfnstem
    global sheetsdir
    global sheetspath
    global msraptablesdir
    global msraptablespath
    analysisdir = 'recount-methylation-analysis'
    analysisfilesdir = 'files'
    analysissrcdir = 'src'
    analysissrcpath = os.path.join(analysisdir, analysissrcdir)
    analysisfilespath = os.path.join(analysisdir, analysisfilesdir)
    sheetfnstem = 'rsheet'
    sheetsdir = 'rsheets'
    sheetspath = os.path.join(analysisfilespath, sheetsdir)
    msraptablesdir = 'msraptables'
    msraptablespath = os.path.join(analysisfilespath, msraptablesdir)

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
    global grnidat_expcatch
    global redidat_expcatch
    gsepatt = '^GSE.*'
    softallpatt = '.*\.soft.*'
    compsoftpatt = '.*\.gz$'
    expsoftpatt = '.*\.soft$'
    jsonfnpattern = '.*json$'
    msrapoutfnpattern = '.*\.msrapout$'
    gsequerystr = 'gse_edirectquery' 
    gsmquerystr = 'gsm_edirectquery'
    grnidat_expcatch = '.*Grn\.idat$'
    redidat_expcatch = '.*Red\.idat$'

    # global eqfiltdict
    # eqfiltdict = utilities.get_queryfilt_dict()

