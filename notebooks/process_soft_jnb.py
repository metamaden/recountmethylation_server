#!/usr/bin/env python3

"""

Authors: Sean Maden, Abhi Nellore

To process the experiment soft files, Recount Methylation begins by downloading
the soft metadata file for each experiment. Each downloaded soft file thus has a 
corresponding GSE experiment ID. 

Experiment soft file downloads are managed along with GSM/sample idat downloads
from the celery job queue manager. GSE soft files are downloaded first, followed
by the valid GSM sample idats. Files are downloaded to a temp dir then validated 
against existing files in the 'gse_soft' subdir of 'recount-methylation-files'. 
After successful download, followed by validation showing a new soft file or 
that an existing soft file was updated, a call to add a new GSE soft file record 
to RMDB is made. This record includes the relevant file and download metadata.

The filenames of experiment soft files, like other Recount Methylation files,
are versioned using NTP timestamps in the filename, along with the experiment id
and the name of the original file that was downloaded from the GEO FTP call. 

"""

import process_soft

"""

To get from the GSE soft file to the mapped GSM metadata, several steps are 
taken: 1. expand GSE soft files; 2. extract valid GSM subsections from GSE soft
files; 3. convert extracted GSM soft data to JSON format; 4. pass the valid 
JSON-formatted GSM soft metadata to MetaSRA-pipeline , using the specified
fork of this application.

"""

"""

First, expand GSE soft files. This unpacks downloaded experiment soft files with
the option of removing compressed files after expansion:

"""

process_soft.expand_soft()



"""

Then extract GSM sample metadata. This will create new sample files in the 
gsm_soft files directory:

"""

process_soft.extract_gsm_soft()

"""

Next, convert GSM sample metadata from the soft file format to valid JSON format. 
This function passes sample soft files with an R function to convert to JSON.

"""

process_soft.gsm_soft2json()

"""

Finally map the JSON-formatted GSM metadata using MetaSRA-pipeline. 

Note, because the mapping process can take quite a while, depending on file size 
and format, we highly recommend deploying the pipeline across several background 
or detached screen sessions using the wrapper function. 

You can experiment with the number of screens and samples per screen that works 
best on your system. 

"""

process_soft.msrap_screens(nscreensi=200, nmaxscreens=35, qcprint=True)














