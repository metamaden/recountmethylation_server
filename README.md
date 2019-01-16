# recount-methylation-server
Server code for hosting and maintaining Recount-Methylation database.

## Instructions to Set Up and Run Recount Methylation Server and Database Instance

This sections explains how to run the server and handle a local Recount Methylation instance.

### Details About Server Functioning

A Recount Methylation instance includes directories for files (at 'recount-methylation-files') and server source code (at 'recount-methylation-server/src'). It is highly recommended you do not modify the file directory names or paths to files directories, as these are handled and referenced automatically by the server.

To version files in a Recount Methylation instance, we append [NTP](https://en.wikipedia.org/wiki/Network_Time_Protocol) timestamps to file names, along with sample (GSM) or experiment (GSE) IDs and the original name of the file that was retrieved from GEO. NTP timestamps reflect universal time when code was run, and these are used to detect latest file versions by the server. For this reason, please do not modify these timestamps or downloaded file names directly!

You may see a single timestamp in a filename (e.g. for idat files in 'idats' directory, soft files in 'gse_soft' or 'gsm_soft' directories, etc.) or multiple timestamps (e.g. in 'gsm_msrap_outfiles' or 'gsm_json' directory). In the latter cases, the first timestamp reflects the time the file was processed, and the second is the timestamp on the sample GSM soft file that was run. This enables tracking what version of a sample GSM soft file was processed to generate a derived JSON or MetaSRA-pipeline mapped metadata file.

![alt text](https://github.com/metamaden/recount-methylation-server/blob/master/server_workflow.jpg "Recount Methylation Server Process")

### Installing the Server and Other Dependencies

You can clone the latest server repository from GitHub. We also recommend cloning our [fork of MetaSRA-pipeline](https://github.com/metamaden/MetaSRA-pipeline), which has been slightly modified to run with the server. Both repositories should be cloned to the same directory, called 'recount-methylation'. On first run, the server will create a new subdirectory and tree for downloaded files here, called 'recount-methylation-files.' Before attempting to run the server, please install and check all server dependencies, as detailed below.

### Dependencies and Resources Used by The Server

The server requires use of a [MongoDB](https://www.mongodb.com/) document database (for handling the documents and downloads metadata), [Celery](http://www.celeryproject.org/) distributed task queue manager (Python software for server job queue handling), [RabbitMQ](https://www.rabbitmq.com/) broker (for messaging required by Celery), and [EDirect](https://dataguide.nlm.nih.gov/edirect/install.html) query utility from NCBI (for scheduled queries to GEO). The server runs with Python 3##, while MetaSRA-pipeline runs with Python 2##, and various library dependencies are required for each (see affiliated ReadMe documentation for each). Some bash dependencies are required, including [GNU Screen](https://www.gnu.org/software/screen/), to make full use of source code. Other dependencies include R with the ['jsonlite'](https://cran.r-project.org/web/packages/jsonlite/index.html) and ['minfi'](https://bioconductor.org/packages/release/bioc/html/minfi.html) libraries installed. 

Please ensure your system has all of these programs installed and working before attempting to run the server. You can check the script at 'recount-methylation-server/src/settings.py', which specifies global variables and values used throughout the source code, to ensure paths to these various dependency resources are correct and working. Please refrain from modifying any default file directory names, paths, values, etc. in 'settings.py' unnecessarily, as this can impede regular server function.

### Steps to Run the Server with Screen

We highly recommend running the server in the background using GNU screen. For instance, from the top level of the 'recount-methylation' directory, the server can be run from command line using:

```{bash}
screen -S rmserverpy # creates new screen
python3 ./recount-methylation-server/src/server.py # runs server in new screen
# Ctrl-A + Ctrl-D escapes the screen
```

Please wait while the server runs. It may take several days, depending on the system and connection, to complete the download for the compilation of interest (e.g. >35,000 samples and experiments with HM450 idat files available).

### Steps to Process Recount Methylation Files in Python 3

Once you are ready and have enough files downloaded, you can begin to process downloaded files. From the top level of 'recount-methylation' directory, you should now see a 'recount-methylation-files' subdirectory containing array and experiment metadata files. For the following, we will assume your Recount Methylation instance includes all or most files from experiments that use HM450k and have available sample idats in supplment. 

We recommend starting an interactive Python 3 session for file processing, with the following steps to set up the session:

```{python}
import os
import sys
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
from process_soft import expand_soft, extract_gsm_soft, gsm_soft2json, msrap_screens
```

#### Processing GSE Soft files

GSE soft files contain metadata for samples. To extract sample-level (that is, GSM ID-level) metadata from the GSE soft files, first expand the available GSE soft files by running the following in python. Note, you can specify whether to remove compressed versions of successfully expanded files setting the 'rmcompressed' argument to 'True':

```{python}
expand_soft()
```

Now extract the sample metadata from the expanded experiment soft files, using:

```{python}
extract_gsm_soft()
```

Next, in preparation for mapping the sample metadata using MetaSRA-pipeline, convert extracted GSM soft files to JSON format with the following:

```{python}
gsm_soft2json()
```

This will convert the sample soft data to valid JSON format, while expanding most commonly nested metadata terms, by running an R function that uses the 'jsonlite' library.

Finally, after GSM soft data has been extracted and converted to JSON format, you are ready to map samples using MetaSRA-pipeline. Please consult the pipeline's affiliated ReadMe and setup instructions, to ensure it is set up and running properly. To expedite sample mapping, we implement multiple screens in the background automatically using the msrap_screens() function. 

It is important to note the optimum screen count and samples per screen will depend greatly on the system and available compute resources. Often, it is worth experimenting with different numbers of parallel screen sessions to determine a configuration that maps samples as quickly as possible. 

As an example, to automatically deploy N = 30 screens each running an instance of MetaSRA-pipeline to processively map 200 JSON sample files, use:

```{python}
msrap_screens(nscreensi=200,nmaxscreens=30,qcprint=True)
```

After a few minutes, you should start seeing new ".msrapout" files generated in the "gsm_msrap_outfiles" directory. 

## Processing and Analyzing Recount Methylation Samples and Files

### Processing Overview

The preprocess workflow utilizes information on the latest versions of files, stored in the Recount Methylation MongoDB (RMDB), a database of document compilations containing file metadata. RMDB file records are validated by searching the root files directory ('./recount-methylation-files'), and metadata and idat files are paired under their corresponding shared GSM ID. Finally, this information is used to form a data sheet including valid file names and paths, and a compressed/'flattened' version of the MetaSRA-pipeline mapped sample metadata. 

![alt text](https://github.com/metamaden/recount-methylation-server/blob/master/preprocess_workflow.jpg "Recount Methylation Server Process")

### Processing Tutorial

TBD


