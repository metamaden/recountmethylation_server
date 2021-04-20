#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.join("recount-methylation-server","src"))
sys.path.insert(0, os.path.join("recount-methylation-analysis","src"))
from utilities import gettime_ntp, get_queryfilt_dict, getlatest_filepath
import settings
settings.init()

"""

Processing samples is straightforward with Recount Methylation. Note that these 
directions assume all samples were downloaded with 'server.py' and recorded in 
RMDB, as described in the notebook 'run_server'. It is also assumed 
sample metadata has been mapped using MetaSRA-pipeline, as detailed in the 
notebook 'preprocess_soft'. If this is not the case, please consult the relevant
documentation.

First, load the requisite functions from '~./recount-methylation-analysis/src'.

"""

from process_idats import expand_idats
from rsheet import rmdb_fpaths, compile_rsheet
from preprocess_mdat import preprocess_mdat, monitor_and_compile_mdata

"""

Initially, it is necessary to expand compressed idats that were downloaded, which
can be managed using the function 'expand_idats'. This step should be undertaken
before attempting to compile an rsheet for the first time.

"""

expand_idats()

"""

Next, compile the rsheet of downloaded samples. This calls RMDB for valid sample
files, filtering on samples with MetaSRA-pipeline mapped metadata available.

As a side-effect, new hard linkes (or 'hlink' files) to idats are created, and 
any old hlinks removed. This may cause the function to take up to several 
minutes to complete.

"""

gsmfpaths = rmdb_fpaths()
compile_rsheet(gsmfpaths, qcprint=True)


"""

With the rsheet compiled and fresh hard links to idats established, you may 
proceed to preprocess the idat samples. This will produce several new filetypes
under the 'mdata' subdirectory in the './recount-methylation-analysis/files/' 
folder. 

Note, to mitigate processing times, memory constraints, and other 
possible resource-related bottlenecks, samples are processed in subsets or
batches. These subsets are then compiled into large compilations tables in the
'~/mdata/compilations/' subdirectory. Note that compilations can become quite 
large (see below for details), and will vary in size depending on the number of
decimals for rounding specified in the 'mdat.R' script (defaults to '5').

Datatypes created include:
1. rawgrn : Raw Green channel intensities.
2. rawred : Raw Red channel intensities.
3. detp : Detection p-values table
3. rawbeta : Raw (unnormalized) Beta-values (derived from rawgrn and rawred)
2. noobbeta : Noob-normalized Beta-values

Additional subdirectories in '~/mdata' include:
* compilations : Folder to contain all large compilation tables collating subset
tables of each datatype.
* logs : Log dictionary objects detailing exitstatuses and standard errors for
processes generated from using 'preprocess_mdat()'.

"""

preprocess_mdat(nsampproc=10, nprocmax=10, statint=2)


"""

To start forming compilations, ensure at least one subset table is available 
for each datatype (detp, noobbeta, rawbeta, rawgrn, and rawred). It is 
recommended you run this function in a new screen or tab, concurrent with 
'preprocess_mdat()', as it is intended to monitor and collate new subset tables 
as they are generated.

"""

monitor_and_compile_mdata()

"""

Be advised, compilation tables can become quite large. We suggest ensuring you 
have enough disk space for 5 tables varying from 50-100 Gb, depending on the 
rounding factor you specified in 'mdat.R' (defaults to '5').

"""

















