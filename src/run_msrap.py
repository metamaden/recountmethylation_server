#!/usr/bin/env python3

""" process_soft.py

    Authors: Sean Maden, Abhi Nellore
    
    Run the MetaSRA-pipeline. Included functions allow for batch analysis of 
    sample JSON files in the form of composite files that are passed to the
    pipeline. The outputs include the new composite JSON files, the new
    composite files with pipeline mappings, and the sample/GSM-level pipeline
    outputs.
    
    Notes:

    Functions:
        * write_cjson: Writes the composite JSON files to be used in the 
            pipeline.
        * get_gsm_outputs: Writes GSM-specific metadata files using composite
            JSON files and pipeline output files with matching timestamps.
        * run_msrap_compjson: Performs the full workflow, including composite 
            JSON generation from the filtered JSON files (write_cjson function),
            pipeline metadata mapping, and sample/GSM-level metadata output
            writes (get_gsm_outputs function).
"""

import os, sys, re, gzip, shutil, subprocess, filecmp, tempfile, pickle, time
from datetime import datetime; from random import shuffle
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from utilities import gettime_ntp, getlatest_filepath, get_queryfilt_dict
from utilities import monitor_processes; import settings; settings.init()

def write_cjson(jffnv, newfilefn = "cjson", 
    jsonfiltpath = settings.gsmjsonfiltpath, 
    msrap_destpath = settings.gsmmsrapoutpath, tempdname = "cjsontemp",
    ts = gettime_ntp()):
    """ Write a composite JSON file with multiple samples

    Arguments: 
        * jffnv (list): Vector of filtered JSON filenames.
        * newfilefn (str): File name stem of new file to write
        * msrap_destpath (str): Path to MetaSRA-pipeline output files.
        * jsonfiltpath (str): Path to filtered GSM JSON files.
        * tempdname (str): Name of dir, at jsonfiltpath, to contain composite 
            JSON 
        files.
        * ts (int): Timestamp of output and input files.

    Returns:
        * Path to new composite JSON file.

    """
    temppath_read = os.path.join(jsonfiltpath)
    if not os.path.exists(temppath_read):
        os.makedirs(temppath_read)
    temppath_write = os.path.join(msrap_destpath, tempdname)
    if not os.path.exists(temppath_write):
        os.makedirs(temppath_write)
    ll = []; fnl = []
    for fn in jffnv:
        fpath = os.path.join(jsonfiltpath, fn)
        if os.path.exists(fpath):
            with open(fpath, "r") as openjson:
                linesi = openjson.readlines()
                if len(linesi) > 0:
                    ll.append(linesi)
                    fnl.append(fpath)
    newfn = ".".join([newfilefn, ts])
    wite_fpath = os.path.join(temppath_write, newfn)
    if len(ll) > 0:
        print("Read data for " + str(len(ll)) + " files. Writing data")
        lform = []
        with open(wite_fpath, "w") as opencj:
            opencj.write("[\n") # first line
            for fi, file in enumerate(ll):
                ld = []
                for line in file:
                    if line == "}\n":
                        opencj.write("\t{\n"); ld = ld[1::]
                        jfname = os.path.basename(fnl[fi])
                        gsmid = '"'+jfname.split(".")[1]+'"'
                        lpath = ":".join(['"gsm"', gsmid])
                        opencj.write("\t\t" + lpath + ",\n") # sample id
                        for ii, ldi in enumerate(ld):
                            lf = ldi.split(":")
                            if lf[0] == '  !Sample_source_name_ch1':
                                lf = ['source'] + lf[1::]
                            elif lf[0] == '  !Sample_title':
                                lf = ['title'] + lf[1::]
                            else:
                                lf = lf[1::]
                            lf = ['"' + i + '"' for i in lf];lf = ':'.join(lf)
                            if ii == len(ld) - 1:
                                lf = lf + "\n" # comma for values before last
                            else:
                                lf = lf + ",\n"
                            opencj.write("\t\t" + lf)
                        if fi == len(ll) - 1:
                            opencj.write("\t}\n") # no comma for final entry
                        else:
                            opencj.write("\t},\n")
                    else:
                        ldi = line
                        ldi = ldi.replace(']', '');ldi = ldi.replace('[', '')
                        ldi = ldi.replace('"', '');ldi = ldi.replace('\n', '')
                        ldi = ldi.replace(',', '')
                        if not ldi == "":
                            ld.append(ldi)
            opencj.write("]") # last line              
    return wite_fpath

def get_gsm_outputs(cjfn = "cjson", newfn = "msrap.cjson", 
    tempdname = "cjsontemp", jsonfiltpath = settings.gsmjsonfiltpath, 
    msrap_destpath = settings.gsmmsrapoutpath):
    """ Get the GSM-specific data from pipeline output files.

    Get the GSM-specific data from pipeline outputs run using composite JSON 
    files. Detects JSON composite files (cjfn) and output files (newfn) with
    matching timestamps. For matching files, the GSM IDs contained in the 
    composite JSON files are used to make the GSM-specific output files at the
    top level of msrap_destpath. Timestamps of new sample/GSM metadata files 
    match the valid detected composite file pairs.
    
    Arguments:
        * cjfn (str): File name string of composite JSON files containing the 
            GSM IDs.
        * newfn (str): File name string of the composite output metadata.
        * tempdname (str): Name of dir, at jsonfiltpath, to contain composite 
            JSON files.
        * jsonfiltpath (str): Path to filtered GSM JSON files.
        * msrap_destpath (str): Path to MetaSRA-pipeline output files.

    Returns:
        * NULL produces gsm metadata files at top level of msrap_destpath.

    """
    pathread_mdout = os.path.join(msrap_destpath, tempdname)
    if not os.path.exists(pathread_mdout):
        print("Path to composite metadata files doesn't exist: " + 
            str(pathread_mdout) + ". Returning...")
        return NULL
    re_cjson = re.compile("^" + cjfn + ".*")
    prl_cjson = os.listdir(pathread_mdout)
    prl_cjson = [fn for fn in prl_cjson if re_cjson.findall(fn)]
    re_msrap = re.compile("^" + newfn + ".*")
    prl_msrap = os.listdir(pathread_mdout)
    prl_msrap = [fn for fn in prl_msrap if re_msrap.findall(fn)]
    tsl_msrap = [fn.split(".")[2] for fn in prl_msrap];lx = []
    for fn in prl_cjson:
        if fn.split(".")[1] in tsl_msrap:
            print("Detected file match. Getting GSM IDs...")
            msrap_fmatch = [mn for mn in prl_msrap 
                                if mn.split(".")[2] == fn.split(".")[1]][0]
            ts = msrap_fmatch.split(".")[2] # use cfile ts
            rpath1 = os.path.join(pathread_mdout, fn)
            with open(rpath1, "r") as rf:
                for line in rf:
                    if line.split(":")[0] == '\t\t"gsm"':
                        gsmid = line.split(":")[1]
                        gsmid = gsmid.replace('"', '').replace('\n', '')
                        gsmid = gsmid.replace(",", "");llf.append(gsmid)
                        llf.append(gsmid)
            rpath2=os.path.join(pathread_mdout,msrap_fmatch);gsmct=0;ld=[]
            with open(rpath2, "r") as rf:
                for line in rf:
                    if not line in ["    }\n", "    },\n"]:
                        ld.append(line)
                    else:
                        gsmfname = ".".join(["msrapout", 
                            llf[gsmct].replace(",", ""), ts])
                        wpath = os.path.join(msrap_destpath, gsmfname)
                        print("Writing new file ", wpath)
                        with open(wpath, "w") as wf:
                            wf.write("[\n");wf.write("    {\n")
                            ld = ld[1::]
                            for line in ld:
                                wf.write(line)
                            wf.write("    }\n");wf.write("]");ld = []
                        print("Finished writing file num "+
                            str(gsmct)); gsmct += 1            
    return NULL

def run_msrap_compjson(json_flist=[], njint = 100, jsonpatt=".*json.filt$", 
    gsm_jsonpath = settings.gsmjsonfiltpath, tempdname = "cjsontemp",
    msrap_destpath = settings.gsmmsrapoutpath, newfnpattern = "msrap.cjson"):
    """ Run MetaSRA-pipeline on composite JSON files
        
        Runs the MetaSRA-pipeline on composite JSON files containing njint 
        samples' JSON-formatted metadata. The composite JSON files and the 
        composite metadata outputs are both written to tempfname at 
        msrap_destpath. After mapping, get_gsm_outputs() is called to make the
        GSM-specific files, which are output to the top level of msrap_destpath.
        
        Arguments:
            * json_flist (list, optional) : List of JSON filename(s) to process. 
                If not provided, automatically targets all JSON files at 
                gsm_jsondir.
            * njint (int): Number of JSON files per composite file to process.
            * jsonpatt (str): File name pattern for valid filtered JSON files.
            * gsm_jsonpath (str): Path to the filtered GSM JSON files directory.
            * tempdname (str): Dir, located at msrap_destpath, where composite
                JSON files and outputs are to be written.
            * msrap_destpath (str): Path where mapped metadata output files will
                be written.     
            * newfnpattern (str): File name pattern for mapped metadata output.
        
        Returns:
            * NULL, produces the composite file pairs and GSM metadata files.

    """
    eqfiltdict=get_queryfilt_dict()
    validgsmlist = [gsmid for gselist in list(eqfiltdict.values()) 
        for gsmid in gselist]
    msrap_runpath = settings.msraprunscriptpath
    msrap_oldgsm = os.listdir(msrap_runpath)
    msrap_oldgsm = [fn.split(".")[1] for fn in msrap_oldgsm]
    if not (json_flist and len(json_flist))>0:
        json_flist = os.listdir(gsm_jsonpath)
    print("Filtering GSM JSON filenames on pattern, existing msrap files...")
    gsm_json_fn_list = list(filter(re.compile(jsonpatt).match, json_flist))
    gsm_json_fn_list = [fn for fn in gsm_json_fn_list 
                        if not fn.split(.)[1] in msrap_oldgsm]
    cjsonpath = os.path.join(msrap_destpath, tempdname)
    os.makedirs(cjsonpath, exist_ok=True); msrap_statlist = []
    msrap_fn = settings.msrapfnstem; process_list = []
    rl = [r for r in range(0, len(gsm_json_fn_list), njint)]
    print("Running pipeline for composite JSON files...")
    for r in rl:
        ts = gettime_ntp() # use new ts for each new composite file pair
        jsonflist = gsm_json_fn_list[r:r+njint]
        cjreadpath = write_cjson(jffnv = jsonflist, jsonfiltpath = gsm_jsonpath, 
            msrap_destpath = msrap_destpath, ts = ts, tempdname = tempdname)
        newfn = ".".join([newfnpattern, ts])
        cjwritepath = os.path.join(cjsonpath, newfn)
        cmdlist = ['python2', msrap_runpath, "--fnvread", cjreadpath, 
            "--fnvwrite", cjwritepath]
        process_list.append(subprocess.call(cmdlist, shell=False))
        print("Finished index "+str(r))
    print("Extracting GSM data from composite JSON results...")
    get_gsm_outputs()
    return NULL

if __name__ == "__main__":
    run_msrap_compjson()
