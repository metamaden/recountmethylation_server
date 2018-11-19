import subprocess
import os

def msrap_prepare_json(gsm_json_filelist,gsm_json_dir='gsm_soft_json',
    dest_dir='msrap_torun',dest_filename='new_prepared_file'):
    """ Prepare GSM JSON metadata files for input to MetaSRA-pipeline
        Arguments
            * gsm_json_filelist : list of files to prepare
            * gsm_json_dir : directory of GSM JSON files to read
            * dest_dir : location to store new prepared file
            * dest_filename : prepared file name, to write
        Returns
            * null or error
    """
    os.makedirs(dest_dir, exist_ok=True)
    # manually add brackets '[' or ']', exclude from JSON concatenation
    with open(os.path.join(dest_dir,dest_filename),"w+") as preparedfile:
        preparedfile.write("[\n")
        for num, jsonfn in enumerate(gsm_json_filelist,0):
            with open(os.path.join(gsm_json_dir,jsonfn)) as jsonfile:
                for line in jsonfile:
                    if not ( line[0] == "]" or line[0] == "[" ):
                        if line == "  }\n" and num < (len(gsm_json_filelist)-1):
                            preparedfile.write("},\n")
                        else:
                            preparedfile.write(line)
        preparedfile.write("]")

def run_metasrapipeline(gsm_json_file,json_file_dir='msrap_torun',
    msrap_fn='new_msrap_outfile',
    msrap_dest_dir='msrap_output'):
    """ Run MetaSRA-pipeline on a GSM json file or list of JSON files
        Arguments
            * gsm_json_file : individual or group GSM JSON file
            * json_file_dir : dir of JSON file to map
            * msrap_fn : name of file to write new mapped MSRA-p. output
            * msrap_dest_dir : dest. dir. of final mapped output
        Returns
            * null or error, generating a new MetaSRA file as side effect 
    """
    os.makedirs(msrap_dest_dir, exist_ok=True)
    try:
        cmdlist = ['python',os.path.join('MetaSRA-pipeline','run_pipeline.py'),
        os.path.join(json_file_dir,gsm_json_file),'>',
        os.path.join(msrap_dest_dir,msrap_fn)]
        subprocess.call(cmdlist,shell=False)
    except subprocess.CalledProcessError as e:
        raise e

""" examples

gsm_json_filelist = ['GSM3177404.soft.json','GSM3177405.soft.json']
gsm_json_dir = 'gsm_soft_json'
l = []
for num, jsonfn in enumerate(gsm_json_filelist,0):
    with open(os.path.join(gsm_json_dir,jsonfn)) as jsonfile:
        for line in jsonfile:
            print(num)
            print(line)
            l.append(line)

msrap_prepare_json(['GSM3177404.soft.json','GSM3177405.soft.json'])

run_metasrapipeline('new_prepared_file')


"""