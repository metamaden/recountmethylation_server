#!/usr/bin/env bash
set -e

# DESCRIPTION:
# Server-side helper functions to make and update recount-methylation database

# DEPENDENCIES:
# 1. Edirect utility, url: <https://www.ncbi.nlm.nih.gov/books/NBK179288/>

# NOTES:
# Run this script from 1 dir above 'GSE'


gsmNewTally () {
	# Purpose: run a fresh query on HM450 idat experiments,
	
	echo "Getting new E-direct query..."
	# write results of new query and use to make necessary subdir changes to ./GSM/..
	esearch -db gds -query "GPL13534[ACCN] AND idat[suppFile] AND gse[ETYP]" |   
	efetch -format docsum |   
	xtract -pattern DocumentSummary -element Id Accession > ./tempfiles/all_gse_newtally

	echo "Done with new E-direct query."
}

gsmMakeDirtree () {
	# Purpose: update the directory tree, make new dirs/subdirs as necessary
	echo "Making GSE dir tree..."
	newtally=$(cat ./tempfiles/all_gse_newtally) # read in the newtally file
	allGSE=$(grep -o "GSE[0-9]*" ./tempfiles/all_gse_newtally) # strip GSE IDs from ntfile
	allGSEarray=($(echo "$allGSE" | tr ',' '\n')) # format GSE IDs into array

	# iterate over GSEs, and check subdir tree for existence of GSM subdirs
	for gseid in "${allGSEarray[@]}" ;
		do  
			mkdir -p ./GSE/$gseid # make the experiment dir if nonexistant
			ntlinei=$(grep "$gseid" ./tempfiles/all_gse_newtally) # grab ntline for the GSE
			ntline_gsm=$(echo "$ntlinei" | grep -o "GSM[0-9]*") # grab GSMs for GSE
			ntline_gsm_array=($(echo "$ntline_gsm" | tr ',' '\n')) # format GSMs as array
			# iterate over the gsm IDs to make subdirs
			for gsmi in "${ntline_gsm_array[@]}" ;
				do mkdir -p ./GSE/$gseid/$gsmi # make GSM subdir if nonexistant
			done
	done

	echo "Done making GSE dir tree. Returning allGSEarray object to tempfiles."
	echo ${allGSEarray[@]} > ./tempfiles/allGSEarray # save allGSEarray for later use

}

gsmToDL () {
	# Note this downloads GSMs corresponding to missing GSE IDs (entire missing experiments)
	# Need to modify to accomodate cases where single samples/GSM IDs need updating
	echo "Checking for GSMs to download."
	gsmArray=() # init the array obj to be returned
	allGSE=$(grep -o "GSE[0-9]*" ./tempfiles/all_gse_newtally) # read allGSEarray file
	allGSEarray=($(echo "$allGSE" | tr ',' '\n')) # format allGSEarray file as array
	
	# check GSM subdirs for 'success' dl files, and add to list if nonexistant
	for GSEI in "${allGSEarray[@]}" ;
		do
			ntlinei=$(grep "$GSEI" ./tempfiles/all_gse_newtally) ; # grab ntfile expt mappings
			ntline_gsm_array=$(echo "$ntlinei" | grep -o "GSM[0-9]*") # filter ntfile line on GSM IDs
			ntline_gsm_array=($(echo "$ntline_gsm_array" | tr ',' '\n')) # format GSM IDs as array

			# check each GSM subdir for 'success' file
			echo "Checking GSM subdir's for success files..."
			echo "" > ./tempfiles/gsmArray # instantiate the gsmArray file 
			for GSMI in "${ntline_gsm_array[@]}" ;
				do 
					outfilenamei=$(echo $GSMI'_idat-transfer_SUCCESS') # 'success' filename expected

					# if 'success' file nonexistant, append GSM ID to list for new download
					if [ ! -e  './GSE/'$GSEI'/'$GSMI'/'$outfilenamei ] ; 
						then
							echo "Detected missing data file, appending info to gsmArray..."
							gsmArray=(${gsmArray[@]} $GSMI) ;
							echo $GSMI >> ./tempfiles/gsmArray # append new GSM id to end of growing array file
							echo "Current gsmArray length : "${#gsmArray[@]}
					fi
			done
	done
	echo "Finished checking for file download stats, returning file list to download to tempfiles/gsmArray..."
	#echo ${gsmArray[@]} > ./tempfiles/gsmArray # save gsmArray download list for later
}

#-------------
# BASH SCRIPT
#-------------
gsmNewTally
gsmMakeDirtree
gsmToDL
