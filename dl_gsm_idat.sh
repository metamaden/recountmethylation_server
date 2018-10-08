#!/usr/bin/env bash
set -e

#------------
# FUNCTIONS
#------------

dlGSM () {
	# PURPOSE: download 1 gsm and write a logfile ('...SUCCESS' or 'NOTAVAILABLE')

	GSMI=$1 # grab the variable, a valid GSM ID, entered with the function
	ntlinei=$(grep "$GSMI" ./tempfiles/all_gse_newtally) # read ntfile, filtering on GSM
	GSEI=$(echo "$ntlinei" | grep -o "GSE[0-9]*") # grab GSM's GSE ID from ntfile
	FTPI=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/${GSMI:0:${#GSMI}-3}nnn/$GSMI/suppl/ # format the FTP address properly (refer to GEOdb) 
	IDAT_FILES_I=$(curl -i $FTPI | awk '{print $9}' | grep .*idat) # grab the idat filenames for use in FTP address by stripping curl output from pointing at GEO containing dir
	echo "idat files from curl: "$IDAT_FILES_I # status message for troubleshooting
	NUM_IDATS=$(wc -w <<< "$IDAT_FILES_I") # count the idat files available
	echo "num idats found :"$NUM_IDATS # status message for troubleshooting
	dlpath="./GSE/$GSEI/$GSMI/" # form the dl dir path for GSM idat files

	echo "starting the idat dl" # status message for troubleshooting
	# loop over detected idat files, download them and write a status file
	if [[ "$NUM_IDATS" == 0 ]] ;
		then 
			echo "no idats found, writing file and exiting loop" # status message
			outfilenamei=$GSMI'_idat-transfer_NOIDATFILES'
			echo "" > $dlpath'/'$outfilenamei ; # write status file, idats unavailable
		else 
			echo "num idats checked, initiating dl" # status message
			IDAT_FTP_LIST=$(for idati in $IDAT_FILES_I; do echo $FTPI$idati; done) # form the final FTP addresses to dl idats
			ESLIST="" # init the exit status obj 
			# download the GSM idat files, at last
			for idatftpi in $IDAT_FTP_LIST ; 
				do 
					echo "trying dl for: "$idatftpi # status message
					cd $dlpath && { curl -O $idatftpi ; ESI=$(echo $?) ; cd -; } # traverse GSM subdirs to store idats correctly, storing exit status
					ESLIST+="$ESI" ;
			done

			echo "dl done, checking exit status obj and writing file" # status message
			# write 'success' outfile on expected exit status string '00'
			if [ "$ESLIST"=="00" ] ;
				outfilenamei=$GSMI'_idat-transfer_SUCCESS'
				then echo $ESLIST > $dlpath'/'$outfilenamei ; # write the exit status message and name with GSM id and 'SUCCESS'
				echo "Finished downloading data for sample "$1
			fi 
	fi	
}

#-------------
# BASH SCRIPT
#-------------
dlGSM $1
