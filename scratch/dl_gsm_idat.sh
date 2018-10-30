#!/usr/bin/env bash
set -e

#------------
# FUNCTIONS
#------------

dlGSM () {
	# PURPOSE: download 1 gsm and write a logfile ('...SUCCESS' or 'NOTAVAILABLE')
	GSMI=$1 # grab the variable, a valid GSM ID, entered with the function
	ntlinei=$(grep "$GSMI" ./tempfiles/all_gse_newtally) # read ntfile, filtering on GSM
	GSEI=($(echo "$ntlinei" | grep -o "GSE[0-9]*")) # grab GSM's GSE ID from ntfile
	FTPI=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/${GSMI:0:${#GSMI}-3}nnn/$GSMI/suppl/ # format the FTP address properly (refer to GEOdb)
	IDAT_FILES_I=$(curl -i $FTPI | awk '{print $9}' | grep .*idat) # get idat filenames from path 
	NUM_IDATS=$(wc -w <<< "$IDAT_FILES_I") # tally idats found
	dlpath="./idats/"
	datpath="./datidat/"
	touch "$datpath$GSMI" # gsm datidat file, to coerce for mongo
	
	echo "num idats found: $NUM_IDATS" # return tally
	echo "idat files from curl: $IDAT_FILES_I" # return filenames found
	echo "starting the idat dl"
	# iterate on detected files
	if [[ "$NUM_IDATS" == 0 ]] ;
		then 
			echo "no idats found, writing file and exiting loop" # status message
			outfilenamei=$GSMI'_idat-transfer_NOIDATFILES'
			echo -n > $dlpath$outfilenamei ; # write status file, idats unavailable
		else 
			echo "num idats checked, initiating dl" # status message
			# get ftp list for dl
			IDAT_FTP_LIST=$(for idati in $IDAT_FILES_I; do echo $FTPI$idati; done) 
			ftparr=($IDAT_FTP_LIST)
			ESLIST="" # exist status list
			fileslist=($IDAT_FILES_I)
			idatct=0
			for idatftpi in $IDAT_FTP_LIST ; 
				do 
					datlist=()
					echo "trying dl for: "$idatftpi # status message
					# dl file, and store exit status
					cd $dlpath && { curl -R -O $idatftpi ; ESI=$(echo $?) ; cd -; } # 
					datlist+=($GSMI)
					datlist+=(${GSEI[0]})
					datlist+=(${fileslist[idatct]})
					datlist+=(${ftparr[idatct]})
					datlist+=($ESI)
					datlist+=($(date -r "$dlpath${fileslist[idatct]}" +'%Y-%m-%d_%H:%M:%S'))
					echo ${datlist[@]} >> "$datpath$GSMI" # add new line to datidat file
					ESLIST+="$ESI" 
					idatct=$((idatct+1));
			done
			echo "dl done, checking exit status obj and writing file"
			# write logfile for this dl
			if [ "$ESLIST"=="00" ] ;
				outfilenamei=$GSMI'_idat-transfer_SUCCESS'
				then echo $ESLIST > "$dlpath$outfilenamei"; 
				echo "Finished downloading data for sample "$1
			fi 
	fi	
}

#-------------
# BASH SCRIPT
#-------------
dlGSM $1
