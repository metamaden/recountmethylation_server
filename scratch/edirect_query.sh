#!/usr/bin/env bash
set -e

gse_query () {
	# Purpose: run a fresh query on HM450 idat experiments,
	echo "Getting new GSE e-direct query..."
	if [[ $# -eq 0 ]] ; then
		dldest='./gsequery'
	else
		dldest=$1
	fi
	esearch -db gds -query "GPL13534[ACCN] AND idat[suppFile] AND gse[ETYP]" |   
	efetch -format docsum |   
	xtract -pattern DocumentSummary -element Id Accession > $dldest
	echo "Query finished."
}

gsm_query () {
	# Purpose: run a fresh query on HM450 idat experiments,
	echo "Getting new GSM e-direct query..."
	if [[ $# -eq 0 ]] ; then
		dldest='./gsmquery'
	else
		dldest=$1
	fi
	esearch -db gds -query "GPL13534[ACCN] AND idat[suppFile] AND gsm[ETYP]" |   
	efetch -format docsum |   
	xtract -pattern DocumentSummary -element Id Accession > $dldest
	echo "Query finished."
}

edirect_query () {
	# Purpose conduct GSE and GSM queries
	gse_query $1
	gsm_query $1
}

edirect_query