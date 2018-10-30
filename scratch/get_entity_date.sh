#!/usr/bin/env bash
set -e

# functions
get_date() {
	# get entity last update date
	# arg1 : ftp address
	# arg2 : label
	newdate=$(curl -s -I $1 | awk '/Last-Modified/{ date=""; for(i=2;i<=NF;++i) date=(date " " $i); print date;}')
	returnlist=($2,$1,$newdate)
	echo ${returnlist[@]}
}

# script
get_date $1 $2