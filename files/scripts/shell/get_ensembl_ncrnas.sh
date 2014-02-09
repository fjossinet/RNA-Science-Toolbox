#!/bin/bash

DOWNLOAD_DIR=$1

if [ -d $DOWNLOAD_DIR ]
then
	echo "$DOWNLOAD_DIR already exists"
else
	mkdir -p $DOWNLOAD_DIR
	cd $DOWNLOAD_DIR
	wget -qr -nd -A "*.ncrna.fa.gz" "ftp://ftp.ensembl.org/pub/current_fasta/"
	gzip -d *.ncrna.fa.gz
fi

